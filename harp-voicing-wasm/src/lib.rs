pub mod abc;
pub mod chord_display;
pub mod chord_symbol;
pub mod harp;
pub mod music;
pub mod voicing;

use serde::Deserialize;
use wasm_bindgen::prelude::*;

#[derive(Deserialize)]
struct PrevVoicingJson {
    rh_midi: Vec<u8>,
    lh_midi: Vec<u8>,
}

#[wasm_bindgen]
pub struct HarpVoicer {
    key_root: u8,
}

#[wasm_bindgen]
impl HarpVoicer {
    #[wasm_bindgen(constructor)]
    pub fn new(key_root: u8) -> Self {
        HarpVoicer { key_root }
    }

    /// Set the key root (pitch class 0-11, where C=0).
    pub fn set_key_root(&mut self, key_root: u8) {
        self.key_root = key_root;
    }

    /// Get voicing suggestions as a JSON string.
    ///
    /// Returns an array of voicing objects, each with `rh`, `lh` (arrays of
    /// finger-note objects), and `score`.
    pub fn suggest_voicings(
        &self,
        melody_midi: u8,
        chord_symbol: &str,
        max_results: usize,
    ) -> String {
        let chord = match chord_symbol::parse_chord_symbol(chord_symbol) {
            Ok(c) => c,
            Err(e) => return format!(r#"{{"error":"{}"}}"#, e),
        };
        let voicings = voicing::enumerate_voicings(
            melody_midi,
            &chord,
            self.key_root,
            max_results,
        );
        serde_json::to_string(&voicings).unwrap_or_else(|e| format!(r#"{{"error":"{}"}}"#, e))
    }

    /// Identify a chord from MIDI notes, returning the ASCII chord symbol.
    pub fn identify_chord(&self, midi_notes: &[u8]) -> String {
        use std::collections::HashSet;

        let chord_templates: &[((&[u8], &str), (bool, &str), &[u8])] = &[
            ((&[0,4,7], ""),     (false, ""),   &[0,4,7]),
            ((&[0,3,7], "m"),    (true,  ""),   &[0,3,7]),
            ((&[0,4,7,10], "7"), (false, "7"),  &[0,4,7,10]),
            ((&[0,3,7,10], "m7"),(true,  "7"),  &[0,3,7,10]),
            ((&[0,4,7,11], "maj7"),(false,"M7"),&[0,4,7,11]),
            ((&[0,3,6], "dim"),  (true,  "-"),  &[0,3,6]),
            ((&[0,4,8], "aug"),  (false, "+"),  &[0,4,8]),
            ((&[0,5,7], "sus4"), (false, "~4"), &[0,5,7]),
            ((&[0,2,7], "sus2"), (false, "~2"), &[0,2,7]),
            ((&[0,3,6,9], "dim7"),(true, "-7"), &[0,3,6,9]),
            ((&[0,3,6,10],"m7b5"),(true, "%7"), &[0,3,6,10]),
        ];

        let pcs: HashSet<u8> = midi_notes.iter().map(|&m| music::pitch_class(m)).collect();
        if pcs.len() < 2 { return String::new(); }

        let roman_strs = ["I","II","III","IV","V","VI","VII"];

        // Dyad: just show the bass note's scale degree
        if pcs.len() == 2 {
            let bass_pc = midi_notes.iter().min().map(|&m| music::pitch_class(m)).unwrap_or(0);
            let deg = music::MAJOR_SCALE.iter().position(|&s| (self.key_root + s) % 12 == bass_pc);
            return match deg {
                Some(d) => roman_strs[d].to_string(),
                None => String::new(),
            };
        }

        let bass_pc = midi_notes.iter().min().map(|&m| music::pitch_class(m));
        let mut best_root = 0u8;
        let mut best_idx = 0usize;
        let mut best_score = -999i32;

        for root in 0..12u8 {
            let intervals: HashSet<u8> = pcs.iter().map(|&pc| (pc + 12 - root) % 12).collect();
            for (idx, &((tmpl, _), _, _)) in chord_templates.iter().enumerate() {
                let tmpl_set: HashSet<u8> = tmpl.iter().copied().collect();
                let matched = intervals.intersection(&tmpl_set).count() as i32;
                let extra = intervals.difference(&tmpl_set).count() as i32;
                let missing = tmpl_set.difference(&intervals).count() as i32;
                let mut score = matched * 3 - extra * 2 - missing;
                if bass_pc == Some(root) { score += 4; }
                if tmpl.len() == 3 { score += 2; }
                score -= idx as i32 / 10;
                if score > best_score {
                    best_score = score;
                    best_root = root;
                    best_idx = idx;
                }
            }
        }

        let degree = music::MAJOR_SCALE.iter().position(|&s| (self.key_root + s) % 12 == best_root);
        let deg = match degree { Some(d) => d, None => return String::new() };

        let (_, (lower, suffix), inv_order) = &chord_templates[best_idx];
        let roman = if *lower { roman_strs[deg].to_lowercase() } else { roman_strs[deg].to_string() };

        let inv = match bass_pc {
            Some(bpc) => {
                let bi = (bpc + 12 - best_root) % 12;
                inv_order.iter().position(|&x| x == bi).unwrap_or(0)
            }
            None => 0,
        };

        let mut result = format!("{}{}", roman, suffix);
        if inv > 0 { result.push_str(&format!("^{}", inv)); }
        result
    }

    /// Snap a chromatic MIDI note to the nearest diatonic string in the current key.
    pub fn snap_to_diatonic(&self, midi: u8) -> u8 {
        music::snap_to_diatonic(midi, self.key_root)
    }

    /// Voice from SATB parts, returning a JSON voicing object.
    /// alto, tenor, bass are optional (pass 255 to skip).
    pub fn voice_from_satb(
        &self,
        soprano: u8,
        alto: u8,
        tenor: u8,
        bass: u8,
        prev_json: &str,
    ) -> String {
        let alto_opt = if alto == 255 { None } else { Some(alto) };
        let tenor_opt = if tenor == 255 { None } else { Some(tenor) };
        let bass_opt = if bass == 255 { None } else { Some(bass) };

        let prev = if prev_json.is_empty() {
            None
        } else {
            serde_json::from_str::<PrevVoicingJson>(prev_json).ok().map(|p| {
                voicing::PrevVoicing {
                    rh_midi: p.rh_midi,
                    lh_midi: p.lh_midi,
                }
            })
        };

        let identify = |notes: &[u8], kr: u8| -> String {
            let voicer = HarpVoicer { key_root: kr };
            voicer.identify_chord(notes)
        };

        match voicing::voice_from_satb(
            soprano, alto_opt, tenor_opt, bass_opt,
            self.key_root,
            prev.as_ref(),
            &identify,
        ) {
            Some(v) => serde_json::to_string(&v).unwrap_or_else(|e| format!(r#"{{"error":"{}"}}"#, e)),
            None => "null".to_string(),
        }
    }

    /// Process an ABC file (SATB or lead sheet) and return JSON hymns with voicings.
    pub fn process_abc(&self, abc_text: &str) -> String {
        let hymns = abc::parse_lead_sheets(abc_text);
        serde_json::to_string(&hymns)
            .unwrap_or_else(|e| format!(r#"{{"error":"{}"}}"#, e))
    }

    /// Convert ASCII chord symbol to traditional notation.
    pub fn chord_to_traditional(&self, ascii: &str) -> String {
        chord_display::chord_to_traditional(ascii)
    }

    /// Convert ASCII chord symbol to spoken form.
    pub fn chord_to_say(&self, ascii: &str) -> String {
        chord_display::chord_to_say(ascii)
    }

    /// Show the pitch classes for a chord in the current key.
    pub fn chord_to_example(&self, ascii: &str) -> String {
        chord_display::chord_to_example(ascii, self.key_root)
    }

    /// Show the interval structure of a chord.
    pub fn chord_to_intervals(&self, ascii: &str) -> String {
        chord_display::chord_to_intervals(ascii)
    }

    /// Parse a chord symbol and return its structure as JSON.
    pub fn parse_chord(&self, chord_symbol: &str) -> String {
        match chord_symbol::parse_chord_symbol(chord_symbol) {
            Ok(c) => {
                let root_pc = music::degree_to_pc(c.root_degree, self.key_root);
                let tones: Vec<u8> = c.intervals.iter().map(|&iv| (root_pc + iv) % 12).collect();
                let tone_names: Vec<&str> = tones
                    .iter()
                    .map(|&pc| music::PC_NAMES[pc as usize])
                    .collect();
                format!(
                    r#"{{"root_degree":{},"root_pc":{},"intervals":{:?},"tones":{:?},"tone_names":{:?},"inversion":{}}}"#,
                    c.root_degree, root_pc, c.intervals, tones, tone_names, c.inversion
                )
            }
            Err(e) => format!(r#"{{"error":"{}"}}"#, e),
        }
    }
}
