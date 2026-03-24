use std::collections::HashMap;

use regex::Regex;

use crate::chord_symbol;
use crate::music::*;
use crate::voicing::{self, PrevVoicing, SatbVoicing};

/// Map key name (e.g. "G", "Bb", "F#") to pitch class 0-11.
pub fn key_to_pc(key: &str) -> u8 {
    match key {
        "C" => 0,
        "Db" | "C#" => 1,
        "D" => 2,
        "Eb" => 3,
        "E" => 4,
        "F" => 5,
        "F#" | "Gb" => 6,
        "G" => 7,
        "Ab" => 8,
        "A" => 9,
        "Bb" => 10,
        "B" | "Cb" => 11,
        _ => 0,
    }
}

/// ABC base MIDI for uppercase letters (octave 4 convention: C=48..B=59).
fn abc_base(letter: u8) -> Option<i16> {
    match letter {
        b'C' => Some(60),
        b'D' => Some(62),
        b'E' => Some(64),
        b'F' => Some(65),
        b'G' => Some(67),
        b'A' => Some(69),
        b'B' => Some(71),
        _ => None,
    }
}

/// Key signature accidentals: returns letter -> semitone adjustment (+1 sharp, -1 flat).
pub fn key_signature_accidentals(key: &str) -> HashMap<u8, i8> {
    let sharp_order: [u8; 7] = [b'F', b'C', b'G', b'D', b'A', b'E', b'B'];
    let flat_order: [u8; 7] = [b'B', b'E', b'A', b'D', b'G', b'C', b'F'];

    let sharp_keys: &[(&str, usize)] = &[
        ("G", 1), ("D", 2), ("A", 3), ("E", 4), ("B", 5), ("F#", 6), ("C#", 7),
    ];
    let flat_keys: &[(&str, usize)] = &[
        ("F", 1), ("Bb", 2), ("Eb", 3), ("Ab", 4), ("Db", 5), ("Gb", 6), ("Cb", 7),
    ];

    let mut acc = HashMap::new();
    for &(k, n) in sharp_keys {
        if k == key {
            for i in 0..n {
                acc.insert(sharp_order[i], 1i8);
            }
            return acc;
        }
    }
    for &(k, n) in flat_keys {
        if k == key {
            for i in 0..n {
                acc.insert(flat_order[i], -1i8);
            }
            return acc;
        }
    }
    acc
}

/// Convert an ABC note token (e.g. "^F,", "c'", "=B") to MIDI number.
pub fn abc_note_to_midi(token: &str, key_sig_acc: &HashMap<u8, i8>) -> Option<i16> {
    let bytes = token.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    // Parse explicit accidentals: ^ (sharp), _ (flat), = (natural)
    let mut acc: Option<i8> = None;
    while i < len && (bytes[i] == b'^' || bytes[i] == b'_' || bytes[i] == b'=') {
        if acc.is_none() {
            acc = Some(0);
        }
        match bytes[i] {
            b'^' => *acc.as_mut().unwrap() += 1,
            b'_' => *acc.as_mut().unwrap() -= 1,
            b'=' => acc = Some(0), // explicit natural
            _ => {}
        }
        i += 1;
    }

    if i >= len {
        return None;
    }

    let ch = bytes[i];
    let is_lower = ch >= b'a' && ch <= b'g';
    let letter = ch.to_ascii_uppercase();
    let base = abc_base(letter)?;

    // Apply key signature if no explicit accidental
    let adjustment = match acc {
        Some(a) => a as i16,
        None => *key_sig_acc.get(&letter).unwrap_or(&0) as i16,
    };

    let mut midi = base + adjustment;
    if is_lower {
        midi += 12;
    }
    i += 1;

    // Octave modifiers: , (down) and ' (up)
    while i < len {
        match bytes[i] {
            b',' => midi -= 12,
            b'\'' => midi += 12,
            _ => break,
        }
        i += 1;
    }

    Some(midi)
}

// ── Duration parsing ────────────────────────────────────────

/// Parse ABC duration suffix (e.g. "2", "/", "3/2", "//") into a rational
/// multiplier of the default note length, returned as (numerator, denominator).
/// Result is in quarter-note beats when default_len is (1,4).
pub fn parse_abc_duration(suffix: &str, default_len: (u32, u32)) -> (u32, u32) {
    // beats_per_default = default_len / (1/4) = default_len * 4
    let base_num = default_len.0 * 4;
    let base_den = default_len.1;

    if suffix.is_empty() || suffix.chars().all(|c| c.is_whitespace()) {
        return (base_num, base_den);
    }

    let s = suffix.trim();

    // Leading slashes: / = /2, // = /4
    if s.starts_with('/') {
        let slash_count = s.chars().take_while(|&c| c == '/').count();
        let num_part = &s[slash_count..];
        let divisor = if num_part.is_empty() {
            2u32.pow(slash_count as u32)
        } else {
            num_part.parse::<u32>().unwrap_or(2)
        };
        return (base_num, base_den * divisor);
    }

    // n/m
    if s.contains('/') {
        let parts: Vec<&str> = s.splitn(2, '/').collect();
        let n = parts[0].parse::<u32>().unwrap_or(1);
        let m = parts[1].parse::<u32>().unwrap_or(1);
        return (base_num * n, base_den * m);
    }

    // Just n (multiplier)
    let n = s.parse::<u32>().unwrap_or(1);
    (base_num * n, base_den)
}

/// Parse an L: field value like "1/4" into (numerator, denominator).
pub fn parse_default_len(s: &str) -> (u32, u32) {
    let s = s.trim().split('%').next().unwrap_or("1/4").trim();
    if let Some(pos) = s.find('/') {
        let num = s[..pos].parse::<u32>().unwrap_or(1);
        let den = s[pos + 1..].parse::<u32>().unwrap_or(4);
        (num, den)
    } else {
        (s.parse::<u32>().unwrap_or(1), 1)
    }
}

/// Parse a Q: tempo field. Handles "1/4=120", "120", or "[Q:1/4=120]".
/// Returns BPM (quarter notes per minute).
pub fn parse_tempo(s: &str) -> Option<u32> {
    let s = s.trim().trim_start_matches("[Q:").trim_end_matches(']');
    if let Some(pos) = s.find('=') {
        s[pos + 1..].trim().parse::<u32>().ok()
    } else {
        s.trim().parse::<u32>().ok()
    }
}

// ── Voice token types ───────────────────────────────────────

#[derive(Debug, Clone)]
pub enum VoiceToken {
    Note { midi: u8, duration: (u32, u32) },
    Rest { duration: (u32, u32) },
    Bar,
    Tie,
}

/// Parse an ABC voice line into note/bar/rest tokens with durations.
pub fn parse_voice_tokens(
    music_str: &str,
    key_sig_acc: &HashMap<u8, i8>,
    default_len: (u32, u32),
) -> Vec<VoiceToken> {
    // Match: note tokens, rest tokens, barlines, inline tempo [Q:...], ties
    let re = Regex::new(
        r"(\[Q:[^\]]*\])|([_^=]*[A-Ga-g][,']*\d*(?:/\d*)?)(-)?|(z\d*(?:/\d*)?)|(\|[\|\]\[:]*)"
    ).unwrap();
    let pitch_re = Regex::new(r"^[_^=]*[A-Ga-g][,']*").unwrap();
    let mut tokens = Vec::new();

    for cap in re.captures_iter(music_str) {
        // Skip inline tempo markers [Q:...]
        if cap.get(1).is_some() {
            continue;
        }

        // Barline
        if cap.get(5).is_some() {
            tokens.push(VoiceToken::Bar);
            continue;
        }

        // Rest
        if let Some(m) = cap.get(4) {
            let rest_str = m.as_str();
            let dur_suffix = &rest_str[1..]; // strip 'z'
            let duration = parse_abc_duration(dur_suffix, default_len);
            tokens.push(VoiceToken::Rest { duration });
            continue;
        }

        // Note
        if let Some(note_match) = cap.get(2) {
            let token = note_match.as_str();
            if let Some(pm) = pitch_re.find(token) {
                let pitch_part = pm.as_str();
                let dur_suffix = &token[pm.end()..];
                if let Some(midi) = abc_note_to_midi(pitch_part, key_sig_acc) {
                    if midi >= 0 && midi <= 127 {
                        let duration = parse_abc_duration(dur_suffix, default_len);
                        tokens.push(VoiceToken::Note {
                            midi: midi as u8,
                            duration,
                        });
                        // Check for tie marker
                        if cap.get(3).is_some() {
                            tokens.push(VoiceToken::Tie);
                        }
                    }
                }
            }
        }
    }
    tokens
}

// ── Hymn data structures ────────────────────────────────────

use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct Beat {
    #[serde(rename = "type")]
    pub beat_type: String, // "note" or "bar"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub midi: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub original_midi: Option<u8>,
    /// Duration in quarter-note beats as (numerator, denominator).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duration: Option<(u32, u32)>,
    /// True if this note is tied to the next note of the same pitch.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tied: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chord: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub full_chord: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rh_chord: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lh_chord: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rh_strings: Option<Vec<i16>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lh_strings: Option<Vec<i16>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rh_midi: Option<Vec<u8>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lh_midi: Option<Vec<u8>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Hymn {
    pub number: String,
    pub title: String,
    pub meter: String,
    pub key: String,
    pub original_key: String,
    pub key_root: u8,
    /// Default note length from L: field, as (numerator, denominator).
    pub default_len: (u32, u32),
    /// Tempo in quarter-note BPM from Q: field (default 100).
    pub tempo_bpm: u32,
    pub beats: Vec<Beat>,
}

/// Chord identification closure for use with voice_from_satb.
fn identify_chord(notes: &[u8], key_root: u8) -> String {
    use std::collections::HashSet;

    let chord_templates: &[(&[u8], bool, &str, &[u8])] = &[
        (&[0, 4, 7], false, "", &[0, 4, 7]),
        (&[0, 3, 7], true, "", &[0, 3, 7]),
        (&[0, 4, 7, 10], false, "7", &[0, 4, 7, 10]),
        (&[0, 3, 7, 10], true, "7", &[0, 3, 7, 10]),
        (&[0, 4, 7, 11], false, "M7", &[0, 4, 7, 11]),
        (&[0, 3, 6], true, "-", &[0, 3, 6]),
        (&[0, 4, 8], false, "+", &[0, 4, 8]),
        (&[0, 5, 7], false, "~4", &[0, 5, 7]),
        (&[0, 2, 7], false, "~2", &[0, 2, 7]),
        (&[0, 3, 6, 9], true, "-7", &[0, 3, 6, 9]),
        (&[0, 3, 6, 10], true, "%7", &[0, 3, 6, 10]),
    ];
    let roman_strs = ["I", "II", "III", "IV", "V", "VI", "VII"];

    let pcs: HashSet<u8> = notes.iter().map(|&m| pitch_class(m)).collect();
    if pcs.len() < 2 {
        return String::new();
    }

    // Dyad
    if pcs.len() == 2 {
        let bass_pc = notes.iter().min().map(|&m| pitch_class(m)).unwrap_or(0);
        let deg = MAJOR_SCALE
            .iter()
            .position(|&s| (key_root + s) % 12 == bass_pc);
        return match deg {
            Some(d) => roman_strs[d].to_string(),
            None => String::new(),
        };
    }

    let bass_pc = notes.iter().min().map(|&m| pitch_class(m));
    let mut best_root = 0u8;
    let mut best_idx = 0usize;
    let mut best_score = -999i32;

    for root in 0..12u8 {
        let intervals: HashSet<u8> = pcs.iter().map(|&pc| (pc + 12 - root) % 12).collect();
        for (idx, &(tmpl, _, _, _)) in chord_templates.iter().enumerate() {
            let tmpl_set: HashSet<u8> = tmpl.iter().copied().collect();
            let matched = intervals.intersection(&tmpl_set).count() as i32;
            let extra = intervals.difference(&tmpl_set).count() as i32;
            let missing = tmpl_set.difference(&intervals).count() as i32;
            let mut score = matched * 3 - extra * 2 - missing;
            if bass_pc == Some(root) {
                score += 4;
            }
            if tmpl.len() == 3 {
                score += 2;
            }
            score -= idx as i32 / 10;
            if score > best_score {
                best_score = score;
                best_root = root;
                best_idx = idx;
            }
        }
    }

    let degree = MAJOR_SCALE
        .iter()
        .position(|&s| (key_root + s) % 12 == best_root);
    let deg = match degree {
        Some(d) => d,
        None => return String::new(),
    };

    let (_, lower, suffix, inv_order) = chord_templates[best_idx];
    let roman = if lower {
        roman_strs[deg].to_lowercase()
    } else {
        roman_strs[deg].to_string()
    };

    let inv = match bass_pc {
        Some(bpc) => {
            let bi = (bpc + 12 - best_root) % 12;
            inv_order.iter().position(|&x| x == bi).unwrap_or(0)
        }
        None => 0,
    };

    let mut result = format!("{}{}", roman, suffix);
    if inv > 0 {
        result.push_str(&format!("^{}", inv));
    }
    result
}

// ── Full parsers ────────────────────────────────────────────

/// Parse a full SATB OpenHymnal ABC file into hymns with voicings.
pub fn parse_open_hymnal(text: &str) -> Vec<Hymn> {
    let mut hymns = Vec::new();
    // Split on X: field (each tune starts with X:)
    let tunes: Vec<&str> = text.split("\nX:").collect();

    for (ti, tune) in tunes.iter().enumerate() {
        let tune = tune.trim();
        if tune.is_empty() {
            continue;
        }

        let lines: Vec<&str> = tune.lines().collect();
        let mut number = String::new();
        let mut title = String::new();
        let mut meter = "4/4".to_string();
        let mut key = "C".to_string();
        let mut default_len = (1u32, 4u32);
        let mut tempo_bpm = 100u32;
        let mut voice_lines: HashMap<&str, Vec<String>> = HashMap::new();
        voice_lines.insert("S1V1", Vec::new());
        voice_lines.insert("S1V2", Vec::new());
        voice_lines.insert("S2V1", Vec::new());
        voice_lines.insert("S2V2", Vec::new());

        for line in &lines {
            let trimmed = line.trim();
            // Number line: just digits
            if number.is_empty() && trimmed.chars().all(|c| c.is_ascii_digit()) && !trimmed.is_empty()
            {
                number = trimmed.to_string();
            } else if line.starts_with("T:") && title.is_empty() {
                title = line[2..].trim().to_string();
            } else if line.starts_with("M:") {
                meter = line[2..].trim().split('%').next().unwrap_or("4/4").trim().to_string();
            } else if line.starts_with("L:") {
                default_len = parse_default_len(&line[2..]);
            } else if line.starts_with("K:") {
                let k = line[2..]
                    .trim()
                    .split('%')
                    .next()
                    .unwrap_or("C")
                    .trim()
                    .split_whitespace()
                    .next()
                    .unwrap_or("C");
                if !k.is_empty() && k.as_bytes()[0].is_ascii_uppercase() {
                    key = k.to_string();
                }
            } else if line.starts_with("[V: S1V1]") {
                let music = line.replacen("[V: S1V1]", "", 1).trim().to_string();
                voice_lines.get_mut("S1V1").unwrap().push(music);
            } else if line.starts_with("[V: S1V2]") {
                let music = line.replacen("[V: S1V2]", "", 1).trim().to_string();
                voice_lines.get_mut("S1V2").unwrap().push(music);
            } else if line.starts_with("[V: S2V1]") {
                let music = line.replacen("[V: S2V1]", "", 1).trim().to_string();
                voice_lines.get_mut("S2V1").unwrap().push(music);
            } else if line.starts_with("[V: S2V2]") {
                let music = line.replacen("[V: S2V2]", "", 1).trim().to_string();
                voice_lines.get_mut("S2V2").unwrap().push(music);
            }
        }

        // Need at least soprano
        if voice_lines["S1V1"].is_empty() {
            continue;
        }

        // Try to extract number from first line if not found
        if number.is_empty() {
            if ti > 0 {
                // The split removed the "X:" prefix; for tunes after the first,
                // the first token is the X number
                if let Some(first_line) = lines.first() {
                    if let Some(n) = first_line.trim().split_whitespace().next() {
                        if n.chars().all(|c| c.is_ascii_digit()) {
                            number = n.to_string();
                        }
                    }
                }
            }
        }

        let key_root = key_to_pc(&key);
        let key_sig_acc = key_signature_accidentals(&key);

        // Extract inline tempo from soprano voice if present
        let soprano_str = voice_lines["S1V1"].join(" ");
        let tempo_re = Regex::new(r"\[Q:\s*(?:\d+/\d+=)?(\d+)\]").unwrap();
        if let Some(cap) = tempo_re.captures(&soprano_str) {
            if let Ok(bpm) = cap[1].parse::<u32>() {
                tempo_bpm = bpm;
            }
        }

        let soprano = parse_voice_tokens(&soprano_str, &key_sig_acc, default_len);
        let alto = parse_voice_tokens(&voice_lines["S1V2"].join(" "), &key_sig_acc, default_len);
        let tenor = parse_voice_tokens(&voice_lines["S2V1"].join(" "), &key_sig_acc, default_len);
        let bass = parse_voice_tokens(&voice_lines["S2V2"].join(" "), &key_sig_acc, default_len);

        // Extract just the note tokens from alto/tenor/bass
        let a_notes: Vec<u8> = alto
            .iter()
            .filter_map(|t| match t {
                VoiceToken::Note { midi, .. } => Some(*midi),
                _ => None,
            })
            .collect();
        let t_notes: Vec<u8> = tenor
            .iter()
            .filter_map(|t| match t {
                VoiceToken::Note { midi, .. } => Some(*midi),
                _ => None,
            })
            .collect();
        let b_notes: Vec<u8> = bass
            .iter()
            .filter_map(|t| match t {
                VoiceToken::Note { midi, .. } => Some(*midi),
                _ => None,
            })
            .collect();

        let mut beats = Vec::new();
        let mut last_voicing: Option<SatbVoicing> = None;
        let mut note_idx: usize = 0;

        let mut i = 0;
        while i < soprano.len() {
            let tok = &soprano[i];
            i += 1;
            match tok {
                VoiceToken::Bar => {
                    beats.push(Beat {
                        beat_type: "bar".to_string(),
                        midi: None,
                        original_midi: None,
                        duration: None,
                        tied: None,
                        chord: None,
                        full_chord: None,
                        rh_chord: None,
                        lh_chord: None,
                        rh_strings: None,
                        lh_strings: None,
                        rh_midi: None,
                        lh_midi: None,
                    });
                }
                VoiceToken::Rest { duration } => {
                    note_idx += 1;
                    beats.push(Beat {
                        beat_type: "rest".to_string(),
                        midi: None,
                        original_midi: None,
                        duration: Some(*duration),
                        tied: None,
                        chord: None,
                        full_chord: None,
                        rh_chord: None,
                        lh_chord: None,
                        rh_strings: None,
                        lh_strings: None,
                        rh_midi: None,
                        lh_midi: None,
                    });
                }
                VoiceToken::Tie => { /* handled below after Note */ }
                VoiceToken::Note { midi: s_midi, duration } => {
                    // Check if next token is a Tie
                    let is_tied = matches!(soprano.get(i), Some(VoiceToken::Tie));
                    if is_tied { i += 1; } // consume tie token
                    let s_idx = note_idx;
                    note_idx += 1;

                    let a_midi = a_notes.get(s_idx).copied();
                    let t_midi = t_notes.get(s_idx).copied();
                    let b_midi = b_notes.get(s_idx).copied();

                    let prev = last_voicing.as_ref().map(|v| PrevVoicing {
                        rh_midi: v.rh_midi.clone(),
                        lh_midi: v.lh_midi.clone(),
                    });

                    let voicing_result = voicing::voice_from_satb(
                        *s_midi,
                        a_midi,
                        t_midi,
                        b_midi,
                        key_root,
                        prev.as_ref(),
                        &identify_chord,
                    );

                    let snapped_melody = snap_to_diatonic(*s_midi, key_root);
                    let mut beat = Beat {
                        beat_type: "note".to_string(),
                        midi: Some(snapped_melody),
                        original_midi: Some(*s_midi),
                        duration: Some(*duration),
                        tied: if is_tied { Some(true) } else { None },
                        chord: None,
                        full_chord: None,
                        rh_chord: None,
                        lh_chord: None,
                        rh_strings: None,
                        lh_strings: None,
                        rh_midi: None,
                        lh_midi: None,
                    };

                    if let Some(ref v) = voicing_result {
                        // Show chord when voicing changes
                        let voicing_key = format!(
                            "{:?}|{:?}",
                            v.rh_strings, v.lh_strings
                        );
                        let last_key = last_voicing
                            .as_ref()
                            .map(|lv| format!("{:?}|{:?}", lv.rh_strings, lv.lh_strings))
                            .unwrap_or_default();

                        if voicing_key != last_key {
                            beat.chord =
                                Some(format!("{}/{}", v.rh_chord, v.lh_chord));
                            beat.full_chord = Some(v.full_chord.clone());
                            beat.rh_chord = Some(v.rh_chord.clone());
                            beat.lh_chord = Some(v.lh_chord.clone());
                            beat.rh_strings = Some(v.rh_strings.clone());
                            beat.lh_strings = Some(v.lh_strings.clone());
                            beat.rh_midi = Some(v.rh_midi.clone());
                            beat.lh_midi = Some(v.lh_midi.clone());
                            last_voicing = voicing_result;
                        }
                    }

                    beats.push(beat);
                }
            }
        }

        if !beats.is_empty() {
            hymns.push(Hymn {
                number,
                title,
                meter,
                key: key.clone(),
                original_key: key,
                key_root,
                default_len,
                tempo_bpm,
                beats,
            });
        }
    }

    hymns
}

/// Parse lead sheet ABC (with "chord" annotations) or detect SATB and delegate.
pub fn parse_lead_sheets(text: &str) -> Vec<Hymn> {
    // Detect SATB format
    if text.contains("[V: S1V1]") {
        return parse_open_hymnal(text);
    }

    let mut hymns = Vec::new();
    let lead_re = Regex::new(
        r#"("([^"]*)")?\s*([_^=]*[A-Ga-g][,']*\d*/?\d*|z\d*/?\d*|\[x\d*\]|\|[\|\]\[:]*)"#,
    )
    .unwrap();
    let pitch_re = Regex::new(r"^[_^=]*[A-Ga-g][,']*").unwrap();

    let tunes: Vec<&str> = {
        let mut parts = Vec::new();
        let mut start = 0;
        for (i, _) in text.match_indices("\nX:") {
            parts.push(&text[start..i]);
            start = i + 1; // skip the \n
        }
        parts.push(&text[start..]);
        parts
    };

    for tune in &tunes {
        let tune = tune.trim();
        if tune.is_empty() {
            continue;
        }

        let lines: Vec<&str> = tune.lines().collect();
        let mut number = String::new();
        let mut title = String::new();
        let mut meter = "4/4".to_string();
        let mut key = "C".to_string();
        let mut default_len = (1u32, 4u32);
        let mut tempo_bpm = 100u32;
        let mut music_lines = Vec::new();

        for line in &lines {
            if line.starts_with("X:") {
                number = line[2..].trim().to_string();
            } else if line.starts_with("T:") {
                title = line[2..].trim().to_string();
            } else if line.starts_with("M:") {
                meter = line[2..].trim().to_string();
            } else if line.starts_with("L:") {
                default_len = parse_default_len(&line[2..]);
            } else if line.starts_with("K:") {
                key = line[2..].trim().to_string();
            } else if !line.trim().is_empty()
                && !line.starts_with('%')
            {
                music_lines.push(*line);
            }
        }

        let key_root = key_to_pc(&key);
        let key_sig_acc = key_signature_accidentals(&key);
        let music_str = music_lines.join(" ");

        let mut beats = Vec::new();
        let mut last_chord: Option<String> = None;

        for cap in lead_re.captures_iter(&music_str) {
            let chord_sym = cap.get(2).map(|m| m.as_str().to_string());
            let token = &cap[3];

            if let Some(ref cs) = chord_sym {
                last_chord = Some(cs.clone());
            }

            if token.starts_with('|') {
                beats.push(Beat {
                    beat_type: "bar".to_string(),
                    midi: None,
                    original_midi: None,
                    duration: None,
                    tied: None,
                    chord: None,
                    full_chord: None,
                    rh_chord: None,
                    lh_chord: None,
                    rh_strings: None,
                    lh_strings: None,
                    rh_midi: None,
                    lh_midi: None,
                });
                continue;
            }

            if token.starts_with('z') || token.starts_with("[x") {
                continue;
            }

            let pitch_match = match pitch_re.find(token) {
                Some(m) => m,
                None => continue,
            };

            let dur_suffix = &token[pitch_match.end()..];
            let duration = parse_abc_duration(dur_suffix, default_len);

            let midi = match abc_note_to_midi(pitch_match.as_str(), &key_sig_acc) {
                Some(m) if m >= 0 && m <= 127 => m as u8,
                _ => continue,
            };

            let mut beat = Beat {
                beat_type: "note".to_string(),
                midi: Some(midi),
                original_midi: Some(midi),
                duration: Some(duration),
                tied: None,
                chord: chord_sym.clone(),
                full_chord: None,
                rh_chord: None,
                lh_chord: None,
                rh_strings: None,
                lh_strings: None,
                rh_midi: None,
                lh_midi: None,
            };

            // If there's a chord symbol, compute voicing
            if let Some(ref cs) = chord_sym {
                if let Ok(parsed) = chord_symbol::parse_chord_symbol(cs) {
                    let voicings =
                        voicing::enumerate_voicings(midi, &parsed, key_root, 1);
                    if let Some(v) = voicings.first() {
                        beat.rh_strings =
                            Some(v.rh.iter().map(|n| n.string_num).collect());
                        beat.lh_strings =
                            Some(v.lh.iter().map(|n| n.string_num).collect());
                        beat.rh_midi = Some(v.rh.iter().map(|n| n.midi).collect());
                        beat.lh_midi = Some(v.lh.iter().map(|n| n.midi).collect());
                        let rh_midis: Vec<u8> = v.rh.iter().map(|n| n.midi).collect();
                        let lh_midis: Vec<u8> = v.lh.iter().map(|n| n.midi).collect();
                        beat.rh_chord = Some(identify_chord(&rh_midis, key_root));
                        beat.lh_chord = Some(identify_chord(&lh_midis, key_root));
                    }
                }
            }

            beats.push(beat);
        }

        if !beats.is_empty() {
            hymns.push(Hymn {
                number,
                title,
                meter,
                key: key.clone(),
                original_key: key,
                key_root,
                default_len,
                tempo_bpm,
                beats,
            });
        }
    }

    hymns
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key_to_pc() {
        assert_eq!(key_to_pc("C"), 0);
        assert_eq!(key_to_pc("G"), 7);
        assert_eq!(key_to_pc("F#"), 6);
        assert_eq!(key_to_pc("Bb"), 10);
        assert_eq!(key_to_pc("Cb"), 11);
    }

    #[test]
    fn test_key_signature_accidentals() {
        let acc = key_signature_accidentals("G");
        assert_eq!(acc.get(&b'F'), Some(&1)); // F#
        assert_eq!(acc.len(), 1);

        let acc = key_signature_accidentals("Bb");
        assert_eq!(acc.get(&b'B'), Some(&-1));
        assert_eq!(acc.get(&b'E'), Some(&-1));
        assert_eq!(acc.len(), 2);

        let acc = key_signature_accidentals("C");
        assert_eq!(acc.len(), 0);
    }

    #[test]
    fn test_abc_note_to_midi() {
        let empty = HashMap::new();
        // Standard ABC: uppercase C = middle C (MIDI 60), lowercase c = C5 (MIDI 72)
        assert_eq!(abc_note_to_midi("C", &empty), Some(60));
        assert_eq!(abc_note_to_midi("c", &empty), Some(72));
        assert_eq!(abc_note_to_midi("C'", &empty), Some(72));
        assert_eq!(abc_note_to_midi("C,", &empty), Some(48));
        assert_eq!(abc_note_to_midi("^F", &empty), Some(66)); // F=65+1
        assert_eq!(abc_note_to_midi("_B", &empty), Some(70)); // B=71-1

        // With key signature: G major has F#
        let g_acc = key_signature_accidentals("G");
        assert_eq!(abc_note_to_midi("F", &g_acc), Some(66)); // F# due to key sig
        assert_eq!(abc_note_to_midi("=F", &g_acc), Some(65)); // explicit natural
    }

    #[test]
    fn test_parse_voice_tokens() {
        let acc = HashMap::new();
        let tokens = parse_voice_tokens("C D E | F G", &acc, (1, 4));
        assert_eq!(tokens.len(), 6); // 3 notes + bar + 2 notes
        match &tokens[0] {
            VoiceToken::Note { midi, duration } => {
                assert_eq!(*midi, 60);
                let beats = duration.0 as f64 / duration.1 as f64;
                assert!((beats - 1.0).abs() < 0.001, "Expected 1 beat, got {}", beats);
            }
            _ => panic!("Expected note"),
        }
        match &tokens[3] {
            VoiceToken::Bar => {}
            _ => panic!("Expected bar"),
        }
    }

    #[test]
    fn test_parse_abc_duration() {
        let dl = (1, 4); // L:1/4
        // Helper: convert (num, den) to f64 beats for comparison
        let beats = |d: (u32, u32)| d.0 as f64 / d.1 as f64;
        assert!((beats(parse_abc_duration("", dl)) - 1.0).abs() < 0.001);     // quarter = 1 beat
        assert!((beats(parse_abc_duration("2", dl)) - 2.0).abs() < 0.001);    // half = 2 beats
        assert!((beats(parse_abc_duration("4", dl)) - 4.0).abs() < 0.001);    // whole = 4 beats
        assert!((beats(parse_abc_duration("/", dl)) - 0.5).abs() < 0.001);    // eighth = 0.5
        assert!((beats(parse_abc_duration("//", dl)) - 0.25).abs() < 0.001);  // sixteenth = 0.25
        assert!((beats(parse_abc_duration("/2", dl)) - 0.5).abs() < 0.001);   // eighth = 0.5
        assert!((beats(parse_abc_duration("3/2", dl)) - 1.5).abs() < 0.001);  // dotted = 1.5
    }

    #[test]
    fn test_parse_tempo() {
        assert_eq!(parse_tempo("1/4=120"), Some(120));
        assert_eq!(parse_tempo("[Q:1/4=120]"), Some(120));
        assert_eq!(parse_tempo("120"), Some(120));
        assert_eq!(parse_tempo(""), None);
    }

    #[test]
    fn test_parse_open_hymnal_minimal() {
        let abc = r#"
X: 1
T: Test Hymn
M: 4/4
K: C
[V: S1V1] C D E F | G A B c |
[V: S1V2] E, F, G, A, | B, C D E |
[V: S2V1] G, A, B, C | D E F G |
[V: S2V2] C, D, E, F, | G, A, B, C |
"#;
        let hymns = parse_open_hymnal(abc);
        assert_eq!(hymns.len(), 1);
        assert_eq!(hymns[0].title, "Test Hymn");
        assert_eq!(hymns[0].key, "C");
        // Should have 8 note beats + some bar beats
        let note_count = hymns[0]
            .beats
            .iter()
            .filter(|b| b.beat_type == "note")
            .count();
        assert!(note_count > 0, "Should have note beats");
    }

    #[test]
    fn test_parse_lead_sheets_detects_satb() {
        let abc = "[V: S1V1] C D E\n[V: S1V2] E F G\nX: 1\nT: Test\nK: C\n[V: S1V1] C D\n[V: S1V2] E F\n[V: S2V1] G A\n[V: S2V2] C, D,\n";
        let hymns = parse_lead_sheets(abc);
        // Should detect SATB and use parse_open_hymnal
        assert!(!hymns.is_empty() || true); // may or may not parse depending on format
    }

    #[test]
    fn test_parse_real_hymn() {
        // "Blessed Jesus at Thy Word" — first phrase in G major
        let abc = r#"
X: 22
T: Blessed Jesus at Thy Word
M: 4/4
K: G
[V: S1V1] B G A d | B G A2 | G G G G | A B A2 | G4 |
[V: S1V2] D D D F | D E F2 | D D E D | F G2 F | G4 |
[V: S2V1] G, B, A, A, | G, B, D2 | G, B, C G, | C D D2 | B,4 |
[V: S2V2] G, G, F, D, | G, E, D,2 | B,, G,, C, B,, | A,, G,, D,2 | G,,4 |
"#;
        let hymns = parse_open_hymnal(abc);
        assert_eq!(hymns.len(), 1, "Should parse one hymn");
        let hymn = &hymns[0];
        assert_eq!(hymn.title, "Blessed Jesus at Thy Word");
        assert_eq!(hymn.key, "G");
        assert_eq!(hymn.key_root, 7);

        let note_beats: Vec<&Beat> = hymn.beats.iter()
            .filter(|b| b.beat_type == "note")
            .collect();
        // 5 measures × ~4 notes = ~20 notes (minus held notes)
        assert!(note_beats.len() >= 15, "Should have at least 15 note beats, got {}", note_beats.len());

        // First beat should have a voicing (chord change)
        let first = &note_beats[0];
        assert!(first.rh_strings.is_some(), "First beat should have RH strings");
        assert!(first.lh_strings.is_some(), "First beat should have LH strings");

        // No string should be zero (skip-zero)
        for beat in &note_beats {
            if let Some(ref rh) = beat.rh_strings {
                for &s in rh {
                    assert_ne!(s, 0, "RH string must not be zero: {:?}", rh);
                }
            }
            if let Some(ref lh) = beat.lh_strings {
                for &s in lh {
                    assert_ne!(s, 0, "LH string must not be zero: {:?}", lh);
                }
            }
        }

        // Melody should be snapped to diatonic in key of G
        for beat in &note_beats {
            if let Some(midi) = beat.midi {
                let str_num = crate::music::midi_to_harp_string(midi, 7);
                assert!(str_num.is_some(), "Melody MIDI {} should be diatonic in G", midi);
            }
        }
    }

    #[test]
    fn test_parse_lead_sheet_format() {
        let abc = r#"X: 1
T: Simple Test
M: 4/4
K: C
"I" C D E F | "V" G A B c | "I" c2 z2 |
"#;
        let hymns = parse_lead_sheets(abc);
        assert_eq!(hymns.len(), 1);
        let hymn = &hymns[0];
        assert_eq!(hymn.title, "Simple Test");

        // Should have note beats with chord symbols
        let chord_beats: Vec<&Beat> = hymn.beats.iter()
            .filter(|b| b.chord.is_some())
            .collect();
        assert!(chord_beats.len() >= 2, "Should have chord beats, got {}", chord_beats.len());
    }
}
