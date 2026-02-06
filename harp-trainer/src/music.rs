/// Note names in chromatic order
pub const NOTE_NAMES: [&str; 12] = [
    "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B",
];

/// A detected musical note
#[derive(Debug, Clone, PartialEq)]
pub struct Note {
    /// MIDI note number (69 = A4)
    pub midi: u8,
    /// Note name (e.g. "C#")
    pub name: &'static str,
    /// Octave number
    pub octave: i8,
    /// How many cents off from the exact pitch (-50 to +50)
    pub cents_off: f32,
}

/// Convert a frequency in Hz to a Note, or None if out of range
pub fn freq_to_note(freq: f32) -> Option<Note> {
    if freq < 20.0 || freq > 5000.0 {
        return None;
    }
    let midi_float = 69.0 + 12.0 * (freq / 440.0).log2();
    let midi_rounded = midi_float.round();
    let midi = midi_rounded as u8;
    let cents_off = (midi_float - midi_rounded) * 100.0;
    let note_idx = (midi % 12) as usize;
    let octave = (midi as i8 / 12) - 1;

    Some(Note {
        midi,
        name: NOTE_NAMES[note_idx],
        octave,
        cents_off,
    })
}

/// Semitone distance between two MIDI notes (0-11)
pub fn interval(a: u8, b: u8) -> u8 {
    ((b as i16 - a as i16).rem_euclid(12)) as u8
}

#[derive(Debug, Clone)]
pub struct ChordType {
    pub name: &'static str,
    pub intervals: &'static [u8],
}

pub const CHORD_TYPES: &[ChordType] = &[
    ChordType { name: "Major",    intervals: &[0, 4, 7] },
    ChordType { name: "Minor",    intervals: &[0, 3, 7] },
    ChordType { name: "Dim",      intervals: &[0, 3, 6] },
    ChordType { name: "Maj7",     intervals: &[0, 4, 7, 11] },
    ChordType { name: "Min7",     intervals: &[0, 3, 7, 10] },
    ChordType { name: "Dom7",     intervals: &[0, 4, 7, 10] },
    ChordType { name: "Half-dim", intervals: &[0, 3, 6, 10] },
];

/// Given a set of detected MIDI pitch classes (0-11), find all matching chord types
pub fn match_chords(pitch_classes: &[u8]) -> Vec<(&'static str, &'static str)> {
    let mut results = Vec::new();
    for root_pc in pitch_classes {
        for chord in CHORD_TYPES {
            let all_match = pitch_classes.iter().all(|&pc| {
                let iv = interval(*root_pc, pc);
                chord.intervals.contains(&iv)
            });
            if all_match {
                let root_name = NOTE_NAMES[*root_pc as usize];
                results.push((root_name, chord.name));
            }
        }
    }
    results
}

// ── Harp pedal system ──────────────────────────────────────────────

/// The 7 harp pedal letters in standard order (left foot: D,C,B | right foot: E,F,G,A)
pub const PEDAL_LETTERS: [char; 7] = ['D', 'C', 'B', 'E', 'F', 'G', 'A'];

/// Order in which sharps are added to key signatures
pub const SHARP_ORDER: [usize; 7] = [5, 0, 4, 1, 6, 3, 2]; // F,C,G,D,A,E,B indices into PEDAL_LETTERS→ mapped to letter indices
// Actually let's use letter indices into the 7 diatonic letters C=0,D=1,E=2,F=3,G=4,A=5,B=6

/// Pedal position
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PedalPos {
    Flat,
    Natural,
    Sharp,
}

/// Base semitones for each diatonic letter: C D E F G A B
const BASE_SEMITONES: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];

/// Letter index for each pedal letter (D=1, C=0, B=6, E=2, F=3, G=4, A=5)
pub const PEDAL_LETTER_INDEX: [usize; 7] = [1, 0, 6, 2, 3, 4, 5];

/// Harp pedal state — 7 pedals, each flat/natural/sharp
#[derive(Debug, Clone)]
pub struct PedalState {
    /// Positions indexed by PEDAL_LETTERS order: D, C, B, E, F, G, A
    pub positions: [PedalPos; 7],
}

impl Default for PedalState {
    fn default() -> Self {
        Self {
            positions: [PedalPos::Natural; 7],
        }
    }
}

impl PedalState {
    /// Get the semitone (0-11) for a given letter index (C=0..B=6)
    pub fn semitone_for_letter(&self, letter_idx: usize) -> u8 {
        let base = BASE_SEMITONES[letter_idx];
        // Find which pedal controls this letter
        let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
        match self.positions[pedal_idx] {
            PedalPos::Flat => (base + 11) % 12,
            PedalPos::Natural => base,
            PedalPos::Sharp => (base + 1) % 12,
        }
    }

    /// Get all 7 available pitch classes from current pedal configuration
    pub fn available_pitch_classes(&self) -> [u8; 7] {
        let mut pcs = [0u8; 7];
        for i in 0..7 {
            pcs[i] = self.semitone_for_letter(i);
        }
        pcs
    }

    /// Cycle a pedal through flat → natural → sharp → flat
    pub fn cycle(&mut self, pedal_idx: usize) {
        self.positions[pedal_idx] = match self.positions[pedal_idx] {
            PedalPos::Flat => PedalPos::Natural,
            PedalPos::Natural => PedalPos::Sharp,
            PedalPos::Sharp => PedalPos::Flat,
        };
    }

    /// Count sharps and flats; return (sharps, flats, is_standard_key)
    pub fn key_sig_info(&self) -> (u8, u8, bool) {
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6]; // F,C,G,D,A,E,B
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3]; // B,E,A,D,G,C,F

        let mut sharps = 0u8;
        let mut flats = 0u8;
        for i in 0..7 {
            let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == i).unwrap();
            match self.positions[pedal_idx] {
                PedalPos::Sharp => sharps += 1,
                PedalPos::Flat => flats += 1,
                PedalPos::Natural => {}
            }
        }

        if sharps > 0 && flats > 0 {
            return (sharps, flats, false);
        }

        let mut is_standard = true;
        if sharps > 0 {
            for i in 0..7 {
                let letter_idx = sharp_order_letters[i];
                let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
                let expected = if (i as u8) < sharps { PedalPos::Sharp } else { PedalPos::Natural };
                if self.positions[pedal_idx] != expected {
                    is_standard = false;
                    break;
                }
            }
        }
        if flats > 0 {
            for i in 0..7 {
                let letter_idx = flat_order_letters[i];
                let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
                let expected = if (i as u8) < flats { PedalPos::Flat } else { PedalPos::Natural };
                if self.positions[pedal_idx] != expected {
                    is_standard = false;
                    break;
                }
            }
        }

        (sharps, flats, is_standard)
    }

    /// Derive the key root (pitch class 0-11) from a standard pedal configuration
    pub fn key_root(&self) -> u8 {
        let (sharps, flats, is_standard) = self.key_sig_info();
        if !is_standard {
            return 0; // default to C
        }
        let sharp_roots: [u8; 8] = [0, 7, 2, 9, 4, 11, 6, 1];
        let flat_roots: [u8; 8] = [0, 5, 10, 3, 8, 1, 6, 11];
        if sharps > 0 {
            sharp_roots[sharps as usize]
        } else if flats > 0 {
            flat_roots[flats as usize]
        } else {
            0
        }
    }

    /// Add a sharp (or remove a flat) to move around the circle of fifths
    pub fn add_sharp(&mut self) {
        let (sharps, flats, _) = self.key_sig_info();
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6]; // F,C,G,D,A,E,B
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3]; // B,E,A,D,G,C,F
        if flats > 0 {
            let letter_idx = flat_order_letters[(flats - 1) as usize];
            let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
            self.positions[pedal_idx] = PedalPos::Natural;
        } else if sharps < 7 {
            let letter_idx = sharp_order_letters[sharps as usize];
            let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
            self.positions[pedal_idx] = PedalPos::Sharp;
        }
    }

    /// Add a flat (or remove a sharp) to move around the circle of fifths
    pub fn add_flat(&mut self) {
        let (sharps, flats, _) = self.key_sig_info();
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6]; // F,C,G,D,A,E,B
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3]; // B,E,A,D,G,C,F
        if sharps > 0 {
            let letter_idx = sharp_order_letters[(sharps - 1) as usize];
            let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
            self.positions[pedal_idx] = PedalPos::Natural;
        } else if flats < 7 {
            let letter_idx = flat_order_letters[flats as usize];
            let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
            self.positions[pedal_idx] = PedalPos::Flat;
        }
    }

    /// Key name string (e.g. "C major", "G major", "Custom")
    pub fn key_name(&self) -> String {
        let (sharps, flats, is_standard) = self.key_sig_info();
        if !is_standard {
            return "Custom".to_string();
        }
        let sharp_keys = ["C major", "G major", "D major", "A major", "E major", "B major", "F\u{266f} major", "C\u{266f} major"];
        let flat_keys = ["C major", "F major", "B\u{266d} major", "E\u{266d} major", "A\u{266d} major", "D\u{266d} major", "G\u{266d} major", "C\u{266d} major"];
        if sharps > 0 {
            sharp_keys[sharps as usize].to_string()
        } else if flats > 0 {
            flat_keys[flats as usize].to_string()
        } else {
            "C major".to_string()
        }
    }
}

// ── Staff position & scale degree ──────────────────────────────────

/// Maps pitch class (0-11) to diatonic position within an octave (C=0,D=1,E=2,F=3,G=4,A=5,B=6)
pub const DIATONIC: [u8; 12] = [0, 0, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6];

/// Which pitch classes are sharped (not on a white key)
pub const IS_SHARP: [bool; 12] = [false, true, false, true, false, false, true, false, true, false, true, false];

/// Major scale intervals
pub const MAJOR_SCALE: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];

pub const SOLFEGE: [&str; 7] = ["Do", "Re", "Mi", "Fa", "Sol", "La", "Ti"];

/// Scale degree 1-7 for diatonic notes, 0 for chromatic
pub fn scale_degree(midi: u8, key_root: u8) -> u8 {
    let interval = ((midi % 12) as i8 - key_root as i8).rem_euclid(12) as u8;
    match MAJOR_SCALE.iter().position(|&s| s == interval) {
        Some(idx) => (idx as u8) + 1,
        None => 0,
    }
}

/// Diatonic staff position. Middle C (MIDI 60) = 28, matching the HTML version.
pub fn midi_to_staff_pos(midi: u8) -> i32 {
    let octave = (midi / 12) as i32 - 1;
    let pc = (midi % 12) as usize;
    octave * 7 + DIATONIC[pc] as i32
}

/// Generate drill notes from the harp's currently available scale (pedal-derived).
/// Returns `count` random MIDI notes within [low_midi, high_midi] using only
/// the 7 pitch classes available from the pedal configuration.
pub fn generate_drill_notes(pedals: &PedalState, low_midi: u8, high_midi: u8, count: usize) -> Vec<u8> {
    use rand::seq::SliceRandom;

    // Build list of all playable MIDI notes in range
    let pcs = pedals.available_pitch_classes();
    let mut candidates: Vec<u8> = Vec::new();
    for midi in low_midi..=high_midi {
        let pc = midi % 12;
        if pcs.contains(&pc) {
            candidates.push(midi);
        }
    }

    if candidates.is_empty() {
        return vec![60; count]; // fallback
    }

    let mut rng = rand::thread_rng();
    let mut notes = Vec::with_capacity(count);
    // Avoid consecutive duplicates
    for _ in 0..count {
        loop {
            let &n = candidates.choose(&mut rng).unwrap();
            if notes.last() != Some(&n) {
                notes.push(n);
                break;
            }
        }
    }
    notes
}

// ── Hand position system ────────────────────────────────────────────

/// A hand position: consecutive harp strings for finger placement
#[derive(Debug, Clone)]
pub struct HandPosition {
    pub notes: Vec<u8>,
}

/// Build all MIDI notes available with current pedals in range, sorted ascending.
pub fn available_notes_in_range(pedals: &PedalState, low_midi: u8, high_midi: u8) -> Vec<u8> {
    let pcs = pedals.available_pitch_classes();
    (low_midi..=high_midi).filter(|m| pcs.contains(&(m % 12))).collect()
}

/// Generate a hand position of `fingers` consecutive available harp strings.
/// Biases toward positions containing `weak_notes` (slow response notes).
/// Avoids repeating `last_position`.
pub fn generate_hand_position(
    available: &[u8],
    fingers: usize,
    weak_notes: &[u8],
    last_position: Option<&HandPosition>,
) -> Option<HandPosition> {
    use rand::Rng;

    if available.len() < fingers {
        return None;
    }

    let max_start = available.len() - fingers;
    let mut rng = rand::thread_rng();

    // Build weights for each possible starting position
    let mut weights: Vec<f32> = vec![1.0; max_start + 1];
    for (i, w) in weights.iter_mut().enumerate() {
        let pos_notes = &available[i..i + fingers];
        // Boost positions containing weak notes
        for &weak in weak_notes {
            if pos_notes.contains(&weak) {
                *w += 3.0;
            }
        }
        // Penalize repeating the same position
        if let Some(last) = last_position {
            if pos_notes == last.notes.as_slice() {
                *w = 0.0;
            }
        }
    }

    let total: f32 = weights.iter().sum();
    if total <= 0.0 {
        let start = rng.gen_range(0..=max_start);
        return Some(HandPosition { notes: available[start..start + fingers].to_vec() });
    }

    let mut r = rng.gen_range(0.0..total);
    for (i, &w) in weights.iter().enumerate() {
        r -= w;
        if r <= 0.0 {
            return Some(HandPosition { notes: available[i..i + fingers].to_vec() });
        }
    }

    Some(HandPosition { notes: available[max_start..max_start + fingers].to_vec() })
}
