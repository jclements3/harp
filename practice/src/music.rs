/// Note names in chromatic order
pub const NOTE_NAMES: [&str; 12] = [
    "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B",
];

/// A detected musical note
#[derive(Debug, Clone, PartialEq)]
pub struct Note {
    pub midi: u8,
    pub name: &'static str,
    pub octave: i8,
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

// ── Harp pedal system ──────────────────────────────────────────────

pub const PEDAL_LETTERS: [char; 7] = ['D', 'C', 'B', 'E', 'F', 'G', 'A'];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PedalPos {
    Flat,
    Natural,
    Sharp,
}

const BASE_SEMITONES: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];

pub const PEDAL_LETTER_INDEX: [usize; 7] = [1, 0, 6, 2, 3, 4, 5];

#[derive(Debug, Clone)]
pub struct PedalState {
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
    pub fn semitone_for_letter(&self, letter_idx: usize) -> u8 {
        let base = BASE_SEMITONES[letter_idx];
        let pedal_idx = PEDAL_LETTER_INDEX.iter().position(|&li| li == letter_idx).unwrap();
        match self.positions[pedal_idx] {
            PedalPos::Flat => (base + 11) % 12,
            PedalPos::Natural => base,
            PedalPos::Sharp => (base + 1) % 12,
        }
    }

    pub fn available_pitch_classes(&self) -> [u8; 7] {
        let mut pcs = [0u8; 7];
        for i in 0..7 {
            pcs[i] = self.semitone_for_letter(i);
        }
        pcs
    }

    pub fn cycle(&mut self, pedal_idx: usize) {
        self.positions[pedal_idx] = match self.positions[pedal_idx] {
            PedalPos::Flat => PedalPos::Natural,
            PedalPos::Natural => PedalPos::Sharp,
            PedalPos::Sharp => PedalPos::Flat,
        };
    }

    pub fn key_sig_info(&self) -> (u8, u8, bool) {
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6];
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];

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

    pub fn key_root(&self) -> u8 {
        let (sharps, flats, is_standard) = self.key_sig_info();
        if !is_standard {
            return 0;
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

    pub fn add_sharp(&mut self) {
        let (sharps, flats, _) = self.key_sig_info();
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6];
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];
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

    pub fn add_flat(&mut self) {
        let (sharps, flats, _) = self.key_sig_info();
        let sharp_order_letters: [usize; 7] = [3, 0, 4, 1, 5, 2, 6];
        let flat_order_letters: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];
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

pub const DIATONIC: [u8; 12] = [0, 0, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6];
pub const IS_SHARP: [bool; 12] = [false, true, false, true, false, false, true, false, true, false, true, false];
pub const MAJOR_SCALE: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];
pub const SOLFEGE: [&str; 7] = ["Do", "Re", "Mi", "Fa", "Sol", "La", "Ti"];

pub fn scale_degree(midi: u8, key_root: u8) -> u8 {
    let interval = ((midi % 12) as i8 - key_root as i8).rem_euclid(12) as u8;
    match MAJOR_SCALE.iter().position(|&s| s == interval) {
        Some(idx) => (idx as u8) + 1,
        None => 0,
    }
}

pub fn midi_to_staff_pos(midi: u8) -> i32 {
    let octave = (midi / 12) as i32 - 1;
    let pc = (midi % 12) as usize;
    octave * 7 + DIATONIC[pc] as i32
}

/// Build all MIDI notes available with current pedals in range, sorted ascending.
pub fn available_notes_in_range(pedals: &PedalState, low_midi: u8, high_midi: u8) -> Vec<u8> {
    let pcs = pedals.available_pitch_classes();
    (low_midi..=high_midi).filter(|m| pcs.contains(&(m % 12))).collect()
}

// ── Interval names ─────────────────────────────────────────────────

pub const INTERVAL_NAMES: [&str; 13] = [
    "Unison", "m2", "M2", "m3", "M3", "P4", "Tritone",
    "P5", "m6", "M6", "m7", "M7", "Octave",
];

pub fn interval_name(semitones: u8) -> &'static str {
    if (semitones as usize) < INTERVAL_NAMES.len() {
        INTERVAL_NAMES[semitones as usize]
    } else {
        "?"
    }
}
