/// Note names in chromatic order
const NOTE_NAMES: [&str; 12] = [
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
