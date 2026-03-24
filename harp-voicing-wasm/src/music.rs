/// Diatonic major scale as semitone offsets from root.
pub const MAJOR_SCALE: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];

/// Note names for each pitch class (0-11).
pub const PC_NAMES: [&str; 12] = [
    "C", "Db", "D", "Eb", "E", "F", "F#", "G", "Ab", "A", "Bb", "B",
];

/// Convert scale degree (0-6) + key root pitch class (0-11) to pitch class.
pub fn degree_to_pc(degree: u8, key_root: u8) -> u8 {
    (key_root + MAJOR_SCALE[degree as usize % 7]) % 12
}

/// Pitch class of a MIDI note.
pub fn pitch_class(midi: u8) -> u8 {
    midi % 12
}

/// Convert a MIDI note to a diatonic harp string number (0-indexed from lowest).
/// Returns None if the note isn't on a diatonic string in the given key.
/// Uses key-relative octave calculation so strings are correctly ordered in all keys.
pub fn midi_to_harp_string(midi: u8, key_root: u8) -> Option<i16> {
    let pc = pitch_class(midi);
    let degree = MAJOR_SCALE.iter().position(|&s| (key_root + s) % 12 == pc)?;
    // Calculate octave relative to the key root, not absolute MIDI octave.
    // This ensures notes within the same scale octave get consecutive string numbers
    // even when the scale wraps around a MIDI octave boundary (e.g., key of G: G..F#).
    let semi_above_root = (pc as i16 - key_root as i16 + 12) % 12;
    let midi_of_root_in_same_octave_group = midi as i16 - semi_above_root;
    let root_octave = midi_of_root_in_same_octave_group / 12 - 1;
    Some(root_octave * 7 + degree as i16)
}

/// Convert a harp string number back to MIDI note.
pub fn harp_string_to_midi(string_num: i16, key_root: u8) -> u8 {
    let octave = string_num.div_euclid(7);
    let degree = string_num.rem_euclid(7) as usize;
    let semis_above_root = MAJOR_SCALE[degree];
    ((octave + 1) * 12 + key_root as i16 + semis_above_root as i16) as u8
}

/// Format a MIDI note as a name like "C4", "G#3".
pub fn midi_to_name(midi: u8) -> String {
    let pc = pitch_class(midi);
    let octave = (midi / 12) as i8 - 1;
    format!("{}{}", PC_NAMES[pc as usize], octave)
}

/// Find all MIDI notes with a given pitch class within a range.
pub fn all_octaves(pc: u8, low: u8, high: u8) -> Vec<u8> {
    let mut notes = Vec::new();
    // Start from lowest octave
    let mut midi = pc;
    if midi < low {
        let diff = low - midi;
        midi += ((diff + 11) / 12) * 12;
    }
    while midi <= high {
        notes.push(midi);
        midi += 12;
    }
    notes
}

/// Snap a chromatic MIDI note to the nearest diatonic string in the given key.
/// If the note is already diatonic, returns it unchanged.
/// Otherwise tries 1-2 semitones up/down, preferring upward resolution.
pub fn snap_to_diatonic(midi: u8, key_root: u8) -> u8 {
    if midi_to_harp_string(midi, key_root).is_some() {
        return midi;
    }
    for &offset in &[1i8, -1, 2, -2] {
        let candidate = (midi as i8 + offset) as u8;
        if midi_to_harp_string(candidate, key_root).is_some() {
            return candidate;
        }
    }
    midi // give up
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_degree_to_pc() {
        // Key of C (root=0): I=C(0), ii=D(2), iii=E(4), IV=F(5), V=G(7)
        assert_eq!(degree_to_pc(0, 0), 0); // C
        assert_eq!(degree_to_pc(1, 0), 2); // D
        assert_eq!(degree_to_pc(4, 0), 7); // G
        // Key of G (root=7): I=G(7), ii=A(9), V=D(2)
        assert_eq!(degree_to_pc(0, 7), 7); // G
        assert_eq!(degree_to_pc(4, 7), 2); // D
    }

    #[test]
    fn test_harp_string_roundtrip() {
        // Key of C: MIDI 60 = C4
        let s = midi_to_harp_string(60, 0).unwrap();
        assert_eq!(harp_string_to_midi(s, 0), 60);
        // MIDI 64 = E4
        let s = midi_to_harp_string(64, 0).unwrap();
        assert_eq!(harp_string_to_midi(s, 0), 64);
    }

    #[test]
    fn test_harp_string_key_of_g() {
        // In key of G (root=7): F#5 (78) must be below G5 (79)
        let fsharp5 = midi_to_harp_string(78, 7).unwrap(); // F#5 = degree 6
        let g5 = midi_to_harp_string(79, 7).unwrap();      // G5 = degree 0 of next octave
        assert!(fsharp5 < g5, "F#5 ({}) should be below G5 ({})", fsharp5, g5);
        // Roundtrip
        assert_eq!(harp_string_to_midi(fsharp5, 7), 78);
        assert_eq!(harp_string_to_midi(g5, 7), 79);
    }

    #[test]
    fn test_harp_string_roundtrip_all_keys() {
        // Test roundtrip for all standard keys across a wide range
        for key_root in [0, 2, 4, 5, 7, 9, 11] {
            for midi in 24..=103u8 {
                if let Some(s) = midi_to_harp_string(midi, key_root) {
                    let back = harp_string_to_midi(s, key_root);
                    assert_eq!(back, midi,
                        "Roundtrip failed: key={}, midi={}, string={}, back={}",
                        key_root, midi, s, back);
                }
            }
        }
    }

    #[test]
    fn test_all_octaves() {
        let notes = all_octaves(0, 48, 84); // all C's from C3 to C6
        assert_eq!(notes, vec![48, 60, 72, 84]);
    }

    #[test]
    fn test_midi_to_name() {
        assert_eq!(midi_to_name(60), "C4");
        assert_eq!(midi_to_name(71), "B4");
        assert_eq!(midi_to_name(48), "C3");
    }
}
