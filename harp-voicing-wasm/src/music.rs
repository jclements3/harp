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
pub fn midi_to_harp_string(midi: u8, key_root: u8) -> Option<u8> {
    let pc = pitch_class(midi);
    let octave = (midi / 12) as i16 - 1; // MIDI octave numbering
    let degree = MAJOR_SCALE.iter().position(|&s| (key_root + s) % 12 == pc)?;
    Some((octave as u8) * 7 + degree as u8)
}

/// Convert a harp string number back to MIDI note.
pub fn harp_string_to_midi(string_num: u8, key_root: u8) -> u8 {
    let octave = string_num / 7;
    let degree = string_num % 7;
    let pc = (key_root + MAJOR_SCALE[degree as usize]) % 12;
    (octave + 1) * 12 + pc
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
