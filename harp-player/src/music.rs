pub const MAJOR_SCALE: [i32; 7] = [0, 2, 4, 5, 7, 9, 11];
pub const PC_NAMES: [&str; 12] = ["C","Db","D","Eb","E","F","F#","G","Ab","A","Bb","B"];
pub const HARP_LOW_MIDI: i32 = 24;
pub const HARP_HIGH_MIDI: i32 = 103;
pub const MAX_HAND_SPAN: i32 = 10;

pub fn pitch_class(midi: i32) -> i32 {
    ((midi % 12) + 12) % 12
}

pub fn midi_to_harp_string(midi: i32, key_root: i32) -> Option<i32> {
    let pc = pitch_class(midi);
    let deg = MAJOR_SCALE.iter().position(|&s| (key_root + s) % 12 == pc)?;
    let semi_above_root = ((pc - key_root) % 12 + 12) % 12;
    let midi_of_root = midi - semi_above_root;
    let root_octave = (midi_of_root as f64 / 12.0).floor() as i32 - 1;
    Some(root_octave * 7 + deg as i32)
}

pub fn harp_string_to_midi(string_num: i32, key_root: i32) -> i32 {
    let octave = string_num.div_euclid(7);
    let degree = string_num.rem_euclid(7) as usize;
    let semis = MAJOR_SCALE[degree];
    (octave + 1) * 12 + key_root + semis
}

pub fn midi_to_name(midi: i32) -> String {
    let pc = pitch_class(midi) as usize;
    let octave = midi / 12 - 1;
    format!("{}{}", PC_NAMES[pc], octave)
}

pub fn all_octaves(pc: i32, low: i32, high: i32) -> Vec<i32> {
    let mut notes = Vec::new();
    let mut midi = pc;
    if midi < low {
        midi += ((low - midi + 11) / 12) * 12;
    }
    while midi <= high {
        notes.push(midi);
        midi += 12;
    }
    notes
}

pub fn snap_to_diatonic(midi: i32, key_root: i32) -> i32 {
    if midi_to_harp_string(midi, key_root).is_some() { return midi; }
    for &delta in &[1, -1, 2, -2] {
        if midi_to_harp_string(midi + delta, key_root).is_some() {
            return midi + delta;
        }
    }
    midi
}

/// Absolute harp string number for middle C (MIDI 60) in a given key.
/// Used as the reference point: staff_pos = abs_string - middle_c_string.
pub fn middle_c_string(key_root: i32) -> i32 {
    midi_to_harp_string(60, key_root).unwrap_or(21)
}

pub fn key_to_pc(key: &str) -> i32 {
    match key {
        "C" => 0, "Db" | "C#" => 1, "D" => 2, "Eb" => 3, "E" => 4,
        "F" => 5, "F#" | "Gb" => 6, "G" => 7, "Ab" => 8, "A" => 9,
        "Bb" => 10, "B" | "Cb" => 11, _ => 0,
    }
}

pub fn key_sig_accidentals(key: &str) -> [i32; 7] {
    // Returns accidental for each letter: C=0,D=1,E=2,F=3,G=4,A=5,B=6
    let sharp_order: [usize; 7] = [3, 0, 4, 1, 5, 2, 6]; // F,C,G,D,A,E,B
    let flat_order: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];  // B,E,A,D,G,C,F
    let mut acc = [0i32; 7];

    let sharp_keys: &[(&str, usize)] = &[
        ("G",1),("D",2),("A",3),("E",4),("B",5),("F#",6),("C#",7),
    ];
    let flat_keys: &[(&str, usize)] = &[
        ("F",1),("Bb",2),("Eb",3),("Ab",4),("Db",5),("Gb",6),("Cb",7),
    ];

    for &(k, n) in sharp_keys {
        if k == key {
            for i in 0..n { acc[sharp_order[i]] = 1; }
            return acc;
        }
    }
    for &(k, n) in flat_keys {
        if k == key {
            for i in 0..n { acc[flat_order[i]] = -1; }
            return acc;
        }
    }
    acc
}

/// Relative string number: string 1 = key root, negative below, skip zero
pub fn to_relative_string(string_num: i32, key_root: i32) -> i32 {
    let lowest_root_midi = all_octaves(key_root, HARP_LOW_MIDI, HARP_HIGH_MIDI)
        .into_iter().next().unwrap_or(HARP_LOW_MIDI);
    let ref_str = midi_to_harp_string(lowest_root_midi, key_root).unwrap_or(0);
    let rel = string_num - ref_str + 1;
    if rel <= 0 { rel - 1 } else { rel }
}
