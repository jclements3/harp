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

/// Convert a relative string number (skips zero: ...-2,-1,1,2,3...) to staff position.
/// Staff position 0 = middle C.
pub fn relative_string_to_staff_pos(rel_str: i32, key_root: i32) -> i32 {
    // Undo zero-skip to get linear index (0-based from lowest root)
    let linear = if rel_str > 0 { rel_str - 1 } else { rel_str };

    // ref_str = absolute string number of the lowest root on the harp
    let lowest_root_midi = all_octaves(key_root, HARP_LOW_MIDI, HARP_HIGH_MIDI)
        .into_iter().next().unwrap_or(HARP_LOW_MIDI);
    let ref_str = midi_to_harp_string(lowest_root_midi, key_root).unwrap_or(0);

    // Absolute string = ref_str + linear
    let abs_str = ref_str + linear;

    // Middle C absolute string
    let mc_str = midi_to_harp_string(60, key_root).unwrap_or(21);

    // +7 (one octave up) to match standard harp display register
    abs_str - mc_str + 7
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

// ── Pedal state (for drill challenge generation) ──────────────────

pub const NOTE_NAMES: [&str; 12] = ["C","C#","D","D#","E","F","F#","G","G#","A","A#","B"];

const BASE_SEMI: [u8; 7] = [0, 2, 4, 5, 7, 9, 11];
const PEDAL_LETTER_IDX: [usize; 7] = [1, 0, 6, 2, 3, 4, 5]; // D,C,B,E,F,G,A -> letter index

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PedalPos { Flat, Natural, Sharp }

#[derive(Debug, Clone)]
pub struct PedalState {
    pub positions: [PedalPos; 7],
}

impl Default for PedalState {
    fn default() -> Self { Self { positions: [PedalPos::Natural; 7] } }
}

impl PedalState {
    pub fn available_pitch_classes(&self) -> [u8; 7] {
        let mut pcs = [0u8; 7];
        for i in 0..7 {
            let base = BASE_SEMI[i];
            let pidx = PEDAL_LETTER_IDX.iter().position(|&li| li == i).unwrap();
            pcs[i] = match self.positions[pidx] {
                PedalPos::Flat => (base + 11) % 12,
                PedalPos::Natural => base,
                PedalPos::Sharp => (base + 1) % 12,
            };
        }
        pcs
    }

    pub fn add_sharp(&mut self) {
        let sharp_order: [usize; 7] = [3, 0, 4, 1, 5, 2, 6];
        let flat_order: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];
        let (sharps, flats, _) = self.key_sig_info();
        if flats > 0 {
            let li = flat_order[(flats-1) as usize];
            let pi = PEDAL_LETTER_IDX.iter().position(|&l| l == li).unwrap();
            self.positions[pi] = PedalPos::Natural;
        } else if sharps < 7 {
            let li = sharp_order[sharps as usize];
            let pi = PEDAL_LETTER_IDX.iter().position(|&l| l == li).unwrap();
            self.positions[pi] = PedalPos::Sharp;
        }
    }

    pub fn add_flat(&mut self) {
        let sharp_order: [usize; 7] = [3, 0, 4, 1, 5, 2, 6];
        let flat_order: [usize; 7] = [6, 2, 5, 1, 4, 0, 3];
        let (sharps, flats, _) = self.key_sig_info();
        if sharps > 0 {
            let li = sharp_order[(sharps-1) as usize];
            let pi = PEDAL_LETTER_IDX.iter().position(|&l| l == li).unwrap();
            self.positions[pi] = PedalPos::Natural;
        } else if flats < 7 {
            let li = flat_order[flats as usize];
            let pi = PEDAL_LETTER_IDX.iter().position(|&l| l == li).unwrap();
            self.positions[pi] = PedalPos::Flat;
        }
    }

    fn key_sig_info(&self) -> (u8, u8, bool) {
        let mut sharps = 0u8;
        let mut flats = 0u8;
        for i in 0..7 {
            let pi = PEDAL_LETTER_IDX.iter().position(|&l| l == i).unwrap();
            match self.positions[pi] {
                PedalPos::Sharp => sharps += 1,
                PedalPos::Flat => flats += 1,
                PedalPos::Natural => {}
            }
        }
        (sharps, flats, !(sharps > 0 && flats > 0))
    }
}

pub fn available_notes_in_range(pedals: &PedalState, low: u8, high: u8) -> Vec<u8> {
    let pcs = pedals.available_pitch_classes();
    (low..=high).filter(|m| pcs.contains(&(m % 12))).collect()
}

pub const INTERVAL_NAMES: [&str; 13] = [
    "Unison","m2","M2","m3","M3","P4","Tritone","P5","m6","M6","m7","M7","Octave",
];

pub fn interval_name(semi: u8) -> &'static str {
    if (semi as usize) < INTERVAL_NAMES.len() { INTERVAL_NAMES[semi as usize] } else { "?" }
}

pub fn pedals_from_key(key: &str) -> PedalState {
    let mut p = PedalState::default();
    for &(k,n) in &[("G",1),("D",2),("A",3),("E",4),("B",5),("F#",6)] {
        if k == key { for _ in 0..n { p.add_sharp(); } return p; }
    }
    for &(k,n) in &[("F",1),("Bb",2),("Eb",3),("Ab",4),("Db",5)] {
        if k == key { for _ in 0..n { p.add_flat(); } return p; }
    }
    p
}

/// Relative string number: string 1 = key root, negative below, skip zero
pub fn to_relative_string(string_num: i32, key_root: i32) -> i32 {
    let lowest_root_midi = all_octaves(key_root, HARP_LOW_MIDI, HARP_HIGH_MIDI)
        .into_iter().next().unwrap_or(HARP_LOW_MIDI);
    let ref_str = midi_to_harp_string(lowest_root_midi, key_root).unwrap_or(0);
    let rel = string_num - ref_str + 1;
    if rel <= 0 { rel - 1 } else { rel }
}
