use serde::Serialize;

/// Maximum diatonic span per hand (10 strings).
pub const MAX_HAND_SPAN: i16 = 10;

/// Maximum fingers per hand.
pub const FINGERS_PER_HAND: usize = 4;

/// Harp range: C1 (MIDI 24) to G7 (MIDI 103).
pub const HARP_LOW_MIDI: u8 = 24;
pub const HARP_HIGH_MIDI: u8 = 103;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Hand {
    Right,
    Left,
}

/// A single note assigned to a finger.
#[derive(Debug, Clone, Serialize)]
pub struct FingerNote {
    pub midi: u8,
    pub string_num: i16,
    pub finger: u8, // 1=thumb, 2=index, 3=middle, 4=ring
    pub hand: Hand,
    pub note_name: String,
    pub is_root: bool,
}

/// A complete voicing across both hands.
#[derive(Debug, Clone, Serialize)]
pub struct Voicing {
    pub rh: Vec<FingerNote>,  // sorted high to low (thumb first)
    pub lh: Vec<FingerNote>,  // sorted high to low (thumb first)
    pub score: f32,
}

impl Voicing {
    pub fn total_notes(&self) -> usize {
        self.rh.len() + self.lh.len()
    }

    pub fn all_midis(&self) -> Vec<u8> {
        let mut midis: Vec<u8> = self.rh.iter().chain(self.lh.iter())
            .map(|n| n.midi)
            .collect();
        midis.sort();
        midis
    }

    pub fn lowest_midi(&self) -> Option<u8> {
        self.all_midis().first().copied()
    }

    pub fn highest_midi(&self) -> Option<u8> {
        self.all_midis().last().copied()
    }

    pub fn rh_lowest(&self) -> Option<u8> {
        self.rh.last().map(|n| n.midi)
    }

    pub fn lh_highest(&self) -> Option<u8> {
        self.lh.first().map(|n| n.midi)
    }
}
