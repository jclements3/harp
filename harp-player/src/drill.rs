//! Drill mode: generates random "hymns" and feeds them through the same
//! voicing + rendering pipeline as real hymns.

use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use crate::abc::{Score, ScoreEvent, compute_display_strings};
use crate::chord::identify_chord_names;
use crate::music;
use crate::voicing::{voice_from_satb_fingers, Voicing};

// ── Drill config ──────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DrillMode { Intervals, Chords }

#[derive(Debug, Clone)]
pub struct DrillConfig {
    pub mode: DrillMode,
    pub chord_size: usize,
}

impl Default for DrillConfig {
    fn default() -> Self { Self { mode: DrillMode::Intervals, chord_size: 3 } }
}

// ── Persistent progress ───────────────────────────────────────────

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TypeStats { pub attempts: u32, pub total_time: f32 }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Progress {
    pub level: u32,
    pub target_secs: f32,
    pub best_npm: f32,
    pub total_notes: u64,
    pub total_sessions: u32,
    pub best_streak: u32,
    pub total_time_secs: f32,
    #[serde(default)]
    pub history: Vec<SessionRecord>,
    #[serde(default)]
    pub type_stats: HashMap<String, TypeStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionRecord { pub npm: f32, pub level: u32, pub pace: f32, pub challenges: u32, pub notes: usize }

impl Default for Progress {
    fn default() -> Self {
        Self { level: 0, target_secs: 5.0, best_npm: 0.0, total_notes: 0, total_sessions: 0,
               best_streak: 0, total_time_secs: 0.0, history: Vec::new(), type_stats: HashMap::new() }
    }
}

fn progress_path() -> std::path::PathBuf {
    if cfg!(target_os = "android") {
        std::path::PathBuf::from("/data/data/com.harp.player/practice_progress.json")
    } else {
        let home = std::env::var("HOME").unwrap_or_else(|_| ".".into());
        std::path::PathBuf::from(home).join(".harp_player_drill_progress.json")
    }
}

pub fn load_progress() -> Progress {
    std::fs::read_to_string(progress_path()).ok()
        .and_then(|s| serde_json::from_str(&s).ok()).unwrap_or_default()
}

pub fn save_progress(p: &Progress) {
    if let Ok(j) = serde_json::to_string_pretty(p) { let _ = std::fs::write(progress_path(), j); }
}

// ── Generate a random Score (same format as a hymn) ───────────────

/// Generate a random drill Score — random SATB notes voiced through the
/// same engine as hymns. The result is a Score with ScoreEvents that
/// render_score() can display identically to a real hymn.
pub fn generate_random_score(
    key: &str,
    mode_offset: i32,
    num_measures: usize,
    beats_per_measure: u8,
    rh_fingers: usize,
    lh_fingers: usize,
) -> Score {
    use rand::Rng;
    let mut rng = rand::thread_rng();

    let key_root = {
        let pc = music::key_to_pc(key);
        ((pc - music::MAJOR_SCALE[mode_offset as usize]) % 12 + 12) % 12
    };

    // Build diatonic scale for this key across harp range
    let scale: Vec<i32> = (music::HARP_LOW_MIDI..=music::HARP_HIGH_MIDI)
        .filter(|&m| music::midi_to_harp_string(m, key_root).is_some())
        .collect();

    // Split into soprano range (C4-C6) and bass range (C2-B3)
    let soprano_range: Vec<i32> = scale.iter().copied().filter(|&m| m >= 60 && m <= 84).collect();
    let bass_range: Vec<i32> = scale.iter().copied().filter(|&m| m >= 36 && m <= 59).collect();

    let mut events = Vec::new();
    let mut last_voicing: Option<Voicing> = None;

    for measure in 0..num_measures {
        for beat in 0..beats_per_measure {
            // Pick random soprano note
            let s_midi = soprano_range[rng.gen_range(0..soprano_range.len())];

            // Pick random alto (near soprano, slightly lower)
            let alto_candidates: Vec<i32> = soprano_range.iter().copied()
                .filter(|&m| m <= s_midi && m >= s_midi - 12).collect();
            let a_midi = if !alto_candidates.is_empty() {
                Some(alto_candidates[rng.gen_range(0..alto_candidates.len())])
            } else { None };

            // Pick random tenor/bass from bass range
            let t_midi = if !bass_range.is_empty() {
                Some(bass_range[rng.gen_range(0..bass_range.len())])
            } else { None };
            let b_midi = if !bass_range.is_empty() {
                Some(bass_range[rng.gen_range(0..bass_range.len())])
            } else { None };

            let beats = 1.0;

            let voicing = voice_from_satb_fingers(
                s_midi, a_midi, t_midi, b_midi, key_root,
                last_voicing.as_ref(), rh_fingers, lh_fingers,
            );

            let snapped = music::snap_to_diatonic(s_midi, key_root);

            if let Some(ref v) = voicing {
                let voicing_key = format!("{:?}|{:?}", v.rh_strings, v.lh_strings);
                let last_key = last_voicing.as_ref()
                    .map(|lv| format!("{:?}|{:?}", lv.rh_strings, lv.lh_strings))
                    .unwrap_or_default();
                let is_chord_change = voicing_key != last_key;

                let rh_strs: Vec<i32> = v.rh_strings.iter()
                    .map(|&s| music::to_relative_string(s, key_root)).collect();
                let lh_strs: Vec<i32> = v.lh_strings.iter()
                    .map(|&s| music::to_relative_string(s, key_root)).collect();
                let melody_str = *rh_strs.first().unwrap_or(&1);

                events.push(ScoreEvent::Note {
                    melody_midi: snapped,
                    melody_string: melody_str,
                    rh_midi: v.rh_midi.clone(),
                    lh_midi: v.lh_midi.clone(),
                    rh_strings: rh_strs,
                    lh_strings: lh_strs,
                    beats,
                    chord_name: if is_chord_change {
                        let all: Vec<i32> = v.rh_midi.iter().chain(v.lh_midi.iter()).copied().collect();
                        Some(identify_chord_names(&all, key_root))
                    } else { None },
                    rh_chord: if is_chord_change { Some(identify_chord_names(&v.rh_midi, key_root)) } else { None },
                    lh_chord: if is_chord_change { Some(identify_chord_names(&v.lh_midi, key_root)) } else { None },
                    is_chord_change,
                    display_rh: vec![],
                    display_lh: vec![],
                });
                last_voicing = voicing;
            } else {
                let mel_str = music::midi_to_harp_string(snapped, key_root)
                    .map(|s| music::to_relative_string(s, key_root))
                    .unwrap_or(1);
                events.push(ScoreEvent::Note {
                    melody_midi: snapped,
                    melody_string: mel_str,
                    rh_midi: vec![snapped],
                    lh_midi: vec![],
                    rh_strings: vec![mel_str],
                    lh_strings: vec![],
                    beats,
                    chord_name: Some(music::midi_to_name(snapped)),
                    rh_chord: None, lh_chord: None,
                    is_chord_change: true,
                    display_rh: vec![mel_str],
                    display_lh: vec![],
                });
            }
        }
        events.push(ScoreEvent::Bar);
    }

    let mut score = Score {
        title: "Random Drill".into(),
        number: "0".into(),
        key: key.into(),
        key_root,
        meter_num: beats_per_measure,
        meter_den: 4,
        tempo: 120,
        events,
    };

    compute_display_strings(&mut score, rh_fingers, lh_fingers, &[]);
    score
}
