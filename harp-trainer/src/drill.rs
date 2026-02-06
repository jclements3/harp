use std::collections::HashMap;
use std::time::Instant;

use crate::music::{self, HandPosition, PedalState};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputMode {
    Mic,
    Manual,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoteResult {
    Hit,
    Miss,
    Waiting,
}

// ── Adaptive pacing ────────────────────────────────────────────────

/// Starting difficulty presets: seconds per note
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Difficulty {
    Relaxed,     // 4.0s per note
    Moderate,    // 2.5s per note
    Brisk,       // 1.5s per note
}

impl Difficulty {
    pub fn start_secs(&self) -> f32 {
        match self {
            Difficulty::Relaxed => 4.0,
            Difficulty::Moderate => 2.5,
            Difficulty::Brisk => 1.5,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Difficulty::Relaxed => "Relaxed",
            Difficulty::Moderate => "Moderate",
            Difficulty::Brisk => "Brisk",
        }
    }
}

const STREAK_TO_SPEED_UP: u32 = 4;  // consecutive on-pace hits to speed up
const SPEED_UP_FACTOR: f32 = 0.90;  // reduce target time by 10%
const SLOW_DOWN_FACTOR: f32 = 1.15; // increase target time by 15%
const MIN_SECS_PER_NOTE: f32 = 0.4; // fastest possible pace
const MAX_SECS_PER_NOTE: f32 = 8.0; // slowest possible pace
const ACCURACY_GATE: f32 = 90.0;    // require 90%+ accuracy to speed up
const RANGE_EXPAND_INTERVAL: u32 = 8; // expand range every N positions

/// Session duration options
pub const DURATION_OPTIONS: [(&str, u32); 3] = [
    ("1 min", 60),
    ("2 min", 120),
    ("5 min", 300),
];

// ── Config ─────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct DrillConfig {
    pub duration_secs: u32,
    pub fingers: usize,
    pub low_midi: u8,
    pub high_midi: u8,
    pub difficulty: Difficulty,
    pub input_mode: InputMode,
}

impl Default for DrillConfig {
    fn default() -> Self {
        Self {
            duration_secs: 120,
            fingers: 4,
            low_midi: 48,  // C3
            high_midi: 79, // G5
            difficulty: Difficulty::Relaxed,
            input_mode: InputMode::Manual,
        }
    }
}

// ── Session ────────────────────────────────────────────────────────

pub struct DrillSession {
    pub config: DrillConfig,

    // Hand position state
    pub position: HandPosition,
    pub finger_index: usize,
    pub positions_completed: u32,

    // Scoring
    pub hits: usize,
    pub misses: usize,

    // Timing
    pub start_time: Instant,
    pub note_start_time: Instant,
    pub finished: bool,

    // Adaptive pacing
    pub target_secs: f32,
    pub streak: u32,
    pub level: u32,
    pub last_was_late: bool,

    // Metronome
    last_click_time: Instant,

    // Progressive range
    pub effective_low: u8,
    pub effective_high: u8,

    // Weak note tracking: midi -> list of response times in seconds
    pub note_times: HashMap<u8, Vec<f32>>,

    // Cached available notes for current effective range
    available: Vec<u8>,
    pedals: PedalState,
}

impl DrillSession {
    pub fn new(config: DrillConfig, pedals: &PedalState) -> Self {
        let now = Instant::now();
        let target_secs = config.difficulty.start_secs();

        // Start with a narrow range (~1 octave) centered in the configured range
        let mid = (config.low_midi as u16 + config.high_midi as u16) / 2;
        let half_span = 7u16;
        let effective_low = (mid.saturating_sub(half_span) as u8).max(config.low_midi);
        let effective_high = ((mid + half_span) as u8).min(config.high_midi);

        let available = music::available_notes_in_range(pedals, effective_low, effective_high);
        let fingers = config.fingers.min(available.len().max(1));

        let position = music::generate_hand_position(&available, fingers, &[], None)
            .unwrap_or(HandPosition { notes: vec![60; fingers] });

        Self {
            config,
            position,
            finger_index: 0,
            positions_completed: 0,
            hits: 0,
            misses: 0,
            start_time: now,
            note_start_time: now,
            finished: false,
            target_secs,
            streak: 0,
            level: 0,
            last_was_late: false,
            last_click_time: now,
            effective_low,
            effective_high,
            note_times: HashMap::new(),
            available,
            pedals: pedals.clone(),
        }
    }

    pub fn current_note(&self) -> Option<u8> {
        self.position.notes.get(self.finger_index).copied()
    }

    /// How far through the target time we are (0.0 = just shown, 1.0 = at target, >1.0 = late)
    pub fn note_elapsed_fraction(&self) -> f32 {
        self.note_start_time.elapsed().as_secs_f32() / self.target_secs
    }

    /// Current BPM derived from adaptive target pace
    pub fn current_bpm(&self) -> f32 {
        60.0 / self.target_secs
    }

    /// Seconds remaining in the session
    pub fn time_remaining(&self) -> f32 {
        (self.config.duration_secs as f32 - self.start_time.elapsed().as_secs_f32()).max(0.0)
    }

    fn check_time(&mut self) {
        if self.start_time.elapsed().as_secs_f32() >= self.config.duration_secs as f32 {
            self.finished = true;
        }
    }

    fn record_response(&mut self, midi: u8) {
        let elapsed = self.note_start_time.elapsed().as_secs_f32();
        self.note_times.entry(midi).or_default().push(elapsed);
    }

    fn advance_note(&mut self) {
        let elapsed = self.note_start_time.elapsed().as_secs_f32();
        let on_pace = elapsed <= self.target_secs * 1.1; // 10% grace

        if let Some(midi) = self.current_note() {
            self.record_response(midi);
        }

        self.hits += 1;
        self.finger_index += 1;
        self.last_was_late = !on_pace;

        // Adaptive pacing with accuracy gate
        if on_pace && self.accuracy() >= ACCURACY_GATE {
            self.streak += 1;
            if self.streak >= STREAK_TO_SPEED_UP {
                self.target_secs = (self.target_secs * SPEED_UP_FACTOR).max(MIN_SECS_PER_NOTE);
                self.streak = 0;
                self.level += 1;
            }
        } else if !on_pace {
            self.target_secs = (self.target_secs * SLOW_DOWN_FACTOR).min(MAX_SECS_PER_NOTE);
            self.streak = 0;
        }

        self.note_start_time = Instant::now();

        // Check if we finished the current hand position
        if self.finger_index >= self.position.notes.len() {
            self.positions_completed += 1;
            self.maybe_expand_range();
            self.next_position();
        }

        self.check_time();
    }

    fn next_position(&mut self) {
        let weak = self.weakest_notes(3);
        let last = Some(&self.position);
        let fingers = self.config.fingers.min(self.available.len().max(1));
        self.position = music::generate_hand_position(&self.available, fingers, &weak, last)
            .unwrap_or(HandPosition { notes: vec![60; fingers] });
        self.finger_index = 0;
    }

    fn maybe_expand_range(&mut self) {
        if self.positions_completed % RANGE_EXPAND_INTERVAL != 0 || self.positions_completed == 0 {
            return;
        }

        let pcs = self.pedals.available_pitch_classes();

        // Expand one note lower
        if self.effective_low > self.config.low_midi {
            for midi in (self.config.low_midi..self.effective_low).rev() {
                if pcs.contains(&(midi % 12)) {
                    self.effective_low = midi;
                    break;
                }
            }
        }

        // Expand one note higher
        if self.effective_high < self.config.high_midi {
            for midi in (self.effective_high + 1)..=self.config.high_midi {
                if pcs.contains(&(midi % 12)) {
                    self.effective_high = midi;
                    break;
                }
            }
        }

        // Rebuild available notes cache
        self.available = music::available_notes_in_range(&self.pedals, self.effective_low, self.effective_high);
    }

    /// Get the N notes with slowest average response times
    fn weakest_notes(&self, n: usize) -> Vec<u8> {
        let mut avgs: Vec<(u8, f32)> = self.note_times.iter()
            .filter(|(_, times)| !times.is_empty())
            .map(|(&midi, times)| {
                let avg = times.iter().sum::<f32>() / times.len() as f32;
                (midi, avg)
            })
            .collect();
        avgs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        avgs.into_iter().take(n).map(|(midi, _)| midi).collect()
    }

    /// Check a detected MIDI note against the current target (mic mode).
    pub fn check_note(&mut self, detected_midi: u8) -> NoteResult {
        if self.finished {
            return NoteResult::Waiting;
        }
        let Some(target) = self.current_note() else {
            return NoteResult::Waiting;
        };

        if detected_midi == target {
            self.advance_note();
            NoteResult::Hit
        } else {
            self.misses += 1;
            // On miss, slow down and reset streak
            self.target_secs = (self.target_secs * SLOW_DOWN_FACTOR).min(MAX_SECS_PER_NOTE);
            self.streak = 0;
            NoteResult::Miss
        }
    }

    /// Advance to the next note (manual mode).
    pub fn advance_manual(&mut self) {
        if self.finished {
            return;
        }
        self.advance_note();
    }

    /// Auto-advance when the beat arrives (music scrolls at tempo).
    /// Returns true if the note was auto-advanced (missed the beat).
    pub fn check_auto_advance(&mut self) -> bool {
        if self.finished {
            return false;
        }
        // Auto-advance at the beat (target_secs) — the music keeps moving
        let elapsed = self.note_start_time.elapsed().as_secs_f32();
        if elapsed >= self.target_secs {
            // Count as miss in mic mode
            if self.config.input_mode == InputMode::Mic {
                self.misses += 1;
            }
            self.last_was_late = true;
            self.streak = 0;
            self.target_secs = (self.target_secs * SLOW_DOWN_FACTOR).min(MAX_SECS_PER_NOTE);
            self.finger_index += 1;
            self.note_start_time = Instant::now();

            // Check if we finished the current hand position
            if self.finger_index >= self.position.notes.len() {
                self.positions_completed += 1;
                self.maybe_expand_range();
                self.next_position();
            }

            self.check_time();
            return true;
        }
        false
    }

    /// Should the metronome click now? Ticks at the current adaptive pace.
    pub fn should_click(&mut self) -> bool {
        let elapsed = self.last_click_time.elapsed().as_secs_f32();
        if elapsed >= self.target_secs {
            self.last_click_time = Instant::now();
            true
        } else {
            false
        }
    }

    pub fn notes_per_minute(&self) -> f32 {
        let elapsed = self.start_time.elapsed().as_secs_f32();
        if elapsed < 0.1 {
            return 0.0;
        }
        self.hits as f32 / (elapsed / 60.0)
    }

    pub fn accuracy(&self) -> f32 {
        let total = self.hits + self.misses;
        if total == 0 {
            return 100.0;
        }
        self.hits as f32 / total as f32 * 100.0
    }

    pub fn is_finished(&self) -> bool {
        self.finished
    }

    pub fn progress_text(&self) -> String {
        let remaining = self.time_remaining();
        let mins = remaining as u32 / 60;
        let secs = remaining as u32 % 60;
        format!("{mins}:{secs:02}")
    }

    pub fn range_display(&self) -> String {
        let low_name = music::NOTE_NAMES[self.effective_low as usize % 12];
        let low_oct = (self.effective_low as i32 / 12) - 1;
        let high_name = music::NOTE_NAMES[self.effective_high as usize % 12];
        let high_oct = (self.effective_high as i32 / 12) - 1;
        format!("{low_name}{low_oct}\u{2013}{high_name}{high_oct}")
    }
}
