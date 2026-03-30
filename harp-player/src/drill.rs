use std::collections::HashMap;
use std::time::Instant;

use serde::{Deserialize, Serialize};

use crate::music::{self, PedalState};

// ── Drill mode ────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DrillMode {
    Intervals,
    Chords,
}

impl DrillMode {
    pub fn label(&self) -> &'static str {
        match self {
            DrillMode::Intervals => "Intervals",
            DrillMode::Chords => "Chords",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoteResult {
    Hit,
    Miss,
    Waiting,
}

// ── Difficulty ────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Difficulty {
    Relaxed,
    Moderate,
    Brisk,
}

impl Difficulty {
    pub fn start_secs(&self) -> f32 {
        match self {
            Difficulty::Relaxed => 5.0,
            Difficulty::Moderate => 3.0,
            Difficulty::Brisk => 2.0,
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

// ── Adaptive pacing constants ─────────────────────────────────────

const STREAK_TO_SPEED_UP: u32 = 4;
const SPEED_UP_FACTOR: f32 = 0.90;
const SLOW_DOWN_FACTOR: f32 = 1.15;
const MIN_SECS: f32 = 0.4;
const MAX_SECS: f32 = 12.0;
const ACCURACY_GATE: f32 = 90.0;
const RANGE_EXPAND_INTERVAL: u32 = 6;
const LOOKAHEAD_COUNT: usize = 8;

pub const DURATION_OPTIONS: [(&str, u32); 3] = [
    ("1 min", 60),
    ("2 min", 120),
    ("5 min", 300),
];

// ── Persistent progress (Mavis Beacon style) ──────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionRecord {
    pub npm: f32,
    pub level: u32,
    pub pace: f32,
    pub challenges: u32,
    pub notes: usize,
}

/// Per-interval/chord accuracy stats
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TypeStats {
    pub attempts: u32,
    pub total_time: f32,  // sum of response times
}

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
    /// Lifetime stats per interval/chord label (e.g. "M3" -> {attempts, total_time})
    #[serde(default)]
    pub type_stats: HashMap<String, TypeStats>,
}

impl Default for Progress {
    fn default() -> Self {
        Self {
            level: 0, target_secs: 5.0, best_npm: 0.0, total_notes: 0, total_sessions: 0,
            best_streak: 0, total_time_secs: 0.0, history: Vec::new(), type_stats: HashMap::new(),
        }
    }
}

fn progress_path() -> std::path::PathBuf {
    if cfg!(target_os = "android") {
        std::path::PathBuf::from("/data/data/com.harp.practice/practice_progress.json")
    } else {
        let home = std::env::var("HOME").unwrap_or_else(|_| ".".into());
        std::path::PathBuf::from(home).join(".practice_progress.json")
    }
}

pub fn load_progress() -> Progress {
    std::fs::read_to_string(progress_path())
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_default()
}

pub fn save_progress(p: &Progress) {
    if let Ok(json) = serde_json::to_string_pretty(p) {
        let _ = std::fs::write(progress_path(), json);
    }
}

// ── Chord templates ───────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct ChordTemplate {
    pub name: &'static str,
    pub short: &'static str,
    pub intervals: &'static [u8],
}

pub const CHORD_TEMPLATES: &[ChordTemplate] = &[
    ChordTemplate { name: "Major",    short: "Maj",  intervals: &[0, 4, 7] },
    ChordTemplate { name: "Minor",    short: "m",    intervals: &[0, 3, 7] },
    ChordTemplate { name: "Dim",      short: "dim",  intervals: &[0, 3, 6] },
    ChordTemplate { name: "Maj7",     short: "M7",   intervals: &[0, 4, 7, 11] },
    ChordTemplate { name: "Min7",     short: "m7",   intervals: &[0, 3, 7, 10] },
    ChordTemplate { name: "Dom7",     short: "7",    intervals: &[0, 4, 7, 10] },
    ChordTemplate { name: "Half-dim", short: "\u{00f8}7", intervals: &[0, 3, 6, 10] },
];

// ── Challenge ─────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct Challenge {
    pub notes: Vec<u8>,
    pub label: String,
    pub next_finger: usize,
}

impl Challenge {
    pub fn current_note(&self) -> Option<u8> {
        self.notes.get(self.next_finger).copied()
    }
    pub fn is_complete(&self) -> bool {
        self.next_finger >= self.notes.len()
    }
}

#[derive(Debug, Clone)]
pub struct TimedChallenge {
    pub challenge: Challenge,
    pub start_time: f32,
    pub duration: f32,
}

// ── Config ────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct DrillConfig {
    pub duration_secs: u32,
    pub mode: DrillMode,
    pub low_midi: u8,
    pub high_midi: u8,
    pub difficulty: Difficulty,
    pub chord_size: usize,
}

/// Max hand span in diatonic steps (10 strings on harp)
const MAX_SPAN_STRINGS: u8 = 10;

impl Default for DrillConfig {
    fn default() -> Self {
        Self {
            duration_secs: 120,
            mode: DrillMode::Intervals,
            low_midi: 33,  // A1
            high_midi: 88, // E6
            difficulty: Difficulty::Relaxed,
            chord_size: 3,
        }
    }
}

// ── Session ───────────────────────────────────────────────────────

pub struct DrillSession {
    pub config: DrillConfig,
    pub timeline: Vec<TimedChallenge>,
    pub current_idx: usize,
    pub challenges_completed: u32,
    pub hits: usize,
    pub misses: usize,
    pub start_time: Instant,
    pub finished: bool,

    pub target_secs: f32,
    pub streak: u32,
    pub level: u32,

    pub effective_low: u8,
    pub effective_high: u8,
    pub weak_times: HashMap<String, Vec<f32>>,

    available: Vec<u8>,
    pedals: PedalState,
    next_start: f32,
}

impl DrillSession {
    pub fn new(config: DrillConfig, pedals: &PedalState, progress: &Progress) -> Self {
        let now = Instant::now();

        // Resume from saved progress: use saved pace, or difficulty default if first session
        let target_secs = if progress.total_sessions > 0 {
            progress.target_secs
        } else {
            config.difficulty.start_secs()
        };

        let mid = (config.low_midi as u16 + config.high_midi as u16) / 2;
        let half_span = 7u16;
        let effective_low = (mid.saturating_sub(half_span) as u8).max(config.low_midi);
        let effective_high = ((mid + half_span) as u8).min(config.high_midi);

        let available = music::available_notes_in_range(pedals, effective_low, effective_high);

        let mut session = Self {
            config,
            timeline: Vec::new(),
            current_idx: 0,
            challenges_completed: 0,
            hits: 0,
            misses: 0,
            start_time: now,
            finished: false,
            target_secs,
            streak: 0,
            level: progress.level,
            effective_low,
            effective_high,
            weak_times: HashMap::new(),
            available,
            pedals: pedals.clone(),
            next_start: 0.0,
        };
        session.fill_lookahead();
        session
    }

    pub fn elapsed(&self) -> f32 {
        self.start_time.elapsed().as_secs_f32()
    }

    pub fn current_challenge(&self) -> Option<&TimedChallenge> {
        self.timeline.get(self.current_idx)
    }

    pub fn current_note(&self) -> Option<u8> {
        self.current_challenge().and_then(|tc| tc.challenge.current_note())
    }

    pub fn note_elapsed_fraction(&self) -> f32 {
        if let Some(tc) = self.current_challenge() {
            (self.elapsed() - tc.start_time) / tc.duration
        } else { 0.0 }
    }

    pub fn current_bpm(&self) -> f32 { 60.0 / self.target_secs }

    pub fn time_remaining(&self) -> f32 {
        (self.config.duration_secs as f32 - self.elapsed()).max(0.0)
    }

    fn check_time(&mut self) {
        if self.elapsed() >= self.config.duration_secs as f32 {
            self.finished = true;
        }
    }

    fn fill_lookahead(&mut self) {
        while self.timeline.len() < self.current_idx + LOOKAHEAD_COUNT {
            let weak = self.weakest_labels(3);
            let ch = generate_challenge(&self.config, &self.available, &weak);
            self.timeline.push(TimedChallenge {
                challenge: ch,
                start_time: self.next_start,
                duration: self.target_secs,
            });
            self.next_start += self.target_secs;
        }
    }

    fn complete_current(&mut self, on_pace: bool) {
        if let Some(tc) = self.timeline.get(self.current_idx) {
            let elapsed = self.elapsed() - tc.start_time;
            self.weak_times.entry(tc.challenge.label.clone()).or_default().push(elapsed);
        }

        self.challenges_completed += 1;

        if on_pace && self.accuracy() >= ACCURACY_GATE {
            self.streak += 1;
            if self.streak >= STREAK_TO_SPEED_UP {
                self.target_secs = (self.target_secs * SPEED_UP_FACTOR).max(MIN_SECS);
                self.streak = 0;
                self.level += 1;
            }
        } else if !on_pace {
            self.target_secs = (self.target_secs * SLOW_DOWN_FACTOR).min(MAX_SECS);
            self.streak = 0;
        }

        self.maybe_expand_range();
        self.current_idx += 1;
        self.fill_lookahead();
        self.check_time();
    }

    fn advance_note(&mut self, on_pace: bool) {
        self.hits += 1;
        if let Some(tc) = self.timeline.get_mut(self.current_idx) {
            tc.challenge.next_finger += 1;
            if tc.challenge.is_complete() {
                self.complete_current(on_pace);
            }
        }
    }

    fn maybe_expand_range(&mut self) {
        if self.challenges_completed % RANGE_EXPAND_INTERVAL != 0 || self.challenges_completed == 0 {
            return;
        }
        let pcs = self.pedals.available_pitch_classes();
        if self.effective_low > self.config.low_midi {
            for midi in (self.config.low_midi..self.effective_low).rev() {
                if pcs.contains(&(midi % 12)) { self.effective_low = midi; break; }
            }
        }
        if self.effective_high < self.config.high_midi {
            for midi in (self.effective_high + 1)..=self.config.high_midi {
                if pcs.contains(&(midi % 12)) { self.effective_high = midi; break; }
            }
        }
        self.available = music::available_notes_in_range(&self.pedals, self.effective_low, self.effective_high);
    }

    fn weakest_labels(&self, n: usize) -> Vec<String> {
        let mut avgs: Vec<(String, f32)> = self.weak_times.iter()
            .filter(|(_, t)| !t.is_empty())
            .map(|(l, t)| (l.clone(), t.iter().sum::<f32>() / t.len() as f32))
            .collect();
        avgs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        avgs.into_iter().take(n).map(|(l, _)| l).collect()
    }

    pub fn check_note(&mut self, detected_midi: u8, score_mistakes: bool) -> NoteResult {
        if self.finished { return NoteResult::Waiting; }
        let Some(target) = self.current_note() else { return NoteResult::Waiting; };
        if detected_midi == target {
            if let Some(tc) = self.timeline.get(self.current_idx) {
                let on_pace = (self.elapsed() - tc.start_time) <= tc.duration * 1.1;
                self.advance_note(on_pace);
            }
            NoteResult::Hit
        } else {
            if score_mistakes {
                self.misses += 1;
                self.streak = 0;
            }
            NoteResult::Miss
        }
    }

    pub fn advance_manual(&mut self) {
        if self.finished { return; }
        if let Some(tc) = self.timeline.get(self.current_idx) {
            let since = self.elapsed() - tc.start_time;
            let on_pace = since <= tc.duration * 1.1;
            self.advance_note(on_pace);
        }
    }

    pub fn check_auto_advance(&mut self) -> bool {
        if self.finished { return false; }
        let Some(tc) = self.timeline.get(self.current_idx) else { return false; };
        if self.elapsed() >= tc.start_time + tc.duration {
            self.streak = 0;
            self.target_secs = (self.target_secs * SLOW_DOWN_FACTOR).min(MAX_SECS);
            self.complete_current(false);
            return true;
        }
        false
    }

    pub fn notes_per_minute(&self) -> f32 {
        let e = self.elapsed();
        if e < 0.1 { 0.0 } else { self.hits as f32 / (e / 60.0) }
    }

    pub fn accuracy(&self) -> f32 {
        let total = self.hits + self.misses;
        if total == 0 { 100.0 } else { self.hits as f32 / total as f32 * 100.0 }
    }

    pub fn is_finished(&self) -> bool { self.finished }

    pub fn progress_text(&self) -> String {
        let r = self.time_remaining();
        format!("{}:{:02}", r as u32 / 60, r as u32 % 60)
    }

    pub fn weakest_labels_public(&self, n: usize) -> Vec<(String, f32)> {
        let mut avgs: Vec<(String, f32)> = self.weak_times.iter()
            .filter(|(_, t)| !t.is_empty())
            .map(|(l, t)| (l.clone(), t.iter().sum::<f32>() / t.len() as f32))
            .collect();
        avgs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        avgs.into_iter().take(n).collect()
    }

    pub fn range_display(&self) -> String {
        let ln = music::NOTE_NAMES[self.effective_low as usize % 12];
        let lo = (self.effective_low as i32 / 12) - 1;
        let hn = music::NOTE_NAMES[self.effective_high as usize % 12];
        let ho = (self.effective_high as i32 / 12) - 1;
        format!("{ln}{lo}\u{2013}{hn}{ho}")
    }

    /// Save progress at end of session.
    pub fn save_to_progress(&self) -> Progress {
        let npm = self.notes_per_minute();
        let mut p = load_progress();
        p.level = self.level;
        p.target_secs = self.target_secs;
        p.total_notes += self.hits as u64;
        p.total_sessions += 1;
        p.total_time_secs += self.elapsed();
        if npm > p.best_npm { p.best_npm = npm; }
        if self.streak > p.best_streak { p.best_streak = self.streak; }
        // Merge per-type stats
        for (label, times) in &self.weak_times {
            let entry = p.type_stats.entry(label.clone()).or_default();
            entry.attempts += times.len() as u32;
            entry.total_time += times.iter().sum::<f32>();
        }
        p.history.push(SessionRecord {
            npm, level: self.level, pace: self.target_secs,
            challenges: self.challenges_completed, notes: self.hits,
        });
        if p.history.len() > 50 { p.history.drain(..p.history.len() - 50); }
        save_progress(&p);
        p
    }
}

// ── Challenge generation ──────────────────────────────────────────

fn generate_challenge(config: &DrillConfig, available: &[u8], weak_labels: &[String]) -> Challenge {
    let mut rng = rand::thread_rng();
    match config.mode {
        DrillMode::Intervals => generate_interval(&mut rng, available, weak_labels),
        DrillMode::Chords => generate_chord(&mut rng, available, config.chord_size, weak_labels),
    }
}

/// Check if all notes fit within MAX_SPAN_STRINGS diatonic steps.
fn within_span(notes: &[u8], available: &[u8]) -> bool {
    if notes.len() < 2 { return true; }
    let low = *notes.iter().min().unwrap();
    let high = *notes.iter().max().unwrap();
    // Count how many available (diatonic) notes lie between low and high inclusive
    let span = available.iter().filter(|&&m| m >= low && m <= high).count();
    span <= MAX_SPAN_STRINGS as usize
}

fn generate_interval(rng: &mut impl rand::Rng, available: &[u8], weak_labels: &[String]) -> Challenge {
    if available.len() < 2 {
        return Challenge { notes: vec![60, 64], label: "M3".into(), next_finger: 0 };
    }
    for _ in 0..20 {
        let i = rng.gen_range(0..available.len());
        let j = rng.gen_range(0..available.len());
        if i == j { continue; }
        let low = available[i].min(available[j]);
        let high = available[i].max(available[j]);
        let notes = vec![low, high];
        if !within_span(&notes, available) { continue; }
        let semi = high - low;
        if semi > 17 { continue; } // ~10 strings max
        let label = music::interval_name(semi).to_string();
        if !weak_labels.is_empty() && !weak_labels.contains(&label) && rng.gen_ratio(1, 3) { continue; }
        return Challenge { notes, label, next_finger: 0 };
    }
    let (low, high) = (available[0], available[available.len().min(5)]);
    Challenge { notes: vec![low, high], label: music::interval_name(high - low).into(), next_finger: 0 }
}

fn generate_chord(rng: &mut impl rand::Rng, available: &[u8], chord_size: usize, weak_labels: &[String]) -> Challenge {
    if available.len() < 3 {
        return Challenge { notes: vec![60, 64, 67], label: "C Maj".into(), next_finger: 0 };
    }
    let templates: Vec<&ChordTemplate> = CHORD_TEMPLATES.iter()
        .filter(|t| t.intervals.len() == chord_size || (chord_size == 3 && t.intervals.len() == 3))
        .collect();
    let templates = if templates.is_empty() {
        CHORD_TEMPLATES.iter().filter(|t| t.intervals.len() == 3).collect::<Vec<_>>()
    } else { templates };

    for _ in 0..30 {
        let root = available[rng.gen_range(0..available.len())];
        let tmpl = templates[rng.gen_range(0..templates.len())];
        let mut notes = Vec::new();
        let mut ok = true;
        for &iv in tmpl.intervals {
            let m = root + iv;
            if !available.contains(&m) { ok = false; break; }
            notes.push(m);
        }
        if !ok { continue; }
        notes.sort();
        if !within_span(&notes, available) { continue; }
        let label = format!("{} {}", music::NOTE_NAMES[root as usize % 12], tmpl.short);
        if !weak_labels.is_empty() && !weak_labels.contains(&label) && rng.gen_ratio(1, 3) { continue; }
        return Challenge { notes, label, next_finger: 0 };
    }
    let n = chord_size.min(available.len());
    Challenge { notes: available[..n].to_vec(), label: "?".into(), next_finger: 0 }
}

// ── Generate ScoreEvents from drill timeline ──────────────────────

/// Convert drill challenges into ScoreEvents so render_score() can draw them
/// identically to hymns — with chord names, string numbers, finger assignments.
pub fn timeline_to_events(
    timeline: &[TimedChallenge],
    key_root: i32,
    rh_fingers: usize,
    lh_fingers: usize,
) -> Vec<crate::abc::ScoreEvent> {
    use crate::abc::ScoreEvent;
    use crate::voicing::voice_from_satb_fingers;
    use crate::chord::identify_chord_names;
    use crate::music::{snap_to_diatonic, to_relative_string};

    let mut events = Vec::new();
    let mut prev_voicing: Option<crate::voicing::Voicing> = None;

    for tc in timeline {
        let ch = &tc.challenge;
        let beats = tc.duration; // match scroll speed to drill timing

        // Sort notes: highest = soprano, next = alto, etc.
        let mut sorted = ch.notes.clone();
        sorted.sort();
        sorted.reverse(); // highest first

        let soprano = snap_to_diatonic(sorted[0] as i32, key_root);
        let alto = if sorted.len() > 1 { Some(snap_to_diatonic(sorted[1] as i32, key_root)) } else { None };
        let tenor = if sorted.len() > 2 { Some(snap_to_diatonic(sorted[2] as i32, key_root)) } else { None };
        let bass = if sorted.len() > 3 { Some(snap_to_diatonic(sorted[3] as i32, key_root)) } else if sorted.len() > 1 {
            Some(snap_to_diatonic(*sorted.last().unwrap() as i32, key_root))
        } else { None };

        let voicing = voice_from_satb_fingers(
            soprano, alto, tenor, bass, key_root,
            prev_voicing.as_ref(),
            rh_fingers, lh_fingers,
        );

        if let Some(ref v) = voicing {
            let rh_strs: Vec<i32> = v.rh_strings.iter()
                .map(|&s| to_relative_string(s, key_root)).collect();
            let lh_strs: Vec<i32> = v.lh_strings.iter()
                .map(|&s| to_relative_string(s, key_root)).collect();
            let melody_str = *rh_strs.first().unwrap_or(&1);

            let all_notes: Vec<i32> = v.rh_midi.iter().chain(v.lh_midi.iter()).copied().collect();
            let chord_name = identify_chord_names(&all_notes, key_root);
            let rh_chord = identify_chord_names(&v.rh_midi, key_root);
            let lh_chord = identify_chord_names(&v.lh_midi, key_root);

            events.push(ScoreEvent::Note {
                melody_midi: soprano,
                melody_string: melody_str,
                rh_midi: v.rh_midi.clone(),
                lh_midi: v.lh_midi.clone(),
                rh_strings: rh_strs.clone(),
                lh_strings: lh_strs.clone(),
                beats,
                chord_name: Some(chord_name),
                rh_chord: Some(rh_chord),
                lh_chord: Some(lh_chord),
                is_chord_change: true,
                display_rh: rh_strs,
                display_lh: lh_strs,
            });
            prev_voicing = voicing;
        } else {
            // Fallback: show melody note with chord name
            let mel_str = crate::music::midi_to_harp_string(soprano, key_root)
                .map(|s| to_relative_string(s, key_root))
                .unwrap_or(1);
            events.push(ScoreEvent::Note {
                melody_midi: soprano,
                melody_string: mel_str,
                rh_midi: vec![soprano],
                lh_midi: vec![],
                rh_strings: vec![mel_str],
                lh_strings: vec![],
                beats,
                chord_name: Some(music::midi_to_name(soprano)),
                rh_chord: None,
                lh_chord: None,
                is_chord_change: true,
                display_rh: vec![mel_str],
                display_lh: vec![],
            });
        }

        // Bar after each challenge
        events.push(ScoreEvent::Bar);
    }

    events
}
