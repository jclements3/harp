/// ABC parser that produces harp-voiced notation events with durations.
///
/// Uses the same voicing engine as harp-hymnal to produce RH and LH MIDI notes,
/// then wraps them with duration info for the notation renderer.

use crate::music::*;
use crate::chord::identify_chord_names;
use crate::voicing::{voice_from_satb, voice_from_satb_fingers, Voicing};

#[derive(Debug, Clone)]
pub struct Score {
    pub title: String,
    pub number: String,
    pub key: String,
    pub key_root: i32,
    pub meter_num: u8,
    pub meter_den: u8,
    pub tempo: u16,
    pub events: Vec<ScoreEvent>,
}

#[derive(Debug, Clone)]
pub enum ScoreEvent {
    /// A beat with harp voicing: melody + RH notes on treble, LH notes on bass
    Note {
        melody_midi: i32,
        melody_string: i32,
        rh_midi: Vec<i32>,
        lh_midi: Vec<i32>,
        rh_strings: Vec<i32>,
        lh_strings: Vec<i32>,
        beats: f32,
        chord_name: Option<String>,
        rh_chord: Option<String>,
        lh_chord: Option<String>,
        is_chord_change: bool,
        // Pre-computed display strings (selected for finger count + max span)
        display_rh: Vec<i32>,
        display_lh: Vec<i32>,
    },
    Rest { beats: f32 },
    Bar,
}

// ── ABC parsing internals ──

// Correct ABC base: uppercase C = C4 = MIDI 60 (standard ABC notation).
const ABC_BASE: [(char, i32); 7] = [
    ('C', 60), ('D', 62), ('E', 64), ('F', 65), ('G', 67), ('A', 69), ('B', 71),
];

/// No display offset needed for notation — MIDI values are correct.
pub const DISPLAY_OCTAVE_OFFSET: i32 = 0;

/// String number offset: ABC base 60 produces strings 7 higher than the
/// hymnal's base-48 system. Subtract this from relative string numbers
/// to match the hymnal's display.
pub const STRING_NUMBER_OFFSET: i32 = 7;

fn abc_note_to_midi(token: &str, key_sig_acc: &[i32; 7]) -> Option<i32> {
    let bytes = token.as_bytes();
    let mut acc: Option<i32> = None;
    let mut i = 0;

    while i < bytes.len() && matches!(bytes[i], b'^' | b'_' | b'=') {
        if acc.is_none() { acc = Some(0); }
        match bytes[i] {
            b'^' => *acc.as_mut().unwrap() += 1,
            b'_' => *acc.as_mut().unwrap() -= 1,
            b'=' => acc = Some(0),
            _ => {}
        }
        i += 1;
    }

    if i >= bytes.len() { return None; }
    let ch = bytes[i] as char;
    let is_lower = ch.is_ascii_lowercase();
    let letter = ch.to_ascii_uppercase();

    let base = ABC_BASE.iter().find(|(l, _)| *l == letter)?.1;
    let letter_idx = match letter {
        'C' => 0, 'D' => 1, 'E' => 2, 'F' => 3, 'G' => 4, 'A' => 5, 'B' => 6,
        _ => return None,
    };
    let adjustment = acc.unwrap_or(key_sig_acc[letter_idx]);
    let mut midi = base + adjustment;
    if is_lower { midi += 12; }
    i += 1;

    while i < bytes.len() {
        match bytes[i] {
            b',' => midi -= 12,
            b'\'' => midi += 12,
            _ => break,
        }
        i += 1;
    }

    Some(midi)
}

fn parse_duration(s: &str) -> f32 {
    if s.is_empty() { return 1.0; }
    if let Some(rest) = s.strip_prefix('/') {
        let den: f32 = rest.parse().unwrap_or(2.0);
        return 1.0 / den;
    }
    if s.contains('/') {
        let parts: Vec<&str> = s.splitn(2, '/').collect();
        let num: f32 = parts[0].parse().unwrap_or(1.0);
        let den: f32 = parts[1].parse().unwrap_or(1.0);
        return num / den;
    }
    s.parse().unwrap_or(1.0)
}

#[derive(Debug)]
enum Token {
    Bar,
    Rest(f32),
    Note(i32, f32), // midi, duration_multiplier
}

fn parse_voice_tokens(music: &str, key_sig_acc: &[i32; 7]) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut chars = music.chars().peekable();
    let mut note_buf = String::new();
    let mut dur_buf = String::new();

    while let Some(&ch) = chars.peek() {
        if ch == '|' {
            tokens.push(Token::Bar);
            chars.next();
            while chars.peek().map_or(false, |&c| matches!(c, '|' | ']' | '[' | ':')) {
                chars.next();
            }
        } else if ch == 'z' || ch == 'Z' {
            chars.next();
            dur_buf.clear();
            while chars.peek().map_or(false, |c| c.is_ascii_digit() || *c == '/') {
                dur_buf.push(chars.next().unwrap());
            }
            tokens.push(Token::Rest(parse_duration(&dur_buf)));
        } else if matches!(ch, '^' | '_' | '=') || ch.is_ascii_alphabetic() {
            note_buf.clear();
            dur_buf.clear();
            while chars.peek().map_or(false, |&c| matches!(c, '^' | '_' | '=')) {
                note_buf.push(chars.next().unwrap());
            }
            if let Some(&c) = chars.peek() {
                if c.is_ascii_alphabetic() && !"zZxXwWhH".contains(c) {
                    note_buf.push(chars.next().unwrap());
                    while chars.peek().map_or(false, |&c| c == ',' || c == '\'') {
                        note_buf.push(chars.next().unwrap());
                    }
                    while chars.peek().map_or(false, |c| c.is_ascii_digit() || *c == '/') {
                        dur_buf.push(chars.next().unwrap());
                    }
                    if let Some(midi) = abc_note_to_midi(&note_buf, key_sig_acc) {
                        tokens.push(Token::Note(midi, parse_duration(&dur_buf)));
                    }
                } else {
                    chars.next();
                }
            }
        } else if ch == '"' {
            chars.next();
            while chars.peek().map_or(false, |&c| c != '"') { chars.next(); }
            chars.next();
        } else if ch == 'w' && chars.clone().nth(1) == Some(':') {
            break;
        } else {
            chars.next();
        }
    }

    tokens
}

/// Embedded ABC data
const EMBEDDED_ABC_BYTES: &[u8] = include_bytes!("../../data/OpenHymnal.abc");

/// Select display strings from full voicing for given finger counts.
/// RH: melody on thumb, pick chord tones below to maximize span.
/// LH: start below RH, pick chord tones to maximize span, avoid melody degree.
fn select_display_strings(
    melody_string: i32,
    rh_strings: &[i32],
    lh_strings: &[i32],
    rh_fingers: usize,
    lh_fingers: usize,
) -> (Vec<i32>, Vec<i32>) {
    let to_degree = |s: i32| -> i32 {
        let lin = if s > 0 { s - 1 } else { s };
        ((lin % 7) + 7) % 7
    };
    let mel = melody_string;
    let mel_degree = to_degree(mel);

    // All available chord degrees
    let mut all_pool: Vec<i32> = rh_strings.iter().chain(lh_strings.iter()).copied().collect();
    all_pool.sort_by(|a, b| b.cmp(a));
    all_pool.dedup();
    let mut chord_degrees: Vec<i32> = all_pool.iter().map(|&s| to_degree(s)).collect();
    chord_degrees.sort(); chord_degrees.dedup();

    let step_down = |s: i32| -> i32 { let n = s - 1; if n == 0 { -1 } else { n } };
    let chord_tones_below = |start: i32, span: i32| -> Vec<i32> {
        let low = start - span;
        let mut tones = Vec::new();
        let mut s = start;
        while s >= low {
            if chord_degrees.contains(&to_degree(s)) { tones.push(s); }
            s = step_down(s);
        }
        tones
    };

    let rh_max = rh_fingers.min(4);
    let lh_max = lh_fingers.min(4);

    // RH: melody + chord tones below within 10-string span
    let rh_candidates = chord_tones_below(mel, 10);
    let rh_used = if rh_max >= rh_candidates.len() {
        rh_candidates.clone()
    } else {
        let mut pick = vec![rh_candidates[0]]; // melody
        if rh_max >= 2 { pick.push(*rh_candidates.last().unwrap()); }
        let middle: Vec<i32> = rh_candidates[1..rh_candidates.len()-1].to_vec();
        for &c in &middle {
            if pick.len() >= rh_max { break; }
            if !pick.contains(&c) { pick.push(c); }
        }
        pick.sort_by(|a, b| b.cmp(a));
        pick.truncate(rh_max);
        pick
    };

    // LH: below RH, maximize span, prefer non-melody-degree
    let lh_start = *rh_used.last().unwrap_or(&mel) - 1;
    let lh_start = if lh_start == 0 { -1 } else { lh_start };
    let lh_all = chord_tones_below(lh_start, 10);
    let lh_no_mel: Vec<i32> = lh_all.iter().filter(|&&s| to_degree(s) != mel_degree).copied().collect();
    let lh_source = if lh_no_mel.len() >= lh_max { &lh_no_mel } else { &lh_all };
    let lh_used = if lh_max >= lh_source.len() {
        lh_source.clone()
    } else {
        let mut pick = vec![lh_source[0]];
        if lh_max >= 2 { pick.push(*lh_source.last().unwrap()); }
        let middle: Vec<i32> = lh_source[1..lh_source.len()-1].to_vec();
        for &c in &middle {
            if pick.len() >= lh_max { break; }
            if !pick.contains(&c) { pick.push(c); }
        }
        pick.sort_by(|a, b| b.cmp(a));
        pick.truncate(lh_max);
        pick
    };

    (rh_used, lh_used)
}

/// Compute display_rh/display_lh for all events in a score.
/// If random_pattern is non-empty, use per-chord finger counts.
pub fn compute_display_strings(score: &mut Score, rh_fingers: usize, lh_fingers: usize, random_pattern: &[(usize, usize)]) {
    let mut chord_idx = 0usize;
    for event in &mut score.events {
        if let ScoreEvent::Note { melody_string, rh_strings, lh_strings, is_chord_change, display_rh, display_lh, .. } = event {
            if *is_chord_change && !rh_strings.is_empty() {
                let (rf, lf) = if !random_pattern.is_empty() && chord_idx < random_pattern.len() {
                    random_pattern[chord_idx]
                } else {
                    (rh_fingers, lh_fingers)
                };
                let (rh, lh) = select_display_strings(*melody_string, rh_strings, lh_strings, rf, lf);
                *display_rh = rh;
                *display_lh = lh;
                chord_idx += 1;
            } else {
                // Non-chord-change: just melody
                *display_rh = vec![*melody_string];
                *display_lh = vec![];
            }
        }
    }
}

pub fn load_embedded_scores() -> Vec<Score> {
    let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
    let mut scores = parse_all_with_fingers(&text, 4, 4);
    for s in &mut scores {
        compute_display_strings(s, 4, 4, &[]);
    }
    scores
}

pub fn load_embedded_scores_fingers(rh: usize, lh: usize) -> Vec<Score> {
    let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
    parse_all_with_fingers(&text, rh, lh)
}

pub fn revoice_score(score: &Score, rh: usize, lh: usize) -> Score {
    // Re-parse just this one hymn from the embedded ABC
    let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
    let all = parse_all_with_fingers(&text, rh, lh);
    all.into_iter()
        .find(|s| s.number == score.number)
        .unwrap_or_else(|| score.clone())
}

pub fn revoice_score_transposed(score: &Score, rh: usize, lh: usize, transpose_semitones: i32, new_key: &str) -> Score {
    let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
    let all = parse_all_with_fingers(&text, rh, lh);
    let mut result = all.into_iter()
        .find(|s| s.number == score.number)
        .unwrap_or_else(|| score.clone());

    if transpose_semitones != 0 {
        // Transpose all MIDI values and re-compute string numbers
        let new_root = key_to_pc(new_key);
        for event in &mut result.events {
            if let ScoreEvent::Note { melody_midi, rh_midi, lh_midi, rh_strings, lh_strings, .. } = event {
                *melody_midi += transpose_semitones;
                for m in rh_midi.iter_mut() { *m += transpose_semitones; }
                for m in lh_midi.iter_mut() { *m += transpose_semitones; }
                // Recompute string numbers in the new key
                *rh_strings = rh_midi.iter()
                    .filter_map(|&m| midi_to_harp_string(m, new_root)
                        .map(|s| to_relative_string(s, new_root) - (STRING_NUMBER_OFFSET - 7)))
                    .collect();
                *lh_strings = lh_midi.iter()
                    .filter_map(|&m| midi_to_harp_string(m, new_root)
                        .map(|s| to_relative_string(s, new_root) - (STRING_NUMBER_OFFSET - 7)))
                    .collect();
            }
        }
        result.key = new_key.to_string();
        result.key_root = new_root;
    }

    result
}

fn parse_all_with_fingers(text: &str, max_rh: usize, max_lh: usize) -> Vec<Score> {
    let mut scores = Vec::new();
    let tunes: Vec<&str> = text.split("\nX:").collect();

    for (idx, tune) in tunes.iter().enumerate() {
        if tune.trim().is_empty() { continue; }

        let mut number = String::new();
        let mut title = String::new();
        let mut meter_num: u8 = 4;
        let mut meter_den: u8 = 4;
        let mut key = "C".to_string();
        let mut default_len: f32 = 0.25;
        let mut tempo: u16 = 120;
        let mut voice_lines: [Vec<String>; 4] = [vec![], vec![], vec![], vec![]];
        let mut has_voices = false;

        for line in tune.lines() {
            let trimmed = line.trim();
            if trimmed.parse::<u32>().is_ok() && number.is_empty() {
                number = trimmed.to_string();
            } else if let Some(t) = line.strip_prefix("T:") {
                if title.is_empty() { title = t.trim().to_string(); }
            } else if let Some(m) = line.strip_prefix("M:") {
                let m = m.trim().split('%').next().unwrap_or("4/4").trim();
                let parts: Vec<&str> = m.split('/').collect();
                if parts.len() == 2 {
                    meter_num = parts[0].parse().unwrap_or(4);
                    meter_den = parts[1].parse().unwrap_or(4);
                }
            } else if let Some(l) = line.strip_prefix("L:") {
                let l = l.trim().split('%').next().unwrap_or("1/4").trim();
                let parts: Vec<&str> = l.split('/').collect();
                if parts.len() == 2 {
                    let num: f32 = parts[0].parse().unwrap_or(1.0);
                    let den: f32 = parts[1].parse().unwrap_or(4.0);
                    default_len = num / den;
                }
            } else if let Some(k) = line.strip_prefix("K:") {
                let k = k.trim().split('%').next().unwrap_or("C").trim()
                    .split_whitespace().next().unwrap_or("C");
                if k.starts_with(|c: char| c.is_ascii_uppercase()) {
                    key = k.to_string();
                }
            } else if let Some(rest) = line.strip_prefix("[V: S1V1]") {
                let music = rest.trim();
                if let Some(q) = music.strip_prefix("[Q:") {
                    if let Some(end) = q.find(']') {
                        let qstr = &q[..end];
                        if let Some(eq) = qstr.find('=') {
                            tempo = qstr[eq+1..].trim().parse().unwrap_or(120);
                        }
                        voice_lines[0].push(music[end + 2..].trim().to_string());
                    } else {
                        voice_lines[0].push(music.to_string());
                    }
                } else {
                    voice_lines[0].push(music.to_string());
                }
                has_voices = true;
            } else if let Some(rest) = line.strip_prefix("[V: S1V2]") {
                voice_lines[1].push(rest.trim().to_string()); has_voices = true;
            } else if let Some(rest) = line.strip_prefix("[V: S2V1]") {
                voice_lines[2].push(rest.trim().to_string()); has_voices = true;
            } else if let Some(rest) = line.strip_prefix("[V: S2V2]") {
                voice_lines[3].push(rest.trim().to_string()); has_voices = true;
            }
        }

        if !has_voices || voice_lines[0].is_empty() { continue; }
        if number.is_empty() {
            if idx > 0 {
                if let Some(n) = tune.split_whitespace().next() {
                    if n.parse::<u32>().is_ok() { number = n.to_string(); }
                }
            }
        }

        let key_root = key_to_pc(&key);
        let key_sig_acc = key_sig_accidentals(&key);
        let beats_per_unit = default_len / 0.25;

        // Parse all 4 voices with durations
        let soprano_tokens = parse_voice_tokens(&voice_lines[0].join(" "), &key_sig_acc);
        let alto_tokens = parse_voice_tokens(&voice_lines[1].join(" "), &key_sig_acc);
        let tenor_tokens = parse_voice_tokens(&voice_lines[2].join(" "), &key_sig_acc);
        let bass_tokens = parse_voice_tokens(&voice_lines[3].join(" "), &key_sig_acc);

        // Extract just the MIDI notes (without durations) for inner voices
        let a_notes: Vec<i32> = alto_tokens.iter().filter_map(|t| match t { Token::Note(m, _) => Some(*m), _ => None }).collect();
        let t_notes: Vec<i32> = tenor_tokens.iter().filter_map(|t| match t { Token::Note(m, _) => Some(*m), _ => None }).collect();
        let b_notes: Vec<i32> = bass_tokens.iter().filter_map(|t| match t { Token::Note(m, _) => Some(*m), _ => None }).collect();

        // Walk soprano with durations, apply voicing engine
        let mut events = Vec::new();
        let mut last_voicing: Option<Voicing> = None;
        let mut note_idx = 0usize;

        for tok in &soprano_tokens {
            match tok {
                Token::Bar => { events.push(ScoreEvent::Bar); }
                Token::Rest(dur_mult) => {
                    events.push(ScoreEvent::Rest { beats: dur_mult * beats_per_unit });
                    note_idx += 1;
                }
                Token::Note(s_midi, dur_mult) => {
                    let beats = dur_mult * beats_per_unit;
                    let si = note_idx;
                    note_idx += 1;

                    let a_midi = a_notes.get(si).copied();
                    let t_midi = t_notes.get(si).copied();
                    let b_midi = b_notes.get(si).copied();

                    let voicing = voice_from_satb_fingers(*s_midi, a_midi, t_midi, b_midi, key_root, last_voicing.as_ref(), max_rh, max_lh);
                    let snapped = snap_to_diatonic(*s_midi, key_root);

                    if let Some(ref v) = voicing {
                        // Detect chord change
                        let voicing_key = format!("{:?}|{:?}", v.rh_strings, v.lh_strings);
                        let last_key = last_voicing.as_ref()
                            .map(|lv| format!("{:?}|{:?}", lv.rh_strings, lv.lh_strings))
                            .unwrap_or_default();
                        let is_chord_change = voicing_key != last_key;

                        let rh_strs: Vec<i32> = v.rh_strings.iter()
                            .map(|&s| crate::music::to_relative_string(s, key_root)).collect();
                        let lh_strs: Vec<i32> = v.lh_strings.iter()
                            .map(|&s| crate::music::to_relative_string(s, key_root)).collect();
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
                                let all_notes: Vec<i32> = v.rh_midi.iter().chain(v.lh_midi.iter()).copied().collect();
                                Some(identify_chord_names(&all_notes, key_root))
                            } else { None },
                            rh_chord: if is_chord_change { Some(identify_chord_names(&v.rh_midi, key_root)) } else { None },
                            lh_chord: if is_chord_change { Some(identify_chord_names(&v.lh_midi, key_root)) } else { None },
                            is_chord_change,
                            display_rh: vec![],
                            display_lh: vec![],
                        });
                        last_voicing = voicing;
                    } else {
                        let mel_str = crate::music::midi_to_harp_string(snapped, key_root)
                            .map(|s| crate::music::to_relative_string(s, key_root))
                            .unwrap_or(0);
                        events.push(ScoreEvent::Note {
                            melody_midi: snapped,
                            melody_string: mel_str,
                            rh_midi: vec![snapped],
                            lh_midi: vec![],
                            rh_strings: vec![],
                            lh_strings: vec![],
                            beats,
                            chord_name: None,
                            rh_chord: None,
                            lh_chord: None,
                            is_chord_change: false,
                            display_rh: vec![],
                            display_lh: vec![],
                        });
                    }
                }
            }
        }

        if !events.is_empty() {
            scores.push(Score {
                title, number, key, key_root,
                meter_num, meter_den, tempo,
                events,
            });
        }
    }

    scores
}
