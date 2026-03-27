/// ABC parser that produces harp-voiced notation events with durations.
///
/// Uses the same voicing engine as harp-hymnal to produce RH and LH MIDI notes,
/// then wraps them with duration info for the notation renderer.

use crate::music::*;
use crate::voicing::{voice_from_satb, Voicing};

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
        rh_midi: Vec<i32>,
        lh_midi: Vec<i32>,
        rh_strings: Vec<i32>,
        lh_strings: Vec<i32>,
        beats: f32,
        chord_name: Option<String>,
        rh_chord: Option<String>,
        lh_chord: Option<String>,
        is_chord_change: bool,
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

pub fn load_embedded_scores() -> Vec<Score> {
    let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
    parse_all(&text)
}

pub fn parse_all(text: &str) -> Vec<Score> {
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

                    let voicing = voice_from_satb(*s_midi, a_midi, t_midi, b_midi, key_root, last_voicing.as_ref());
                    let snapped = snap_to_diatonic(*s_midi, key_root);

                    if let Some(ref v) = voicing {
                        // Detect chord change
                        let voicing_key = format!("{:?}|{:?}", v.rh_strings, v.lh_strings);
                        let last_key = last_voicing.as_ref()
                            .map(|lv| format!("{:?}|{:?}", lv.rh_strings, lv.lh_strings))
                            .unwrap_or_default();
                        let is_chord_change = voicing_key != last_key;

                        events.push(ScoreEvent::Note {
                            melody_midi: snapped,
                            rh_midi: v.rh_midi.clone(),
                            lh_midi: v.lh_midi.clone(),
                            rh_strings: v.rh_strings.iter().map(|&s| s - STRING_NUMBER_OFFSET).collect(),
                            lh_strings: v.lh_strings.iter().map(|&s| s - STRING_NUMBER_OFFSET).collect(),
                            beats,
                            chord_name: if is_chord_change { Some(v.full_chord.clone()) } else { None },
                            rh_chord: if is_chord_change { Some(v.rh_chord.clone()) } else { None },
                            lh_chord: if is_chord_change { Some(v.lh_chord.clone()) } else { None },
                            is_chord_change,
                        });
                        last_voicing = voicing;
                    } else {
                        events.push(ScoreEvent::Note {
                            melody_midi: snapped,
                            rh_midi: vec![snapped],
                            lh_midi: vec![],
                            rh_strings: vec![],
                            lh_strings: vec![],
                            beats,
                            chord_name: None,
                            rh_chord: None,
                            lh_chord: None,
                            is_chord_change: false,
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
