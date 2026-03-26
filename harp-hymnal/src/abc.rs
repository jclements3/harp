use crate::music::*;
use crate::voicing::{voice_from_satb, Voicing};


#[derive(Debug, Clone)]
pub enum Beat {
    Bar,
    Note {
        midi: i32,
        chord: Option<String>,
        full_chord: Option<String>,
        rh_chord: Option<String>,
        lh_chord: Option<String>,
        rh_strings: Vec<i32>,
        lh_strings: Vec<i32>,
        rh_midi: Vec<i32>,
        lh_midi: Vec<i32>,
    },
}

#[derive(Debug, Clone)]
pub struct Hymn {
    pub number: String,
    pub title: String,
    #[allow(dead_code)]
    pub meter: String,
    pub key: String,
    pub key_root: i32,
    pub beats: Vec<Beat>,
}

// ── ABC note parser ──

const ABC_BASE: [(char, i32); 7] = [
    ('C', 48), ('D', 50), ('E', 52), ('F', 53), ('G', 55), ('A', 57), ('B', 59),
];

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

#[derive(Debug)]
enum Token {
    Bar,
    Rest,
    Note(i32),
}

fn parse_voice_tokens(music: &str, key_sig_acc: &[i32; 7]) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut chars = music.chars().peekable();
    let mut buf = String::new();

    while let Some(&ch) = chars.peek() {
        if ch == '|' {
            tokens.push(Token::Bar);
            chars.next();
            // skip repeat markers
            while chars.peek().map_or(false, |&c| matches!(c, '|' | ']' | '[' | ':')) {
                chars.next();
            }
        } else if ch == 'z' || ch == 'Z' {
            tokens.push(Token::Rest);
            chars.next();
            while chars.peek().map_or(false, |c| c.is_ascii_digit() || *c == '/') {
                chars.next();
            }
        } else if matches!(ch, '^' | '_' | '=') || ch.is_ascii_alphabetic() {
            buf.clear();
            // accidentals
            while chars.peek().map_or(false, |&c| matches!(c, '^' | '_' | '=')) {
                buf.push(chars.next().unwrap());
            }
            // letter
            if let Some(&c) = chars.peek() {
                if c.is_ascii_alphabetic() && !"zZxXwWhH".contains(c) {
                    buf.push(chars.next().unwrap());
                    // octave markers
                    while chars.peek().map_or(false, |&c| c == ',' || c == '\'') {
                        buf.push(chars.next().unwrap());
                    }
                    // duration
                    while chars.peek().map_or(false, |c| c.is_ascii_digit() || *c == '/') {
                        chars.next();
                    }
                    if let Some(midi) = abc_note_to_midi(&buf, key_sig_acc) {
                        tokens.push(Token::Note(midi));
                    }
                } else {
                    chars.next();
                }
            }
        } else {
            chars.next();
        }
    }
    tokens
}

pub fn parse_abc(text: &str) -> Vec<Hymn> {
    if text.contains("[V: S1V1]") {
        parse_open_hymnal(text)
    } else {
        parse_lead_sheets(text)
    }
}

fn parse_open_hymnal(text: &str) -> Vec<Hymn> {
    let mut hymns = Vec::new();
    let tunes: Vec<&str> = text.split("\nX:").collect();

    for (idx, tune) in tunes.iter().enumerate() {
        if tune.trim().is_empty() { continue; }

        let mut number = String::new();
        let mut title = String::new();
        let mut meter = "4/4".to_string();
        let mut key = "C".to_string();
        let mut voices: [Vec<String>; 4] = [vec![], vec![], vec![], vec![]];

        for line in tune.lines() {
            let trimmed = line.trim();
            if trimmed.parse::<u32>().is_ok() && number.is_empty() {
                number = trimmed.to_string();
            } else if let Some(t) = line.strip_prefix("T:") {
                if title.is_empty() { title = t.trim().to_string(); }
            } else if let Some(m) = line.strip_prefix("M:") {
                meter = m.trim().split('%').next().unwrap_or("4/4").trim().to_string();
            } else if let Some(k) = line.strip_prefix("K:") {
                let k = k.trim().split('%').next().unwrap_or("C").trim()
                    .split_whitespace().next().unwrap_or("C");
                if k.starts_with(|c: char| c.is_ascii_uppercase()) {
                    key = k.to_string();
                }
            } else if let Some(rest) = line.strip_prefix("[V: S1V1]") {
                voices[0].push(rest.trim().to_string());
            } else if let Some(rest) = line.strip_prefix("[V: S1V2]") {
                voices[1].push(rest.trim().to_string());
            } else if let Some(rest) = line.strip_prefix("[V: S2V1]") {
                voices[2].push(rest.trim().to_string());
            } else if let Some(rest) = line.strip_prefix("[V: S2V2]") {
                voices[3].push(rest.trim().to_string());
            }
        }

        if voices[0].is_empty() { continue; }
        if number.is_empty() {
            if idx > 0 {
                // First token might be the X: number
                if let Some(n) = tune.split_whitespace().next() {
                    if n.parse::<u32>().is_ok() { number = n.to_string(); }
                }
            }
        }

        let key_root = key_to_pc(&key);
        let key_sig_acc = key_sig_accidentals(&key);

        let soprano = parse_voice_tokens(&voices[0].join(" "), &key_sig_acc);
        let alto = parse_voice_tokens(&voices[1].join(" "), &key_sig_acc);
        let tenor = parse_voice_tokens(&voices[2].join(" "), &key_sig_acc);
        let bass_tokens = parse_voice_tokens(&voices[3].join(" "), &key_sig_acc);

        let a_notes: Vec<i32> = alto.iter().filter_map(|t| match t { Token::Note(m) => Some(*m), _ => None }).collect();
        let t_notes: Vec<i32> = tenor.iter().filter_map(|t| match t { Token::Note(m) => Some(*m), _ => None }).collect();
        let b_notes: Vec<i32> = bass_tokens.iter().filter_map(|t| match t { Token::Note(m) => Some(*m), _ => None }).collect();

        let mut beats = Vec::new();
        let mut last_voicing: Option<Voicing> = None;
        let mut note_idx = 0usize;

        for tok in &soprano {
            match tok {
                Token::Bar => { beats.push(Beat::Bar); }
                Token::Rest => { note_idx += 1; }
                Token::Note(s_midi) => {
                    let si = note_idx;
                    note_idx += 1;

                    let a_midi = a_notes.get(si).copied();
                    let t_midi = t_notes.get(si).copied();
                    let b_midi = b_notes.get(si).copied();

                    let voicing = voice_from_satb(*s_midi, a_midi, t_midi, b_midi, key_root, last_voicing.as_ref());

                    let snapped = snap_to_diatonic(*s_midi, key_root);

                    if let Some(ref v) = voicing {
                        let voicing_key = format!("{:?}|{:?}", v.rh_strings, v.lh_strings);
                        let last_key = last_voicing.as_ref()
                            .map(|lv| format!("{:?}|{:?}", lv.rh_strings, lv.lh_strings))
                            .unwrap_or_default();

                        if voicing_key != last_key {
                            beats.push(Beat::Note {
                                midi: snapped,
                                chord: Some(format!("{}/{}", v.rh_chord, v.lh_chord)),
                                full_chord: Some(v.full_chord.clone()),
                                rh_chord: Some(v.rh_chord.clone()),
                                lh_chord: Some(v.lh_chord.clone()),
                                rh_strings: v.rh_strings.clone(),
                                lh_strings: v.lh_strings.clone(),
                                rh_midi: v.rh_midi.clone(),
                                lh_midi: v.lh_midi.clone(),
                            });
                            last_voicing = voicing;
                        } else {
                            beats.push(Beat::Note {
                                midi: snapped,
                                chord: None, full_chord: None,
                                rh_chord: None, lh_chord: None,
                                rh_strings: vec![], lh_strings: vec![],
                                rh_midi: vec![], lh_midi: vec![],
                            });
                        }
                    } else {
                        beats.push(Beat::Note {
                            midi: snapped,
                            chord: None, full_chord: None,
                            rh_chord: None, lh_chord: None,
                            rh_strings: vec![], lh_strings: vec![],
                            rh_midi: vec![], lh_midi: vec![],
                        });
                    }
                }
            }
        }

        if !beats.is_empty() {
            hymns.push(Hymn { number, title, meter, key, key_root, beats });
        }
    }

    hymns
}

fn parse_lead_sheets(text: &str) -> Vec<Hymn> {
    let mut hymns = Vec::new();

    for tune in text.split("\nX:") {
        if tune.trim().is_empty() { continue; }

        let mut number = String::new();
        let mut title = String::new();
        let mut meter = "4/4".to_string();
        let mut key = "C".to_string();
        let mut music_lines = Vec::new();

        for line in tune.lines() {
            if let Some(n) = line.strip_prefix("X:") { number = n.trim().to_string(); }
            else if let Some(t) = line.strip_prefix("T:") { title = t.trim().to_string(); }
            else if let Some(m) = line.strip_prefix("M:") { meter = m.trim().to_string(); }
            else if let Some(k) = line.strip_prefix("K:") { key = k.trim().to_string(); }
            else if !line.is_empty() && !line.starts_with('%') && !line.starts_with("L:") {
                music_lines.push(line.to_string());
            }
        }

        let key_root = key_to_pc(&key);
        let key_sig_acc = key_sig_accidentals(&key);
        let tokens = parse_voice_tokens(&music_lines.join(" "), &key_sig_acc);

        let mut beats = Vec::new();
        for tok in &tokens {
            match tok {
                Token::Bar => beats.push(Beat::Bar),
                Token::Rest => {}
                Token::Note(midi) => {
                    beats.push(Beat::Note {
                        midi: *midi,
                        chord: None, full_chord: None,
                        rh_chord: None, lh_chord: None,
                        rh_strings: vec![], lh_strings: vec![],
                        rh_midi: vec![], lh_midi: vec![],
                    });
                }
            }
        }

        if !beats.is_empty() {
            hymns.push(Hymn { number, title, meter, key, key_root, beats });
        }
    }

    hymns
}
