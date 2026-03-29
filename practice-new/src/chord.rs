use crate::music::*;

// ── Chord templates for identification ──

#[allow(dead_code)]
struct ChordTemplate {
    pcs: &'static [i32],
    name: &'static str,
    lower: bool,
    suffix: &'static str,
}

const CHORD_TEMPLATES: &[ChordTemplate] = &[
    ChordTemplate { pcs: &[0,4,7],    name: "",     lower: false, suffix: "" },
    ChordTemplate { pcs: &[0,3,7],    name: "m",    lower: true,  suffix: "" },
    ChordTemplate { pcs: &[0,4,7,10], name: "7",    lower: false, suffix: "7" },
    ChordTemplate { pcs: &[0,3,7,10], name: "m7",   lower: true,  suffix: "7" },
    ChordTemplate { pcs: &[0,4,7,11], name: "maj7", lower: false, suffix: "M7" },
    ChordTemplate { pcs: &[0,3,6],    name: "dim",  lower: true,  suffix: "-" },
    ChordTemplate { pcs: &[0,4,8],    name: "aug",  lower: false, suffix: "+" },
    ChordTemplate { pcs: &[0,5,7],    name: "sus4", lower: false, suffix: "~4" },
    ChordTemplate { pcs: &[0,2,7],    name: "sus2", lower: false, suffix: "~2" },
    ChordTemplate { pcs: &[0,3,6,9],  name: "dim7", lower: true,  suffix: "-7" },
    ChordTemplate { pcs: &[0,3,6,10], name: "m7b5", lower: true,  suffix: "%7" },
];

const ROMAN_STRS: [&str; 7] = ["I","II","III","IV","V","VI","VII"];

pub fn identify_chord(midi_notes: &[i32], key_root: i32) -> String {
    if midi_notes.len() < 2 { return String::new(); }

    let pcs: Vec<i32> = {
        let mut s: Vec<i32> = midi_notes.iter().map(|&m| pitch_class(m)).collect();
        s.sort(); s.dedup(); s
    };

    if pcs.len() < 2 { return String::new(); }

    // 2 pitch classes: just show interval
    if pcs.len() == 2 {
        let bass_pc = pitch_class(*midi_notes.iter().min().unwrap());
        for d in 0..7 {
            if (key_root + MAJOR_SCALE[d]) % 12 == bass_pc {
                return ROMAN_STRS[d].to_string();
            }
        }
        return String::new();
    }

    let bass_pc = pitch_class(*midi_notes.iter().min().unwrap());
    let mut best_root = 0i32;
    let mut best_idx = 0usize;
    let mut best_score = i32::MIN;

    for root in 0..12 {
        let intervals: Vec<i32> = pcs.iter().map(|&p| ((p - root) % 12 + 12) % 12).collect();
        for (idx, tmpl) in CHORD_TEMPLATES.iter().enumerate() {
            let matched = intervals.iter().filter(|iv| tmpl.pcs.contains(iv)).count() as i32;
            let extra = intervals.iter().filter(|iv| !tmpl.pcs.contains(iv)).count() as i32;
            let missing = tmpl.pcs.iter().filter(|t| !intervals.contains(t)).count() as i32;
            let mut score = matched * 3 - extra * 2 - missing;
            if bass_pc == root { score += 4; }
            if tmpl.pcs.len() == 3 { score += 2; }
            score -= (idx / 10) as i32;
            if score > best_score {
                best_score = score;
                best_root = root;
                best_idx = idx;
            }
        }
    }

    let deg = (0..7).find(|&d| (key_root + MAJOR_SCALE[d]) % 12 == best_root);
    let deg = match deg { Some(d) => d, None => return String::new() };

    let tmpl = &CHORD_TEMPLATES[best_idx];
    let roman = if tmpl.lower {
        ROMAN_STRS[deg].to_lowercase()
    } else {
        ROMAN_STRS[deg].to_string()
    };
    let mut result = format!("{}{}", roman, tmpl.suffix);

    // Inversion
    let bass_iv = ((bass_pc - best_root) % 12 + 12) % 12;
    if let Some(inv) = tmpl.pcs.iter().position(|&p| p == bass_iv) {
        if inv > 0 {
            result.push_str(&format!("^{}", inv));
        }
    }

    result
}

/// Identify chord using note names (e.g., "G", "Dm7", "A/C#") instead of Roman numerals.
pub fn identify_chord_names(midi_notes: &[i32], key_root: i32) -> String {
    if midi_notes.len() < 2 { return String::new(); }

    let pcs: Vec<i32> = {
        let mut s: Vec<i32> = midi_notes.iter().map(|&m| pitch_class(m)).collect();
        s.sort(); s.dedup(); s
    };

    if pcs.len() < 2 { return String::new(); }

    let bass_pc = pitch_class(*midi_notes.iter().min().unwrap());
    let mut best_root = 0i32;
    let mut best_idx = 0usize;
    let mut best_score = i32::MIN;

    for root in 0..12 {
        let intervals: Vec<i32> = pcs.iter().map(|&p| ((p - root) % 12 + 12) % 12).collect();
        for (idx, tmpl) in CHORD_TEMPLATES.iter().enumerate() {
            let matched = intervals.iter().filter(|iv| tmpl.pcs.contains(iv)).count() as i32;
            let extra = intervals.iter().filter(|iv| !tmpl.pcs.contains(iv)).count() as i32;
            let missing = tmpl.pcs.iter().filter(|t| !intervals.contains(t)).count() as i32;
            let mut score = matched * 3 - extra * 2 - missing;
            if bass_pc == root { score += 4; }
            if tmpl.pcs.len() == 3 { score += 2; }
            score -= (idx / 10) as i32;
            if score > best_score {
                best_score = score;
                best_root = root;
                best_idx = idx;
            }
        }
    }

    let root_name = PC_NAMES[best_root as usize % 12];
    let tmpl = &CHORD_TEMPLATES[best_idx];
    let quality = tmpl.name;

    // Just root + quality, no slash bass — the bass note is visible on the staff
    format!("{}{}", root_name, quality)
}

// ── Chord symbol parser (Roman numeral notation) ──

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ParsedChord {
    pub root_degree: usize,
    pub root_accidental: i32,
    pub intervals: Vec<i32>,
    pub inversion: usize,
    pub symbol: String,
}

#[allow(dead_code)]
pub fn parse_chord_symbol(s: &str) -> Option<ParsedChord> {
    let bytes = s.as_bytes();
    let mut pos = 0;

    let mut root_accidental = 0i32;
    if pos < bytes.len() && bytes[pos] == b'b' { root_accidental = -1; pos += 1; }
    else if pos < bytes.len() && bytes[pos] == b'#' { root_accidental = 1; pos += 1; }

    let (degree, is_minor, consumed) = parse_roman(&s[pos..])?;
    pos += consumed;

    let mut intervals = if is_minor { vec![0, 3, 7] } else { vec![0, 4, 7] };
    let rest = &s[pos..];

    if rest.starts_with('+') {
        intervals = if is_minor { vec![0, 3, 8] } else { vec![0, 4, 8] };
        pos += 1;
    } else if rest.starts_with('-') && !rest.starts_with("-7") {
        intervals = vec![0, 3, 6];
        pos += 1;
    }

    let rest = &s[pos..];
    if rest.starts_with("M7") { intervals.push(11); pos += 2; }
    else if rest.starts_with("%7") { intervals = vec![0, 3, 6, 10]; pos += 2; }
    else if rest.starts_with("-7") { intervals = vec![0, 3, 6, 9]; pos += 2; }
    else if rest.starts_with('7') { intervals.push(10); pos += 1; }

    let rest = &s[pos..];
    if rest.starts_with("~2") {
        intervals.retain(|&x| x != 3 && x != 4);
        intervals.push(2); intervals.sort(); intervals.dedup();
        pos += 2;
        if let Some((semi, c)) = parse_add_digit(&s[pos..]) {
            if !intervals.contains(&semi) { intervals.push(semi); intervals.sort(); }
            pos += c;
        }
    } else if rest.starts_with("~4") {
        intervals.retain(|&x| x != 3 && x != 4);
        intervals.push(5); intervals.sort(); intervals.dedup();
        pos += 2;
        if let Some((semi, c)) = parse_add_digit(&s[pos..]) {
            if !intervals.contains(&semi) { intervals.push(semi); intervals.sort(); }
            pos += c;
        }
    } else if let Some((semi, c)) = parse_add_digit(rest) {
        if !intervals.contains(&semi) { intervals.push(semi); intervals.sort(); }
        pos += c;
    }

    let mut inversion = 0usize;
    let rest = &s[pos..];
    if rest.starts_with('^') && rest.len() > 1 {
        if let Some(d) = rest.as_bytes().get(1).and_then(|&b| (b as char).to_digit(10)) {
            inversion = d as usize;
        }
    }

    if inversion >= intervals.len() { return None; }

    Some(ParsedChord {
        root_degree: degree,
        root_accidental,
        intervals,
        inversion,
        symbol: s.to_string(),
    })
}

fn parse_roman(s: &str) -> Option<(usize, bool, usize)> {
    let upper = s.to_uppercase();
    let candidates: &[(&str, usize)] = &[
        ("VII",6),("VI",5),("IV",3),("V",4),("III",2),("II",1),("I",0),
    ];
    for &(numeral, degree) in candidates {
        if upper.starts_with(numeral) {
            let chunk = &s[..numeral.len()];
            let is_minor = chunk.chars().all(|c| c.is_lowercase());
            return Some((degree, is_minor, numeral.len()));
        }
    }
    None
}

fn parse_add_digit(s: &str) -> Option<(i32, usize)> {
    if s.is_empty() || s.starts_with('^') { return None; }
    match s.as_bytes()[0] {
        b'2' => Some((2, 1)),
        b'4' => Some((5, 1)),
        b'6' => Some((9, 1)),
        b'8' => Some((12, 1)),
        b'9' => Some((14, 1)),
        b'A' => Some((16, 1)),
        _ => None,
    }
}
