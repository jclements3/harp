use crate::chord_symbol;
use crate::music::{MAJOR_SCALE, PC_NAMES};

/// Convert ASCII chord symbol to traditional notation.
/// e.g. "vii%7" → "viiø7", "vii-" → "vii°"
pub fn chord_to_traditional(ascii: &str) -> String {
    ascii
        .replace("%7", "\u{00f8}7") // ø7
        .replace("-7", "\u{00b0}7") // °7
        .replace('-', "\u{00b0}")   // °
        .replace("M7", "maj7")
        .replace("~2", "sus2")
        .replace("~4", "sus4")
        // Inversion: ^1 → (1st inv.)
        .replace("^1", " (1st inv.)")
        .replace("^2", " (2nd inv.)")
        .replace("^3", " (3rd inv.)")
}

const DEGREE_NAMES: [&str; 7] = ["one", "two", "three", "four", "five", "six", "seven"];
const INV_WORDS: [&str; 4] = ["", "first", "second", "third"];

/// Convert ASCII chord symbol to spoken form.
/// e.g. "V7^2" → "five dominant seven second inversion"
pub fn chord_to_say(ascii: &str) -> String {
    let mut s = ascii.to_string();
    let mut inv = String::new();

    // Extract inversion
    if let Some(pos) = s.find('^') {
        if pos + 1 < s.len() {
            let digit = s.as_bytes()[pos + 1];
            if digit.is_ascii_digit() {
                let idx = (digit - b'0') as usize;
                if idx < INV_WORDS.len() {
                    inv = format!(" {} inversion", INV_WORDS[idx]);
                }
                s = format!("{}{}", &s[..pos], &s[pos + 2..]);
            }
        }
    }

    let upper = s.to_uppercase();
    let mut degree = String::new();
    let mut rest = s.clone();

    let romans: &[(&str, usize)] = &[
        ("VII", 7),
        ("VI", 6),
        ("IV", 4),
        ("V", 5),
        ("III", 3),
        ("II", 2),
        ("I", 1),
    ];

    for &(r, d) in romans {
        if upper.starts_with(r)
            || upper.starts_with(&format!("B{}", r))
            || upper.starts_with(&format!("#{}", r))
        {
            let mut prefix = String::new();
            if s.starts_with('b') {
                prefix = "flat ".to_string();
                s = s[1..].to_string();
            } else if s.starts_with('#') {
                prefix = "sharp ".to_string();
                s = s[1..].to_string();
            }
            degree = format!("{}{}", prefix, DEGREE_NAMES[d - 1]);
            rest = s[r.len()..].to_string();
            break;
        }
    }

    let is_minor = s
        .chars()
        .next()
        .map(|c| c.is_lowercase() && c.is_alphabetic())
        .unwrap_or(false);

    let quality = if rest.contains("-7") {
        " diminished seven"
    } else if rest.contains("%7") {
        " half-diminished"
    } else if rest.contains("M7") {
        " major seven"
    } else if rest.contains('7') {
        if is_minor {
            " minor seven"
        } else {
            " dominant seven"
        }
    } else if rest.contains('-') {
        " diminished"
    } else if rest.contains('+') {
        " augmented"
    } else if rest.contains("~4") {
        " sus four"
    } else if rest.contains("~2") {
        " sus two"
    } else if is_minor {
        " minor"
    } else {
        " major"
    };

    format!("{}{}{}", degree, quality, inv)
}

/// Show the actual pitch classes for a chord in a given key.
/// e.g. "V" in C → "G B D"
pub fn chord_to_example(ascii: &str, key_root: u8) -> String {
    let parsed = match chord_symbol::parse_chord_symbol(ascii) {
        Ok(c) => c,
        Err(_) => return String::new(),
    };
    let root_pc = ((key_root as i16
        + MAJOR_SCALE[parsed.root_degree as usize % 7] as i16
        + parsed.root_accidental as i16
        + 12)
        % 12) as u8;

    parsed
        .intervals
        .iter()
        .map(|&iv| PC_NAMES[((root_pc + iv) % 12) as usize])
        .collect::<Vec<_>>()
        .join(" ")
}

/// Show the interval structure of a chord.
/// e.g. "V7" → "1 3 5 b7"
pub fn chord_to_intervals(ascii: &str) -> String {
    let parsed = match chord_symbol::parse_chord_symbol(ascii) {
        Ok(c) => c,
        Err(_) => return String::new(),
    };

    let iv_names: &[(u8, &str)] = &[
        (0, "1"),
        (1, "b2"),
        (2, "2"),
        (3, "b3"),
        (4, "3"),
        (5, "4"),
        (6, "b5"),
        (7, "5"),
        (8, "#5"),
        (9, "6"),
        (10, "b7"),
        (11, "7"),
    ];

    parsed
        .intervals
        .iter()
        .map(|&iv| {
            iv_names
                .iter()
                .find(|&&(s, _)| s == iv % 12)
                .map(|&(_, name)| name)
                .unwrap_or("?")
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chord_to_traditional() {
        assert_eq!(chord_to_traditional("vii%7"), "vii\u{00f8}7");
        assert_eq!(chord_to_traditional("vii-"), "vii\u{00b0}");
        assert_eq!(chord_to_traditional("IM7"), "Imaj7");
        assert_eq!(chord_to_traditional("V~4"), "Vsus4");
        assert_eq!(chord_to_traditional("V^1"), "V (1st inv.)");
    }

    #[test]
    fn test_chord_to_say() {
        assert_eq!(chord_to_say("V7"), "five dominant seven");
        assert_eq!(chord_to_say("ii"), "two minor");
        assert_eq!(chord_to_say("I"), "one major");
        assert_eq!(chord_to_say("bVII"), "flat seven major");
        assert_eq!(chord_to_say("V7^2"), "five dominant seven second inversion");
        assert_eq!(chord_to_say("vii-"), "seven diminished");
    }

    #[test]
    fn test_chord_to_example() {
        // V in C = G B D
        assert_eq!(chord_to_example("V", 0), "G B D");
        // I in G = G B D
        assert_eq!(chord_to_example("I", 7), "G B D");
        // ii in C = D F A
        assert_eq!(chord_to_example("ii", 0), "D F A");
    }

    #[test]
    fn test_chord_to_intervals() {
        assert_eq!(chord_to_intervals("I"), "1 3 5");
        assert_eq!(chord_to_intervals("ii"), "1 b3 5");
        assert_eq!(chord_to_intervals("V7"), "1 3 5 b7");
        assert_eq!(chord_to_intervals("IM7"), "1 3 5 7");
    }
}
