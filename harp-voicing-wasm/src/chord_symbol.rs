/// Parsed representation of an ASCII Roman numeral chord symbol.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedChord {
    /// Root as scale degree 0-6 (I=0, ii=1, ... vii=6)
    pub root_degree: u8,
    /// Chromatic adjustment to root: -1 for b, +1 for #, 0 for natural
    pub root_accidental: i8,
    /// Chord tones as semitone intervals above root
    pub intervals: Vec<u8>,
    /// Inversion: 0=root, 1=1st, 2=2nd, 3=3rd
    pub inversion: u8,
    /// Original symbol string
    pub symbol: String,
}

impl ParsedChord {
    /// Number of distinct chord tones.
    pub fn tone_count(&self) -> usize {
        self.intervals.len()
    }
}

/// Parse a chord symbol like "IM7^2", "vii%7", "bVII", "#iv-7", "I~42".
pub fn parse_chord_symbol(s: &str) -> Result<ParsedChord, String> {
    let original = s.to_string();
    let bytes = s.as_bytes();
    let len = bytes.len();
    let mut pos = 0;

    // 1. Optional accidental on root: b or #
    let root_accidental = if pos < len && bytes[pos] == b'b' {
        pos += 1;
        -1i8
    } else if pos < len && bytes[pos] == b'#' {
        pos += 1;
        1i8
    } else {
        0i8
    };

    // 2. Roman numeral — determine degree and whether base is major/minor
    let (degree, is_minor, consumed) = parse_roman(&s[pos..])?;
    pos += consumed;

    // Start with major or minor triad
    let mut intervals: Vec<u8> = if is_minor {
        vec![0, 3, 7] // minor triad
    } else {
        vec![0, 4, 7] // major triad
    };

    // 3. Quality modifiers
    let rest = &s[pos..];

    // Check for augmented +
    if rest.starts_with('+') {
        intervals = if is_minor {
            vec![0, 3, 8] // minor augmented (rare)
        } else {
            vec![0, 4, 8] // augmented
        };
        pos += 1;
    }
    // Check for diminished -
    else if rest.starts_with('-') && !rest.starts_with("-7") {
        intervals = vec![0, 3, 6]; // diminished triad
        pos += 1;
    }

    let rest = &s[pos..];

    // 4. Seventh extensions
    if rest.starts_with("M7") {
        // Major 7th
        intervals.push(11);
        pos += 2;
    } else if rest.starts_with("%7") {
        // Half-diminished 7th: dim triad + minor 7th
        intervals = vec![0, 3, 6, 10];
        pos += 2;
    } else if rest.starts_with("-7") {
        // Diminished 7th: dim triad + dim 7th
        intervals = vec![0, 3, 6, 9];
        pos += 2;
    } else if rest.starts_with('7') {
        // Dominant 7th (on major triad) or minor 7th (on minor triad)
        intervals.push(10);
        pos += 1;
    }

    let rest = &s[pos..];

    // 5. Suspended: ~2 or ~4, optionally followed by add digit
    if rest.starts_with("~2") {
        // Replace 3rd with 2nd
        intervals.retain(|&x| x != 3 && x != 4);
        intervals.push(2);
        intervals.sort();
        intervals.dedup();
        pos += 2;

        // Check for add digit after sus type
        let rest2 = &s[pos..];
        if let Some((add_semitones, consumed)) = parse_add_digit(rest2) {
            if !intervals.contains(&add_semitones) {
                intervals.push(add_semitones);
                intervals.sort();
            }
            pos += consumed;
        }
    } else if rest.starts_with("~4") {
        // Replace 3rd with 4th (5 semitones)
        intervals.retain(|&x| x != 3 && x != 4);
        intervals.push(5);
        intervals.sort();
        intervals.dedup();
        pos += 2;

        // Check for add digit after sus type
        let rest2 = &s[pos..];
        if let Some((add_semitones, consumed)) = parse_add_digit(rest2) {
            if !intervals.contains(&add_semitones) {
                intervals.push(add_semitones);
                intervals.sort();
            }
            pos += consumed;
        }
    } else {
        // 6. Add digit (not after sus, not 7 which was already handled)
        let rest = &s[pos..];
        if let Some((add_semitones, consumed)) = parse_add_digit(rest) {
            if !intervals.contains(&add_semitones) {
                intervals.push(add_semitones);
                intervals.sort();
            }
            pos += consumed;
        }
    }

    // 7. Inversion: ^N
    let mut inversion = 0u8;
    let rest = &s[pos..];
    if rest.starts_with('^') {
        pos += 1;
        if pos < len && bytes[pos].is_ascii_digit() {
            inversion = (bytes[pos] - b'0') as u8;
            pos += 1;
        }
    }

    // Validate inversion doesn't exceed chord tones
    if inversion as usize >= intervals.len() {
        return Err(format!(
            "Inversion ^{} exceeds chord tone count {} in '{}'",
            inversion,
            intervals.len(),
            original
        ));
    }

    Ok(ParsedChord {
        root_degree: degree,
        root_accidental,
        intervals,
        inversion,
        symbol: original,
    })
}

/// Parse a Roman numeral at the start of `s`.
/// Returns (degree 0-6, is_minor, chars_consumed).
fn parse_roman(s: &str) -> Result<(u8, bool, usize), String> {
    let upper = s.to_uppercase();

    // Try longest matches first
    let candidates = [
        ("VII", 6),
        ("VI", 5),
        ("IV", 3),
        ("V", 4),
        ("III", 2),
        ("II", 1),
        ("I", 0),
    ];

    for &(numeral, degree) in &candidates {
        if upper.starts_with(numeral) {
            let consumed = numeral.len();
            // Check case of original to determine major/minor
            let original_chunk = &s[..consumed];
            let is_minor = original_chunk.chars().all(|c| c.is_lowercase());
            return Ok((degree, is_minor, consumed));
        }
    }

    Err(format!("No Roman numeral found at start of '{}'", s))
}

/// Parse an add-note digit: single char 2,4,6,8,9,A (hex for 10).
/// Returns (semitone interval, chars consumed) or None.
fn parse_add_digit(s: &str) -> Option<(u8, usize)> {
    if s.is_empty() {
        return None;
    }
    let ch = s.as_bytes()[0];
    // Don't parse if next char is part of ^inversion
    if ch == b'^' {
        return None;
    }
    match ch {
        b'2' => Some((2, 1)),   // major 2nd
        b'4' => Some((5, 1)),   // perfect 4th
        b'6' => Some((9, 1)),   // major 6th
        b'8' => Some((12, 1)),  // octave (root doubling)
        b'9' => Some((14, 1)),  // major 9th (= 2nd + octave)
        b'A' => Some((16, 1)),  // major 10th (= 3rd + octave)
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_major_triad() {
        let c = parse_chord_symbol("I").unwrap();
        assert_eq!(c.root_degree, 0);
        assert_eq!(c.intervals, vec![0, 4, 7]);
        assert_eq!(c.inversion, 0);
    }

    #[test]
    fn test_minor_triad() {
        let c = parse_chord_symbol("ii").unwrap();
        assert_eq!(c.root_degree, 1);
        assert_eq!(c.intervals, vec![0, 3, 7]);
    }

    #[test]
    fn test_diminished() {
        let c = parse_chord_symbol("vii-").unwrap();
        assert_eq!(c.root_degree, 6);
        assert_eq!(c.intervals, vec![0, 3, 6]);
    }

    #[test]
    fn test_augmented() {
        let c = parse_chord_symbol("I+").unwrap();
        assert_eq!(c.intervals, vec![0, 4, 8]);
    }

    #[test]
    fn test_dominant_7th() {
        let c = parse_chord_symbol("V7").unwrap();
        assert_eq!(c.root_degree, 4);
        assert_eq!(c.intervals, vec![0, 4, 7, 10]);
    }

    #[test]
    fn test_major_7th() {
        let c = parse_chord_symbol("IM7").unwrap();
        assert_eq!(c.intervals, vec![0, 4, 7, 11]);
    }

    #[test]
    fn test_minor_7th() {
        let c = parse_chord_symbol("ii7").unwrap();
        assert_eq!(c.intervals, vec![0, 3, 7, 10]);
    }

    #[test]
    fn test_half_dim_7th() {
        let c = parse_chord_symbol("vii%7").unwrap();
        assert_eq!(c.intervals, vec![0, 3, 6, 10]);
    }

    #[test]
    fn test_dim_7th() {
        let c = parse_chord_symbol("vii-7").unwrap();
        assert_eq!(c.intervals, vec![0, 3, 6, 9]);
    }

    #[test]
    fn test_sus2() {
        let c = parse_chord_symbol("I~2").unwrap();
        assert_eq!(c.intervals, vec![0, 2, 7]);
    }

    #[test]
    fn test_sus4() {
        let c = parse_chord_symbol("I~4").unwrap();
        assert_eq!(c.intervals, vec![0, 5, 7]);
    }

    #[test]
    fn test_sus4_add2() {
        let c = parse_chord_symbol("I~42").unwrap();
        assert_eq!(c.intervals, vec![0, 2, 5, 7]);
    }

    #[test]
    fn test_inversion() {
        let c = parse_chord_symbol("IM7^2").unwrap();
        assert_eq!(c.intervals, vec![0, 4, 7, 11]);
        assert_eq!(c.inversion, 2);
    }

    #[test]
    fn test_flat_root() {
        let c = parse_chord_symbol("bVII").unwrap();
        assert_eq!(c.root_degree, 6);
        assert_eq!(c.root_accidental, -1);
        assert_eq!(c.intervals, vec![0, 4, 7]);
    }

    #[test]
    fn test_sharp_root() {
        let c = parse_chord_symbol("#iv-").unwrap();
        assert_eq!(c.root_degree, 3);
        assert_eq!(c.root_accidental, 1);
        assert_eq!(c.intervals, vec![0, 3, 6]);
    }

    #[test]
    fn test_complex_symbol() {
        let c = parse_chord_symbol("#iv%7^2").unwrap();
        assert_eq!(c.root_degree, 3);
        assert_eq!(c.root_accidental, 1);
        assert_eq!(c.intervals, vec![0, 3, 6, 10]);
        assert_eq!(c.inversion, 2);
    }

    #[test]
    fn test_vi7_with_inversion() {
        let c = parse_chord_symbol("vi7^3").unwrap();
        assert_eq!(c.root_degree, 5);
        assert_eq!(c.intervals, vec![0, 3, 7, 10]);
        assert_eq!(c.inversion, 3);
    }

    #[test]
    fn test_bii_half_dim() {
        let c = parse_chord_symbol("bii%7").unwrap();
        assert_eq!(c.root_degree, 1);
        assert_eq!(c.root_accidental, -1);
        assert_eq!(c.intervals, vec![0, 3, 6, 10]);
    }

    #[test]
    fn test_sus4_inversion() {
        let c = parse_chord_symbol("V~4^2").unwrap();
        assert_eq!(c.root_degree, 4);
        assert_eq!(c.intervals, vec![0, 5, 7]);
        assert_eq!(c.inversion, 2);
    }

    #[test]
    fn test_all_143_leadsheet_symbols() {
        let symbols = [
            "bii-", "bII+", "bii%7", "bii7", "bii%7^2", "biii-", "bIII", "bIII+",
            "bvi-", "bVI+", "bvi-7", "bVII", "bVII^1", "bvii^2", "bVII^2", "bVII~2",
            "bVIIM7", "i", "i-", "I", "I+", "i^1", "I^1", "I^2", "I~2", "I~4",
            "i-7", "I7", "I7^2", "I7^3", "ii", "ii-", "II", "ii^1", "II^1", "ii^2",
            "II^2", "II~2", "II~4", "II~4^2", "ii-7", "ii7", "II7", "ii7^2", "ii7^3",
            "II7^3", "iii", "iii-", "III", "III+", "iii^1", "III^1", "iii^2", "III^2",
            "III~4", "III~4^2", "iii-7", "iii7", "III7", "iii7^2", "III7^2", "iii7^3",
            "III7^3", "IIIM7", "IM7", "IM7^2", "IM7^3", "#iv", "#iv-", "iv", "iv-",
            "#IV+", "IV", "IV+", "#iv^1", "iv^1", "IV^1", "IV^2", "IV~2", "#IV~4",
            "IV~4", "IV~4^2", "#iv%7", "#iv-7", "#iv7", "iv%7", "iv-7", "IV7",
            "#iv%7^2", "IV7^2", "#iv%7^3", "IV7^3", "IVM7", "IVM7^2", "IVM7^3",
            "v", "v-", "V", "V+", "v^1", "V^1", "V^2", "V~2", "V~4", "V~4^2",
            "v-7", "v7", "V7", "V7^2", "V7^3", "vi", "vi-", "VI", "VI+", "vi^1",
            "VI^1", "vi^2", "VI^2", "VI~2", "VI~4", "VI~4^2", "vi-7", "vi7", "VI7",
            "vi7^2", "VI7^2", "vi7^3", "VI7^3", "vii", "vii-", "VII", "VII+",
            "vii^1", "VII^1", "vii^2", "VII~4", "vii%7", "vii-7", "vii7", "VII7",
            "vii%7^2", "vii%7^3", "VM7^2",
        ];
        let mut failures = Vec::new();
        for s in &symbols {
            if let Err(e) = parse_chord_symbol(s) {
                failures.push(format!("{}: {}", s, e));
            }
        }
        assert!(
            failures.is_empty(),
            "Failed to parse {} symbols:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
