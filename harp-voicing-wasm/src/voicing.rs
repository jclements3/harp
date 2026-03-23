use std::collections::HashSet;

use crate::chord_symbol::ParsedChord;
use crate::harp::*;
use crate::music::*;

/// Generate valid voicings for a chord given melody note and key.
/// Returns voicings sorted by score descending, capped at `max_results`.
pub fn enumerate_voicings(
    melody_midi: u8,
    chord: &ParsedChord,
    key_root: u8,
    max_results: usize,
) -> Vec<Voicing> {
    // Resolve chord root pitch class
    let root_pc = ((key_root as i16
        + MAJOR_SCALE[chord.root_degree as usize % 7] as i16
        + chord.root_accidental as i16
        + 12) % 12) as u8;

    // Resolve all chord tone pitch classes
    let chord_pcs: Vec<u8> = chord
        .intervals
        .iter()
        .map(|&iv| (root_pc + iv) % 12)
        .collect();

    // Determine required bass pitch class (for inversions)
    let bass_pc = if chord.inversion > 0 && (chord.inversion as usize) < chord.intervals.len() {
        Some((root_pc + chord.intervals[chord.inversion as usize]) % 12)
    } else {
        Some(root_pc) // root in bass by default
    };

    // Find melody's harp string
    let melody_string = match midi_to_harp_string(melody_midi, key_root) {
        Some(s) => s,
        None => return vec![], // melody not on a diatonic string
    };

    // Generate candidate MIDI notes for each chord tone (all octaves within range)
    let mut candidates: Vec<u8> = Vec::new();
    for &pc in &chord_pcs {
        candidates.extend(all_octaves(pc, HARP_LOW_MIDI, HARP_HIGH_MIDI));
    }
    // Also include octave doublings of root
    candidates.extend(all_octaves(root_pc, HARP_LOW_MIDI, HARP_HIGH_MIDI));
    candidates.sort();
    candidates.dedup();

    // Filter to only diatonic strings in this key
    let candidates: Vec<u8> = candidates
        .into_iter()
        .filter(|&m| midi_to_harp_string(m, key_root).is_some())
        .collect();

    // Split candidates into RH-eligible (below melody, within span) and all others
    let rh_candidates: Vec<u8> = candidates
        .iter()
        .copied()
        .filter(|&m| {
            m < melody_midi && {
                let s = midi_to_harp_string(m, key_root).unwrap();
                melody_string - s < MAX_HAND_SPAN
            }
        })
        .collect();

    // Generate RH configurations: thumb=melody + 0-3 additional notes
    let rh_configs = generate_hand_configs(&rh_candidates, melody_midi, melody_string, key_root);

    let mut voicings = Vec::new();

    for rh in &rh_configs {
        let rh_lowest_string = if rh.len() > 1 {
            midi_to_harp_string(rh[rh.len() - 1], key_root).unwrap()
        } else {
            melody_string
        };

        // LH must be below RH lowest (at least 1 string gap)
        let lh_upper_bound_string = if rh_lowest_string > 0 {
            rh_lowest_string - 1
        } else {
            continue; // no room for LH
        };
        let lh_upper_bound_midi = harp_string_to_midi(lh_upper_bound_string, key_root);

        let lh_candidates: Vec<u8> = candidates
            .iter()
            .copied()
            .filter(|&m| m <= lh_upper_bound_midi && m >= HARP_LOW_MIDI)
            .collect();

        if lh_candidates.is_empty() {
            // RH-only voicing
            let v = build_voicing(rh, &[], key_root, root_pc, &chord_pcs);
            voicings.push(v);
            continue;
        }

        let lh_configs = generate_lh_configs(&lh_candidates, key_root, bass_pc);

        if lh_configs.is_empty() {
            let v = build_voicing(rh, &[], key_root, root_pc, &chord_pcs);
            voicings.push(v);
        }

        for lh in &lh_configs {
            let v = build_voicing(rh, lh, key_root, root_pc, &chord_pcs);
            voicings.push(v);
        }
    }

    // Score and sort
    for v in &mut voicings {
        v.score = score_voicing(v, root_pc, &chord_pcs, bass_pc);
    }
    voicings.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    voicings.dedup_by(|a, b| a.all_midis() == b.all_midis());
    voicings.truncate(max_results);
    voicings
}

/// Generate RH hand configurations: thumb is fixed at melody_midi,
/// pick 0-3 additional notes from candidates (within span, sorted high to low).
fn generate_hand_configs(
    candidates: &[u8],
    thumb_midi: u8,
    thumb_string: u8,
    key_root: u8,
) -> Vec<Vec<u8>> {
    let mut configs = Vec::new();

    // Start with thumb only
    configs.push(vec![thumb_midi]);

    // Pick 1-3 additional notes (sorted descending, limited to closest 6)
    let mut eligible: Vec<u8> = candidates
        .iter()
        .copied()
        .filter(|&m| {
            let s = midi_to_harp_string(m, key_root).unwrap();
            s < thumb_string && thumb_string - s < MAX_HAND_SPAN
        })
        .collect();
    eligible.sort_by(|a, b| b.cmp(a)); // descending MIDI (closest to thumb first)
    eligible.truncate(7); // limit combinatorial explosion

    // 1 additional
    for &a in &eligible {
        configs.push(vec![thumb_midi, a]);
    }

    // 2 additional
    for i in 0..eligible.len() {
        for j in (i + 1)..eligible.len() {
            configs.push(vec![thumb_midi, eligible[i], eligible[j]]);
        }
    }

    // 3 additional (full hand)
    for i in 0..eligible.len() {
        for j in (i + 1)..eligible.len() {
            for k in (j + 1)..eligible.len() {
                let lo_s = midi_to_harp_string(eligible[k], key_root).unwrap();
                if thumb_string - lo_s < MAX_HAND_SPAN {
                    configs.push(vec![
                        thumb_midi,
                        eligible[i],
                        eligible[j],
                        eligible[k],
                    ]);
                }
            }
        }
    }

    configs
}

/// Generate LH hand configurations: 1-4 notes within span, respecting bass requirement.
fn generate_lh_configs(candidates: &[u8], key_root: u8, bass_pc: Option<u8>) -> Vec<Vec<u8>> {
    let mut configs = Vec::new();

    // LH thumb is the highest LH note
    // Try each candidate as LH thumb, then fill below
    for (ti, &thumb) in candidates.iter().enumerate().rev() {
        let thumb_string = match midi_to_harp_string(thumb, key_root) {
            Some(s) => s,
            None => continue,
        };

        let below: Vec<u8> = candidates[..ti]
            .iter()
            .copied()
            .filter(|&m| {
                let s = midi_to_harp_string(m, key_root).unwrap();
                s < thumb_string && thumb_string - s < MAX_HAND_SPAN
            })
            .collect();

        // 1 note (thumb only)
        if bass_ok(&[thumb], bass_pc) {
            configs.push(vec![thumb]);
        }

        // 2 notes
        for &a in &below {
            let hand = vec![thumb, a];
            if bass_ok(&hand, bass_pc) {
                configs.push(hand);
            }
        }

        // 3 notes
        for i in 0..below.len() {
            for j in (i + 1)..below.len() {
                let hand = vec![thumb, below[i], below[j]];
                if bass_ok(&hand, bass_pc) {
                    configs.push(hand);
                }
            }
        }

        // 4 notes
        for i in 0..below.len() {
            for j in (i + 1)..below.len() {
                for k in (j + 1)..below.len() {
                    let hand = vec![thumb, below[i], below[j], below[k]];
                    if bass_ok(&hand, bass_pc) {
                        configs.push(hand);
                    }
                }
            }
        }
    }

    configs
}

/// Check if the lowest note in hand matches the required bass pitch class.
fn bass_ok(hand: &[u8], required_bass_pc: Option<u8>) -> bool {
    match required_bass_pc {
        None => true,
        Some(pc) => {
            if let Some(&lowest) = hand.iter().min() {
                pitch_class(lowest) == pc
            } else {
                false
            }
        }
    }
}

/// Build a Voicing struct from RH and LH MIDI note lists.
fn build_voicing(
    rh: &[u8],
    lh: &[u8],
    key_root: u8,
    root_pc: u8,
    _chord_pcs: &[u8],
) -> Voicing {
    let make_notes = |notes: &[u8], hand: Hand| -> Vec<FingerNote> {
        let mut sorted: Vec<u8> = notes.to_vec();
        sorted.sort_by(|a, b| b.cmp(a)); // high to low
        sorted
            .iter()
            .enumerate()
            .map(|(i, &midi)| FingerNote {
                midi,
                string_num: midi_to_harp_string(midi, key_root).unwrap_or(0),
                finger: (i + 1) as u8, // 1=thumb, 2=index, etc.
                hand,
                note_name: midi_to_name(midi),
                is_root: pitch_class(midi) == root_pc,
            })
            .collect()
    };

    Voicing {
        rh: make_notes(rh, Hand::Right),
        lh: make_notes(lh, Hand::Left),
        score: 0.0,
    }
}

/// Score a voicing for "fullness" — wider spread, complete chord coverage, bass presence.
fn score_voicing(
    voicing: &Voicing,
    root_pc: u8,
    chord_pcs: &[u8],
    _bass_pc: Option<u8>,
) -> f32 {
    let mut score = 0.0f32;
    let midis = voicing.all_midis();
    if midis.is_empty() {
        return 0.0;
    }

    // 1. Total spread
    let span = midis.last().unwrap() - midis.first().unwrap();
    score += span as f32 * 2.0;

    // 2. Bass note is root
    if pitch_class(*midis.first().unwrap()) == root_pc {
        score += 20.0;
    }

    // 3. Root doubled in bass register (below MIDI 60)
    let bass_doublings = midis
        .iter()
        .filter(|&&m| m < 60 && pitch_class(m) == root_pc)
        .count();
    score += bass_doublings as f32 * 10.0;

    // 4. Chord tone coverage
    let covered: HashSet<u8> = midis.iter().map(|&m| pitch_class(m)).collect();
    let coverage = chord_pcs.iter().filter(|&&pc| covered.contains(&pc)).count();
    score += coverage as f32 * 15.0;

    // 5. Total notes — strongly prefer 6-8 notes
    score += voicing.total_notes() as f32 * 8.0;

    // 6. Penalize large gaps
    for w in midis.windows(2) {
        let gap = w[1] - w[0];
        if gap > 12 {
            score -= (gap - 12) as f32 * 3.0;
        }
    }

    // 7. Prefer both hands having notes
    if !voicing.lh.is_empty() && !voicing.rh.is_empty() {
        score += 15.0;
    }

    // 8. Bias load to LH: prefer LH 3-4 notes, RH 2-3 notes
    let lh_count = voicing.lh.len();
    let rh_count = voicing.rh.len();
    if lh_count >= 3 { score += 20.0; }
    if lh_count == 4 { score += 10.0; }
    if rh_count == 2 { score += 10.0; }  // melody + 1 support
    if rh_count == 3 { score += 5.0; }   // melody + 2 support
    if rh_count == 4 && lh_count < 3 { score -= 15.0; } // don't overload RH at LH's expense

    score
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chord_symbol::parse_chord_symbol;

    #[test]
    fn test_basic_voicing() {
        // C major chord, melody = E5 (MIDI 76), key of C
        let chord = parse_chord_symbol("I").unwrap();
        let voicings = enumerate_voicings(76, &chord, 0, 5);
        assert!(!voicings.is_empty(), "Should generate at least one voicing");

        // Check RH thumb is melody note
        let best = &voicings[0];
        assert_eq!(best.rh[0].midi, 76);
        assert_eq!(best.rh[0].finger, 1); // thumb

        // Check no hand overlap
        if let (Some(rh_low), Some(lh_high)) = (best.rh_lowest(), best.lh_highest()) {
            assert!(lh_high < rh_low, "LH must be below RH");
        }
    }

    #[test]
    fn test_span_constraint() {
        let chord = parse_chord_symbol("I").unwrap();
        let voicings = enumerate_voicings(72, &chord, 0, 20);

        for v in &voicings {
            // Check RH span
            if v.rh.len() > 1 {
                let rh_hi = v.rh.first().unwrap().string_num;
                let rh_lo = v.rh.last().unwrap().string_num;
                assert!(
                    rh_hi - rh_lo < MAX_HAND_SPAN,
                    "RH span {} exceeds max {}",
                    rh_hi - rh_lo,
                    MAX_HAND_SPAN
                );
            }
            // Check LH span
            if v.lh.len() > 1 {
                let lh_hi = v.lh.first().unwrap().string_num;
                let lh_lo = v.lh.last().unwrap().string_num;
                assert!(
                    lh_hi - lh_lo < MAX_HAND_SPAN,
                    "LH span {} exceeds max {}",
                    lh_hi - lh_lo,
                    MAX_HAND_SPAN
                );
            }
        }
    }

    #[test]
    fn test_inversion_bass() {
        // V^1 = V chord, 1st inversion (3rd in bass)
        let chord = parse_chord_symbol("V^1").unwrap();
        let voicings = enumerate_voicings(71, &chord, 0, 5);

        if !voicings.is_empty() {
            let best = &voicings[0];
            let all = best.all_midis();
            let bass = *all.first().unwrap();
            // V in C = G-B-D. 1st inversion = B in bass. B = pitch class 11.
            assert_eq!(
                pitch_class(bass),
                11,
                "V^1 should have B in bass, got {}",
                midi_to_name(bass)
            );
        }
    }

    #[test]
    fn test_fullness_scoring() {
        let chord = parse_chord_symbol("I").unwrap();
        let voicings = enumerate_voicings(72, &chord, 0, 10);

        // Best voicing should have more notes than worst
        if voicings.len() >= 2 {
            assert!(voicings[0].score >= voicings[1].score);
            assert!(voicings[0].total_notes() >= 2);
        }
    }
}

#[cfg(test)]
mod debug_tests {
    use super::*;
    use crate::chord_symbol::parse_chord_symbol;

    #[test]
    fn test_debug_voicing_output() {
        // V chord, melody B4 (MIDI 71), key of C
        let chord = parse_chord_symbol("V").unwrap();
        let voicings = enumerate_voicings(71, &chord, 0, 3);
        for (i, v) in voicings.iter().enumerate() {
            println!("Voicing {}: score={:.0}", i, v.score);
            println!("  RH ({} notes):", v.rh.len());
            for n in &v.rh {
                println!("    finger {} midi {} {} root={}", n.finger, n.midi, n.note_name, n.is_root);
            }
            println!("  LH ({} notes):", v.lh.len());
            for n in &v.lh {
                println!("    finger {} midi {} {} root={}", n.finger, n.midi, n.note_name, n.is_root);
            }
        }
        assert!(voicings.len() > 0);
        assert!(voicings[0].rh.len() > 1 || voicings[0].lh.len() > 0, 
            "Best voicing should have more than just melody");
    }
}

#[cfg(test)]
mod bias_tests {
    use super::*;
    use crate::chord_symbol::parse_chord_symbol;

    #[test]
    fn test_lh_bias() {
        // Try several melody positions
        for (chord_str, melody) in &[("I", 60), ("I", 64), ("I", 67), ("V", 71), ("IV", 65)] {
            let chord = parse_chord_symbol(chord_str).unwrap();
            let voicings = enumerate_voicings(*melody, &chord, 0, 1);
            if let Some(v) = voicings.first() {
                println!("{} melody={}: RH={} LH={} total={} score={:.0}",
                    chord_str, melody, v.rh.len(), v.lh.len(), v.total_notes(), v.score);
            } else {
                println!("{} melody={}: NO VOICING", chord_str, melody);
            }
        }
    }
}
