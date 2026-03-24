use std::collections::HashSet;

use serde::Serialize;

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
    let melody_string: i16 = match midi_to_harp_string(melody_midi, key_root) {
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
            if m >= melody_midi { return false; }
            match midi_to_harp_string(m, key_root) {
                Some(s) => melody_string - s < MAX_HAND_SPAN,
                None => false,
            }
        })
        .collect();

    // Generate RH configurations: thumb=melody + 0-3 additional notes
    let rh_configs = generate_hand_configs(&rh_candidates, melody_midi, melody_string, key_root);

    let mut voicings = Vec::new();

    for rh in &rh_configs {
        let rh_lowest_string: i16 = if rh.len() > 1 {
            midi_to_harp_string(rh[rh.len() - 1], key_root).unwrap_or(melody_string)
        } else {
            melody_string
        };

        // LH must be below RH lowest (at least 1 string gap)
        if rh_lowest_string <= 0 {
            continue; // no room for LH
        }
        let lh_upper_bound_string = rh_lowest_string - 1;
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
    thumb_string: i16,
    key_root: u8,
) -> Vec<Vec<u8>> {
    let mut configs = Vec::new();

    // Start with thumb only
    configs.push(vec![thumb_midi]);

    // Pick 1-3 additional notes (sorted descending, limited to closest 7)
    let mut eligible: Vec<u8> = candidates
        .iter()
        .copied()
        .filter(|&m| {
            match midi_to_harp_string(m, key_root) {
                Some(s) => s < thumb_string && thumb_string - s < MAX_HAND_SPAN,
                None => false,
            }
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
                let lo_s = midi_to_harp_string(eligible[k], key_root).unwrap_or(0);
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
        let thumb_string: i16 = match midi_to_harp_string(thumb, key_root) {
            Some(s) => s,
            None => continue,
        };

        let below: Vec<u8> = candidates[..ti]
            .iter()
            .copied()
            .filter(|&m| {
                match midi_to_harp_string(m, key_root) {
                    Some(s) => s < thumb_string && thumb_string - s < MAX_HAND_SPAN,
                    None => false,
                }
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

/// Compute the reference string for relative numbering: lowest root on the harp.
pub fn ref_string(key_root: u8) -> i16 {
    all_octaves(key_root, HARP_LOW_MIDI, HARP_HIGH_MIDI)
        .first()
        .and_then(|&m| midi_to_harp_string(m, key_root))
        .unwrap_or(0)
}

/// Convert a raw harp string to relative skip-zero numbering.
/// String 1 = lowest root on harp. Below that: -1, -2, ... (no zero).
pub fn to_relative_string(raw: i16, ref_str: i16) -> i16 {
    let rel = raw - ref_str + 1;
    if rel <= 0 { rel - 1 } else { rel }
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
                string_num: midi_to_harp_string(midi, key_root).unwrap_or(0i16),
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

// ── SATB Voicing Engine (ported from HarpHymnal.html) ────────────

/// Result of SATB-based voicing: RH and LH MIDI notes plus chord analysis.
#[derive(Debug, Clone, Serialize)]
pub struct SatbVoicing {
    pub rh_midi: Vec<u8>,
    pub lh_midi: Vec<u8>,
    pub rh_strings: Vec<i16>,
    pub lh_strings: Vec<i16>,
    pub rh_chord: String,
    pub lh_chord: String,
    pub full_chord: String,
}

/// Previous voicing state for voice leading.
pub struct PrevVoicing {
    pub rh_midi: Vec<u8>,
    pub lh_midi: Vec<u8>,
}

/// Voice a chord from SATB notes, expanding to 6-8 note harp voicing split across two hands.
/// Implements scoring rules R1-R8 from the harp arranging system.
pub fn voice_from_satb(
    soprano: u8,
    alto: Option<u8>,
    tenor: Option<u8>,
    bass: Option<u8>,
    key_root: u8,
    prev: Option<&PrevVoicing>,
    identify: &dyn Fn(&[u8], u8) -> String,
) -> Option<SatbVoicing> {
    use std::collections::HashSet;

    // Snap melody to diatonic
    let melody_midi = snap_to_diatonic(soprano, key_root);
    let melody_str = midi_to_harp_string(melody_midi, key_root)?;

    // Collect and snap all SATB notes
    let mut satb = vec![soprano];
    if let Some(a) = alto { satb.push(a); }
    if let Some(t) = tenor { satb.push(t); }
    if let Some(b) = bass { satb.push(b); }
    if satb.len() < 2 { return None; }

    let satb_snapped: Vec<u8> = satb.iter().map(|&m| snap_to_diatonic(m, key_root)).collect();
    let satb_pcs: Vec<u8> = {
        let set: HashSet<u8> = satb_snapped.iter().map(|&m| pitch_class(m)).collect();
        set.into_iter().collect()
    };
    let bass_pc = bass.map(|b| pitch_class(snap_to_diatonic(b, key_root)))
        .unwrap_or_else(|| satb_pcs[0]);

    // Generate all candidate notes (chord tones in all octaves, diatonic only)
    let mut all_cands: HashSet<u8> = HashSet::new();
    for &pc in &satb_pcs {
        for m in all_octaves(pc, HARP_LOW_MIDI, HARP_HIGH_MIDI) {
            if midi_to_harp_string(m, key_root).is_some() {
                all_cands.insert(m);
            }
        }
    }
    let mut candidates: Vec<u8> = all_cands.into_iter().collect();
    candidates.sort();

    // RH: melody + up to 3 notes below within 10-string span
    let mut rh_eligible: Vec<u8> = candidates.iter().copied()
        .filter(|&m| {
            if m >= melody_midi { return false; }
            match midi_to_harp_string(m, key_root) {
                Some(s) => melody_str - s < MAX_HAND_SPAN,
                None => false,
            }
        })
        .collect();
    rh_eligible.sort_by(|a, b| b.cmp(a));
    rh_eligible.truncate(6);

    // Generate RH configs: melody + 0-3 below
    let mut rh_configs: Vec<Vec<u8>> = vec![vec![melody_midi]];
    for i in 0..rh_eligible.len() {
        rh_configs.push(vec![melody_midi, rh_eligible[i]]);
        for j in (i + 1)..rh_eligible.len() {
            rh_configs.push(vec![melody_midi, rh_eligible[i], rh_eligible[j]]);
            for k in (j + 1)..rh_eligible.len() {
                let lo_s = midi_to_harp_string(rh_eligible[k], key_root).unwrap_or(0);
                if melody_str - lo_s < MAX_HAND_SPAN {
                    rh_configs.push(vec![melody_midi, rh_eligible[i], rh_eligible[j], rh_eligible[k]]);
                }
            }
        }
    }

    let melody_pc = pitch_class(melody_midi);
    let mut best_split: Option<(Vec<u8>, Vec<u8>)> = None;
    let mut best_score = i32::MIN;

    for rh in &rh_configs {
        let rh_strs: Vec<Option<i16>> = rh.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
        if rh_strs.iter().any(|s| s.is_none()) { continue; }
        let rh_strs: Vec<i16> = rh_strs.into_iter().map(|s| s.unwrap()).collect();
        let rh_lowest_str = *rh_strs.iter().min().unwrap();

        // Try LH with gap 0-3 below RH
        for gap in 0..=3i16 {
            let lh_upper_str = rh_lowest_str - 1 - gap;
            if lh_upper_str < 0 { continue; }
            let lh_upper_midi = harp_string_to_midi(lh_upper_str, key_root);

            let mut lh_eligible: Vec<u8> = candidates.iter().copied()
                .filter(|&m| {
                    if m > lh_upper_midi { return false; }
                    match midi_to_harp_string(m, key_root) {
                        Some(s) => lh_upper_str - s < MAX_HAND_SPAN,
                        None => false,
                    }
                })
                .collect();
            lh_eligible.sort_by(|a, b| b.cmp(a));
            lh_eligible.truncate(6);

            // Generate LH configs: 2-4 notes, prefer bass note as lowest
            let mut lh_configs: Vec<Vec<u8>> = Vec::new();
            for ti in 0..lh_eligible.len() {
                let thumb = lh_eligible[ti];
                let below = &lh_eligible[(ti + 1)..];

                // 2 notes
                for &a in below {
                    let h = vec![thumb, a];
                    if pitch_class(*h.iter().min().unwrap()) == bass_pc {
                        lh_configs.push(h);
                    }
                }
                // 3 notes
                for i in 0..below.len() {
                    for j in (i + 1)..below.len() {
                        let h = vec![thumb, below[i], below[j]];
                        if pitch_class(*h.iter().min().unwrap()) == bass_pc {
                            lh_configs.push(h);
                        }
                    }
                }
                // 4 notes
                for i in 0..below.len() {
                    for j in (i + 1)..below.len() {
                        for k in (j + 1)..below.len() {
                            let h = vec![thumb, below[i], below[j], below[k]];
                            let lh_s: Vec<Option<i16>> = h.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
                            if lh_s.iter().any(|s| s.is_none()) { continue; }
                            let lh_s: Vec<i16> = lh_s.into_iter().map(|s| s.unwrap()).collect();
                            let span = lh_s.iter().max().unwrap() - lh_s.iter().min().unwrap();
                            if span >= MAX_HAND_SPAN { continue; }
                            if pitch_class(*h.iter().min().unwrap()) == bass_pc {
                                lh_configs.push(h);
                            }
                        }
                    }
                }
            }

            // Fallback: without bass constraint
            if lh_configs.is_empty() {
                for ti in 0..lh_eligible.len() {
                    let thumb = lh_eligible[ti];
                    let below = &lh_eligible[(ti + 1)..];
                    for &a in below { lh_configs.push(vec![thumb, a]); }
                    for i in 0..below.len() {
                        for j in (i + 1)..below.len() {
                            lh_configs.push(vec![thumb, below[i], below[j]]);
                        }
                    }
                }
            }

            for lh in &lh_configs {
                let total = rh.len() + lh.len();
                if total < 5 || total > 8 { continue; }

                // Check no string collisions
                let all_strs: Vec<i16> = rh.iter().chain(lh.iter())
                    .filter_map(|&m| midi_to_harp_string(m, key_root))
                    .collect();
                if all_strs.len() != rh.len() + lh.len() { continue; } // some note not diatonic
                let all_strs_set: HashSet<i16> = all_strs.iter().copied().collect();
                if all_strs_set.len() != all_strs.len() { continue; } // duplicate strings

                // 29-string check: highest string + 6 (max mode offset) must fit
                let highest_str = *all_strs.iter().max().unwrap();
                let lowest_str = *all_strs.iter().min().unwrap();
                if highest_str + 6 > 29 + lowest_str { continue; }

                let mut all_midis: Vec<u8> = rh.iter().chain(lh.iter()).copied().collect();
                all_midis.sort();

                // Score — implements SRS rules R1-R8
                let mut score: i32 = 0;

                // R1: LH carries the load
                if lh.len() == 4 { score += 30; }
                else if lh.len() == 3 { score += 20; }
                else if lh.len() == 2 { score += 5; }
                if rh.len() == 2 { score += 10; }
                else if rh.len() == 3 { score += 8; }
                if rh.len() >= 4 && lh.len() < 3 { score -= 15; }
                if total >= 6 && total <= 8 { score += 20; }
                if total < 5 { score -= 20; }

                // R2: Maximize span per hand
                let rh_span = rh_strs.iter().max().unwrap() - rh_strs.iter().min().unwrap();
                let lh_strs_all: Vec<i16> = lh.iter()
                    .filter_map(|&m| midi_to_harp_string(m, key_root))
                    .collect();
                let lh_span = if lh_strs_all.is_empty() { 0 }
                    else { lh_strs_all.iter().max().unwrap() - lh_strs_all.iter().min().unwrap() };
                score += (rh_span * 3) as i32;
                score += (lh_span * 3) as i32;

                // R3: Melody avoidance — don't double melody pitch class in RH support
                let rh_support: Vec<u8> = rh.iter().copied().filter(|&m| m != melody_midi).collect();
                let melody_doubled = rh_support.iter().any(|&m| pitch_class(m) == melody_pc);
                if melody_doubled { score -= 15; } else { score += 10; }

                // R4: Chord tone coverage
                let covered: HashSet<u8> = all_midis.iter().map(|&m| pitch_class(m)).collect();
                score += satb_pcs.iter().filter(|&&pc| covered.contains(&pc)).count() as i32 * 15;

                // R5: Bass correctness
                if pitch_class(*lh.iter().min().unwrap()) == bass_pc { score += 20; }
                score += lh.iter().filter(|&&m| m < 60 && pitch_class(m) == bass_pc).count() as i32 * 10;

                // R6: Voice leading
                if let Some(pv) = prev {
                    let prev_all: Vec<u8> = pv.rh_midi.iter().chain(pv.lh_midi.iter()).copied().collect();
                    if !prev_all.is_empty() {
                        let mut movement: i32 = 0;
                        for &m in &all_midis {
                            let closest = prev_all.iter()
                                .min_by_key(|&&p| (p as i32 - m as i32).unsigned_abs())
                                .unwrap();
                            movement += (*closest as i32 - m as i32).unsigned_abs() as i32;
                        }
                        score += (20 - movement).max(0);

                        // Bass leap penalty
                        let prev_bass = *pv.lh_midi.iter().min().unwrap_or(&60);
                        let cur_bass = *lh.iter().min().unwrap();
                        let bass_leap = (cur_bass as i32 - prev_bass as i32).unsigned_abs() as i32;
                        if bass_leap > 7 { score -= (bass_leap - 7) * 5; }
                    }
                }

                // R7: Gap preference (prefer gap=1 between hands)
                score -= ((gap - 1).unsigned_abs() * 2) as i32;

                // R8: Even distribution — penalize gaps > octave between adjacent notes
                for w in all_midis.windows(2) {
                    let d = w[1] as i32 - w[0] as i32;
                    if d > 12 { score -= (d - 12) * 5; }
                }

                // Both hands used
                score += 15;

                if score > best_score {
                    best_score = score;
                    best_split = Some((rh.clone(), lh.clone()));
                }
            }
        }
    }

    // Progressive relaxation: try wider gap, fewer notes
    if best_split.is_none() {
        'outer: for rh in &rh_configs {
            let rh_strs: Vec<Option<i16>> = rh.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
            if rh_strs.iter().any(|s| s.is_none()) { continue; }
            let rh_lowest_str = *rh_strs.iter().map(|s| s.as_ref().unwrap()).min().unwrap();

            for gap in 4..=5i16 {
                let lh_upper_str = rh_lowest_str - 1 - gap;
                if lh_upper_str < 0 { continue; }
                let lh_upper_midi = harp_string_to_midi(lh_upper_str, key_root);
                let mut lh_eligible: Vec<u8> = candidates.iter().copied()
                    .filter(|&m| {
                        if m > lh_upper_midi { return false; }
                        match midi_to_harp_string(m, key_root) {
                            Some(s) => lh_upper_str - s < MAX_HAND_SPAN,
                            None => false,
                        }
                    })
                    .collect();
                lh_eligible.sort_by(|a, b| b.cmp(a));
                lh_eligible.truncate(6);

                for ti in 0..lh_eligible.len() {
                    let thumb = lh_eligible[ti];
                    let below = &lh_eligible[(ti + 1)..];
                    let mut configs: Vec<Vec<u8>> = vec![vec![thumb]];
                    for &a in below { configs.push(vec![thumb, a]); }
                    for i in 0..below.len() {
                        for j in (i + 1)..below.len() {
                            configs.push(vec![thumb, below[i], below[j]]);
                        }
                    }
                    for lh in &configs {
                        let total = rh.len() + lh.len();
                        if total < 3 { continue; }
                        let all_strs: HashSet<i16> = rh.iter().chain(lh.iter())
                            .filter_map(|&m| midi_to_harp_string(m, key_root))
                            .collect();
                        if all_strs.len() != rh.len() + lh.len() { continue; }
                        best_split = Some((rh.clone(), lh.clone()));
                        break 'outer;
                    }
                }
            }
        }
    }

    // Last resort: melody only in RH, any 2 notes in LH
    if best_split.is_none() {
        let lh_cands: Vec<u8> = candidates.iter().copied()
            .filter(|&m| m + 12 < melody_midi && midi_to_harp_string(m, key_root).is_some())
            .collect();
        if lh_cands.len() >= 2 {
            best_split = Some((vec![melody_midi], lh_cands[lh_cands.len()-2..].to_vec()));
        } else if lh_cands.len() >= 1 {
            best_split = Some((vec![melody_midi], vec![*lh_cands.last().unwrap()]));
        } else {
            best_split = Some((vec![melody_midi], vec![]));
        }
    }

    let (rh, lh) = best_split?;

    // Compute relative string numbers (string 1 = lowest root on harp)
    let rs = ref_string(key_root);

    let to_rel_strings = |notes: &[u8]| -> Vec<i16> {
        let mut sorted = notes.to_vec();
        sorted.sort_by(|a, b| b.cmp(a));
        sorted.iter().filter_map(|&m| {
            midi_to_harp_string(m, key_root).map(|raw| to_relative_string(raw, rs))
        }).collect()
    };

    let all_notes: Vec<u8> = rh.iter().chain(lh.iter()).copied().collect();
    let mut rh_sorted = rh.clone(); rh_sorted.sort_by(|a, b| b.cmp(a));
    let mut lh_sorted = lh.clone(); lh_sorted.sort_by(|a, b| b.cmp(a));

    Some(SatbVoicing {
        rh_strings: to_rel_strings(&rh),
        lh_strings: to_rel_strings(&lh),
        rh_chord: identify(&rh_sorted, key_root),
        lh_chord: identify(&lh_sorted, key_root),
        full_chord: identify(&all_notes, key_root),
        rh_midi: rh_sorted,
        lh_midi: lh_sorted,
    })
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
                let span = rh_hi - rh_lo;
                assert!(
                    span < MAX_HAND_SPAN,
                    "RH span {} exceeds max {}",
                    span,
                    MAX_HAND_SPAN
                );
            }
            // Check LH span
            if v.lh.len() > 1 {
                let lh_hi = v.lh.first().unwrap().string_num;
                let lh_lo = v.lh.last().unwrap().string_num;
                let span = lh_hi - lh_lo;
                assert!(
                    span < MAX_HAND_SPAN,
                    "LH span {} exceeds max {}",
                    span,
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

#[cfg(test)]
mod satb_tests {
    use super::*;

    fn dummy_identify(notes: &[u8], _key_root: u8) -> String {
        format!("chord({})", notes.len())
    }

    #[test]
    fn test_satb_basic_c_major() {
        // C major SATB: S=E5(76), A=C5(72), T=G4(67), B=C3(48)
        let result = voice_from_satb(76, Some(72), Some(67), Some(48), 0, None, &dummy_identify);
        assert!(result.is_some(), "Should produce a voicing for C major SATB");
        let v = result.unwrap();
        assert!(!v.rh_midi.is_empty(), "RH should have notes");
        assert!(!v.lh_midi.is_empty(), "LH should have notes");
        // Melody should be highest RH note
        assert_eq!(v.rh_midi[0], 76, "RH thumb should be melody E5");
        let total = v.rh_midi.len() + v.lh_midi.len();
        assert!(total >= 5 && total <= 8, "Total notes {} should be 5-8", total);
    }

    #[test]
    fn test_satb_key_of_g() {
        // G major: S=B4(71), A=D4(62), T=G3(55), B=G2(43)
        let result = voice_from_satb(71, Some(62), Some(55), Some(43), 7, None, &dummy_identify);
        assert!(result.is_some(), "Should produce a voicing in key of G");
        let v = result.unwrap();
        assert_eq!(v.rh_midi[0], 71);
    }

    #[test]
    fn test_satb_chromatic_snap() {
        // F#4 (66) is chromatic in key of C, should snap to nearest diatonic
        let result = voice_from_satb(66, Some(64), Some(60), Some(48), 0, None, &dummy_identify);
        assert!(result.is_some());
        let v = result.unwrap();
        // Melody should have been snapped to G4 (67) or F4 (65)
        let melody = v.rh_midi[0];
        assert!(melody == 67 || melody == 65, "Melody {} should be snapped diatonic", melody);
    }

    #[test]
    fn test_snap_to_diatonic() {
        // C major: F# (66) should snap to G (67)
        assert_eq!(snap_to_diatonic(66, 0), 67);
        // C major: C (60) stays
        assert_eq!(snap_to_diatonic(60, 0), 60);
        // G major: F# (66) is diatonic, stays
        assert_eq!(snap_to_diatonic(66, 7), 66);
    }
}

#[cfg(test)]
mod string_offset_tests {
    use super::*;

    fn dummy_identify(notes: &[u8], _key_root: u8) -> String { String::new() }

    #[test]
    fn test_relative_strings_skip_zero() {
        // C major SATB: soprano C5(72), alto G4(67), tenor E4(64), bass C3(48)
        let result = voice_from_satb(72, Some(67), Some(64), Some(48), 0, None, &dummy_identify);
        let v = result.expect("Should produce voicing");
        println!("RH midi: {:?}, strings: {:?}", v.rh_midi, v.rh_strings);
        println!("LH midi: {:?}, strings: {:?}", v.lh_midi, v.lh_strings);
        
        // No string should be zero
        for &s in v.rh_strings.iter().chain(v.lh_strings.iter()) {
            assert_ne!(s, 0, "String number must never be zero (skip-zero). Got rh={:?} lh={:?}", 
                v.rh_strings, v.lh_strings);
        }
        
        // Strings should be ordered: RH descending (thumb highest), LH descending
        for w in v.rh_strings.windows(2) {
            assert!(w[0] > w[1], "RH strings should descend: {:?}", v.rh_strings);
        }
        for w in v.lh_strings.windows(2) {
            assert!(w[0] > w[1], "LH strings should descend: {:?}", v.lh_strings);
        }
        
        // All RH strings above all LH strings
        if let (Some(&rh_min), Some(&lh_max)) = (v.rh_strings.last(), v.lh_strings.first()) {
            // Account for skip-zero: -1 is adjacent to 1
            let rh_real = if rh_min < 0 { rh_min + 1 } else { rh_min };
            let lh_real = if lh_max < 0 { lh_max + 1 } else { lh_max };
            assert!(rh_real > lh_real, "RH ({}) must be above LH ({})", rh_min, lh_max);
        }
    }
    
    #[test]
    fn test_string_1_is_lowest_root() {
        // In C major, C1 (midi 24) should map to string 1
        let result = voice_from_satb(72, Some(67), Some(64), Some(48), 0, None, &dummy_identify);
        let v = result.expect("Should produce voicing");
        
        // The lowest root in C major is C1=midi24, whose raw string is midi_to_harp_string(24, 0)
        let raw_c1 = midi_to_harp_string(24, 0).unwrap();
        // ref_str should be this value
        // C3=48 should have raw string, and relative = raw - ref + 1 (with skip-zero)
        let raw_c3 = midi_to_harp_string(48, 0).unwrap();
        let rel_c3 = raw_c3 - raw_c1 + 1; // should be positive, no skip-zero adjustment
        assert!(rel_c3 > 0, "C3 relative string should be positive: {}", rel_c3);
        
        println!("C1 raw string: {}, C3 raw string: {}, C3 relative: {}", raw_c1, raw_c3, rel_c3);
    }

    #[test]
    fn test_skip_zero_in_non_c_key() {
        // Key of G (root=7): lowest root on harp is G1(31), ref_str = string of G1
        // Notes below G1 exist on harp (e.g. C1=24 in G major has degree 3)
        // and should get negative skip-zero strings
        let ref_midi = all_octaves(7, HARP_LOW_MIDI, HARP_HIGH_MIDI)[0]; // G1 = 31
        assert_eq!(ref_midi, 31);
        let ref_str = midi_to_harp_string(31, 7).unwrap();

        // C1 (24) is degree 3 (F=3 in the G major scale... wait, G=0, A=1, B=2, C=3)
        let raw_c1 = midi_to_harp_string(24, 7).unwrap();
        let rel = raw_c1 - ref_str + 1;
        // rel should be <= 0, triggering skip-zero
        assert!(rel <= 0, "C1 in key of G should have rel <= 0, got {}", rel);
        let skip_zero = if rel <= 0 { rel - 1 } else { rel };
        assert_ne!(skip_zero, 0, "Must never be zero after adjustment");
        assert!(skip_zero < 0, "Below-root note should be negative: {}", skip_zero);
        println!("Key of G: C1 raw={}, ref_str={}, rel={}, skip_zero={}", raw_c1, ref_str, rel, skip_zero);
    }
}
