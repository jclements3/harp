use crate::music::*;
use crate::chord::identify_chord;

#[derive(Debug, Clone)]
pub struct Voicing {
    pub rh_midi: Vec<i32>,
    pub lh_midi: Vec<i32>,
    pub rh_strings: Vec<i32>,
    pub lh_strings: Vec<i32>,
    pub rh_chord: String,
    pub lh_chord: String,
    pub full_chord: String,
}

pub fn voice_from_satb(
    soprano: i32,
    alto: Option<i32>,
    tenor: Option<i32>,
    bass: Option<i32>,
    key_root: i32,
    prev: Option<&Voicing>,
) -> Option<Voicing> {
    let melody_midi = snap_to_diatonic(soprano, key_root);
    let melody_str = midi_to_harp_string(melody_midi, key_root)?;

    let mut satb = vec![soprano];
    if let Some(a) = alto { satb.push(a); }
    if let Some(t) = tenor { satb.push(t); }
    if let Some(b) = bass { satb.push(b); }
    if satb.len() < 2 { return None; }

    let satb_snapped: Vec<i32> = satb.iter().map(|&m| snap_to_diatonic(m, key_root)).collect();
    let mut satb_pcs: Vec<i32> = satb_snapped.iter().map(|&m| pitch_class(m)).collect();
    satb_pcs.sort(); satb_pcs.dedup();

    let bass_pc = bass.map(|b| pitch_class(snap_to_diatonic(b, key_root)))
        .unwrap_or(satb_pcs[0]);

    // All candidate notes
    let mut candidates = Vec::new();
    for &pc in &satb_pcs {
        for m in all_octaves(pc, HARP_LOW_MIDI, HARP_HIGH_MIDI) {
            if midi_to_harp_string(m, key_root).is_some() && !candidates.contains(&m) {
                candidates.push(m);
            }
        }
    }
    candidates.sort();

    // RH eligible: below melody, within span
    let mut rh_eligible: Vec<i32> = candidates.iter()
        .filter(|&&m| {
            if m >= melody_midi { return false; }
            let s = match midi_to_harp_string(m, key_root) { Some(s) => s, None => return false };
            melody_str - s < MAX_HAND_SPAN
        })
        .copied().collect();
    rh_eligible.sort_by(|a, b| b.cmp(a));
    rh_eligible.truncate(6);

    // Generate RH configs
    let mut rh_configs: Vec<Vec<i32>> = vec![vec![melody_midi]];
    for i in 0..rh_eligible.len() {
        rh_configs.push(vec![melody_midi, rh_eligible[i]]);
        for j in (i+1)..rh_eligible.len() {
            rh_configs.push(vec![melody_midi, rh_eligible[i], rh_eligible[j]]);
            for k in (j+1)..rh_eligible.len() {
                let lo_s = midi_to_harp_string(rh_eligible[k], key_root);
                if let Some(s) = lo_s {
                    if melody_str - s < MAX_HAND_SPAN {
                        rh_configs.push(vec![melody_midi, rh_eligible[i], rh_eligible[j], rh_eligible[k]]);
                    }
                }
            }
        }
    }

    let mut best_split: Option<(Vec<i32>, Vec<i32>)> = None;
    let mut best_score = i64::MIN;

    for rh in &rh_configs {
        let rh_strs: Vec<Option<i32>> = rh.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
        if rh_strs.iter().any(|s| s.is_none()) { continue; }
        let rh_strs: Vec<i32> = rh_strs.into_iter().map(|s| s.unwrap()).collect();
        let rh_lowest_str = *rh_strs.iter().min().unwrap();

        for gap in 0..=3i32 {
            let lh_upper_str = rh_lowest_str - 1 - gap;
            if lh_upper_str < 0 { continue; }
            let lh_upper_midi = harp_string_to_midi(lh_upper_str, key_root);

            let mut lh_eligible: Vec<i32> = candidates.iter()
                .filter(|&&m| {
                    if m > lh_upper_midi { return false; }
                    let s = match midi_to_harp_string(m, key_root) { Some(s) => s, None => return false };
                    lh_upper_str - s < MAX_HAND_SPAN
                })
                .copied().collect();
            lh_eligible.sort_by(|a, b| b.cmp(a));
            lh_eligible.truncate(6);

            let lh_configs = generate_lh_configs(&lh_eligible, key_root, bass_pc);

            // Also try without bass constraint
            let lh_configs = if lh_configs.is_empty() {
                generate_lh_configs_no_bass(&lh_eligible, key_root)
            } else {
                lh_configs
            };

            for lh in &lh_configs {
                let total = rh.len() + lh.len();
                if total < 5 || total > 8 { continue; }

                // Check no string collisions
                let all_strs: Vec<Option<i32>> = rh.iter().chain(lh.iter())
                    .map(|&m| midi_to_harp_string(m, key_root)).collect();
                if all_strs.iter().any(|s| s.is_none()) { continue; }
                let all_strs: Vec<i32> = all_strs.into_iter().map(|s| s.unwrap()).collect();
                let mut deduped = all_strs.clone(); deduped.sort(); deduped.dedup();
                if deduped.len() != all_strs.len() { continue; }

                let mut all_midis: Vec<i32> = rh.iter().chain(lh.iter()).copied().collect();
                all_midis.sort();

                let mut score = 0i64;

                // R1: LH carries the load
                match lh.len() {
                    4 => score += 30,
                    3 => score += 20,
                    2 => score += 5,
                    _ => {}
                }
                match rh.len() {
                    2 => score += 10,
                    3 => score += 8,
                    _ => {}
                }
                if rh.len() >= 4 && lh.len() < 3 { score -= 15; }
                if total >= 6 && total <= 8 { score += 20; }

                // R2: span
                let rh_span = rh_strs.iter().max().unwrap() - rh_strs.iter().min().unwrap();
                let lh_strs_all: Vec<i32> = lh.iter()
                    .filter_map(|&m| midi_to_harp_string(m, key_root)).collect();
                let lh_span = if lh_strs_all.is_empty() { 0 }
                    else { lh_strs_all.iter().max().unwrap() - lh_strs_all.iter().min().unwrap() };
                score += (rh_span * 3 + lh_span * 3) as i64;

                // R3: melody avoidance
                let melody_pc = pitch_class(melody_midi);
                let rh_support: Vec<i32> = rh.iter().filter(|&&m| m != melody_midi).copied().collect();
                if rh_support.iter().any(|&m| pitch_class(m) == melody_pc) {
                    score -= 15;
                } else {
                    score += 10;
                }

                // R4: chord tone coverage
                let covered: Vec<i32> = all_midis.iter().map(|&m| pitch_class(m)).collect();
                score += satb_pcs.iter().filter(|pc| covered.contains(pc)).count() as i64 * 15;

                // R5: bass
                if !lh.is_empty() && pitch_class(*lh.iter().min().unwrap()) == bass_pc {
                    score += 20;
                }
                score += lh.iter().filter(|&&m| m < 60 && pitch_class(m) == bass_pc).count() as i64 * 10;

                // R6: voice leading
                if let Some(pv) = prev {
                    let prev_all: Vec<i32> = pv.rh_midi.iter().chain(pv.lh_midi.iter()).copied().collect();
                    let movement: i32 = all_midis.iter().map(|&m| {
                        prev_all.iter().map(|&p| (p - m).abs()).min().unwrap_or(20)
                    }).sum();
                    score += (20 - movement).max(0) as i64;
                    if !pv.lh_midi.is_empty() && !lh.is_empty() {
                        let prev_bass = *pv.lh_midi.iter().min().unwrap();
                        let cur_bass = *lh.iter().min().unwrap();
                        let bass_leap = (cur_bass - prev_bass).abs();
                        if bass_leap > 7 { score -= (bass_leap - 7) as i64 * 5; }
                    }
                }

                // R7: gap preference
                score -= (gap - 1).abs() as i64 * 2;

                // R8: even distribution
                for g in 1..all_midis.len() {
                    let d = all_midis[g] - all_midis[g - 1];
                    if d > 12 { score -= (d - 12) as i64 * 5; }
                }

                score += 15; // both hands

                if score > best_score {
                    best_score = score;
                    best_split = Some((rh.clone(), lh.clone()));
                }
            }
        }
    }

    // Progressive relaxation
    if best_split.is_none() {
        'outer: for rh in &rh_configs {
            let rh_strs: Vec<Option<i32>> = rh.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
            if rh_strs.iter().any(|s| s.is_none()) { continue; }
            let rh_strs: Vec<i32> = rh_strs.into_iter().map(|s| s.unwrap()).collect();
            let rh_lowest_str = *rh_strs.iter().min().unwrap();

            for gap in 4..=5i32 {
                let lh_upper_str = rh_lowest_str - 1 - gap;
                if lh_upper_str < 0 { continue; }
                let lh_upper_midi = harp_string_to_midi(lh_upper_str, key_root);
                let mut lh_eligible: Vec<i32> = candidates.iter()
                    .filter(|&&m| {
                        if m > lh_upper_midi { return false; }
                        let s = match midi_to_harp_string(m, key_root) { Some(s) => s, None => return false };
                        lh_upper_str - s < MAX_HAND_SPAN
                    })
                    .copied().collect();
                lh_eligible.sort_by(|a, b| b.cmp(a));
                lh_eligible.truncate(6);

                for ti in 0..lh_eligible.len() {
                    let thumb = lh_eligible[ti];
                    let below = &lh_eligible[ti+1..];
                    let mut configs: Vec<Vec<i32>> = vec![vec![thumb]];
                    for &a in below { configs.push(vec![thumb, a]); }
                    for i in 0..below.len() {
                        for j in (i+1)..below.len() {
                            configs.push(vec![thumb, below[i], below[j]]);
                        }
                    }
                    for lh in &configs {
                        let total = rh.len() + lh.len();
                        if total < 3 { continue; }
                        let all_strs: Vec<Option<i32>> = rh.iter().chain(lh.iter())
                            .map(|&m| midi_to_harp_string(m, key_root)).collect();
                        let all_strs_flat: Vec<i32> = all_strs.iter().filter_map(|&s| s).collect();
                        let mut deduped = all_strs_flat.clone(); deduped.sort(); deduped.dedup();
                        if deduped.len() == all_strs_flat.len() {
                            best_split = Some((rh.clone(), lh.clone()));
                            break 'outer;
                        }
                    }
                }
            }
        }
    }

    // Last resort
    if best_split.is_none() {
        let lh_cands: Vec<i32> = candidates.iter()
            .filter(|&&m| m < melody_midi - 12 && midi_to_harp_string(m, key_root).is_some())
            .copied().collect();
        if lh_cands.len() >= 2 {
            let n = lh_cands.len();
            best_split = Some((vec![melody_midi], vec![lh_cands[n-2], lh_cands[n-1]]));
        } else if !lh_cands.is_empty() {
            best_split = Some((vec![melody_midi], vec![*lh_cands.last().unwrap()]));
        } else {
            best_split = Some((vec![melody_midi], vec![]));
        }
    }

    let (mut rh, mut lh) = best_split?;
    rh.sort_by(|a, b| b.cmp(a));
    lh.sort_by(|a, b| b.cmp(a));

    let rh_strings: Vec<i32> = rh.iter()
        .filter_map(|&m| midi_to_harp_string(m, key_root).map(|s| to_relative_string(s, key_root)))
        .collect();
    let lh_strings: Vec<i32> = lh.iter()
        .filter_map(|&m| midi_to_harp_string(m, key_root).map(|s| to_relative_string(s, key_root)))
        .collect();

    let all_notes: Vec<i32> = rh.iter().chain(lh.iter()).copied().collect();

    Some(Voicing {
        rh_chord: identify_chord(&rh, key_root),
        lh_chord: identify_chord(&lh, key_root),
        full_chord: identify_chord(&all_notes, key_root),
        rh_midi: rh,
        lh_midi: lh,
        rh_strings,
        lh_strings,
    })
}

fn generate_lh_configs(eligible: &[i32], key_root: i32, bass_pc: i32) -> Vec<Vec<i32>> {
    let mut configs = Vec::new();
    for ti in 0..eligible.len() {
        let thumb = eligible[ti];
        let thumb_str = match midi_to_harp_string(thumb, key_root) { Some(s) => s, None => continue };
        let below: Vec<i32> = eligible[ti+1..].iter()
            .filter(|&&m| {
                let s = match midi_to_harp_string(m, key_root) { Some(s) => s, None => return false };
                s < thumb_str && thumb_str - s < MAX_HAND_SPAN
            })
            .copied().collect();

        let check = |h: &[i32]| -> bool {
            pitch_class(*h.iter().min().unwrap()) == bass_pc
        };

        // 2 notes
        for &a in &below {
            let h = vec![thumb, a];
            if check(&h) { configs.push(h); }
        }
        // 3 notes
        for i in 0..below.len() {
            for j in (i+1)..below.len() {
                let h = vec![thumb, below[i], below[j]];
                if check(&h) { configs.push(h); }
            }
        }
        // 4 notes
        for i in 0..below.len() {
            for j in (i+1)..below.len() {
                for k in (j+1)..below.len() {
                    let h = vec![thumb, below[i], below[j], below[k]];
                    let strs: Vec<Option<i32>> = h.iter().map(|&m| midi_to_harp_string(m, key_root)).collect();
                    if strs.iter().any(|s| s.is_none()) { continue; }
                    let strs: Vec<i32> = strs.into_iter().map(|s| s.unwrap()).collect();
                    if strs.iter().max().unwrap() - strs.iter().min().unwrap() >= MAX_HAND_SPAN { continue; }
                    if check(&h) { configs.push(h); }
                }
            }
        }
    }
    configs
}

fn generate_lh_configs_no_bass(eligible: &[i32], key_root: i32) -> Vec<Vec<i32>> {
    let mut configs = Vec::new();
    for ti in 0..eligible.len() {
        let thumb = eligible[ti];
        let thumb_str = match midi_to_harp_string(thumb, key_root) { Some(s) => s, None => continue };
        let below: Vec<i32> = eligible[ti+1..].iter()
            .filter(|&&m| {
                let s = match midi_to_harp_string(m, key_root) { Some(s) => s, None => return false };
                s < thumb_str && thumb_str - s < MAX_HAND_SPAN
            })
            .copied().collect();
        for &a in &below { configs.push(vec![thumb, a]); }
        for i in 0..below.len() {
            for j in (i+1)..below.len() {
                configs.push(vec![thumb, below[i], below[j]]);
            }
        }
    }
    configs
}
