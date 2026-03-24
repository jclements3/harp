use harp_voicing_wasm::abc::{Beat, Hymn};
use harp_voicing_wasm::music::harp_string_to_midi;
use harp_voicing_wasm::voicing::to_relative_string;
use midly::num::{u4, u7, u15, u24, u28};
use midly::{Format, Header, MetaMessage, MidiMessage, Smf, TrackEvent, TrackEventKind};

const TICKS_PER_QUARTER: u16 = 480;
const HARP_PROGRAM: u7 = u7::new(46); // GM Orchestral Harp

/// Convert a duration rational (num, den) in quarter-note beats to MIDI ticks.
fn duration_to_ticks(dur: (u32, u32)) -> u32 {
    (TICKS_PER_QUARTER as u32 * dur.0) / dur.1
}

/// Convert a BPM to MIDI microseconds per quarter note.
fn bpm_to_tempo(bpm: u32) -> u32 {
    60_000_000 / bpm.max(1)
}

/// Parse a meter string like "4/4" into (numerator, denominator).
fn parse_meter(meter: &str) -> (u8, u8) {
    let parts: Vec<&str> = meter.split('/').collect();
    if parts.len() == 2 {
        let num = parts[0].trim().parse().unwrap_or(4);
        let den = parts[1].trim().parse().unwrap_or(4);
        (num, den)
    } else {
        (4, 4)
    }
}

/// Apply mode offset to MIDI notes via harp string transposition.
fn transpose_for_mode(midi_notes: &[u8], key_root: u8, mode_offset: i16, ref_str: i16) -> Vec<u8> {
    midi_notes
        .iter()
        .filter_map(|&midi| {
            let raw_str = harp_voicing_wasm::music::midi_to_harp_string(midi, key_root)?;
            let rel_str = to_relative_string(raw_str, ref_str);
            let shifted_rel = rel_str + mode_offset;
            // Undo skip-zero
            let shifted_raw = if shifted_rel >= 1 {
                shifted_rel + ref_str - 1
            } else {
                shifted_rel + 1 + ref_str - 1
            };
            let new_midi = harp_string_to_midi(shifted_raw, key_root);
            if new_midi >= 24 && new_midi <= 103 { Some(new_midi) } else { None }
        })
        .collect()
}

const MODE_NAMES: [&str; 7] = [
    "Ionian", "Dorian", "Phrygian", "Lydian", "Mixolydian", "Aeolian", "Locrian",
];

#[derive(Clone, Copy)]
enum Hand { Melody, Right, Left }

/// Generate a complete MIDI file for one hymn in one mode.
pub fn hymn_to_midi(hymn: &Hymn, mode_offset: i16, velocity: u8) -> Vec<u8> {
    let tempo_us = bpm_to_tempo(hymn.tempo_bpm);
    let (meter_num, meter_den) = parse_meter(&hymn.meter);
    let ref_str = harp_voicing_wasm::voicing::ref_string(hymn.key_root);
    let mode_name = MODE_NAMES.get(mode_offset as usize).unwrap_or(&"Unknown");
    let track_name = format!("{} - {} ({})", hymn.number, hymn.title, mode_name);
    let track_name_bytes = track_name.into_bytes();

    // Pre-allocate owned name bytes for each track
    let melody_name = b"Melody".to_vec();
    let rh_name = b"Right Hand".to_vec();
    let lh_name = b"Left Hand".to_vec();

    let vel = u7::new(velocity.min(127));
    let den_power = (meter_den as f64).log2() as u8;

    // -- Build all events as raw bytes using midly's write --
    let mut smf = Smf::new(Header::new(
        Format::Parallel,
        midly::Timing::Metrical(u15::new(TICKS_PER_QUARTER)),
    ));

    // Track 0: Conductor
    {
        let mut t = Vec::new();
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Meta(MetaMessage::TrackName(&track_name_bytes)),
        });
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Meta(MetaMessage::Tempo(u24::new(tempo_us))),
        });
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Meta(MetaMessage::TimeSignature(meter_num, den_power, 24, 8)),
        });
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Meta(MetaMessage::EndOfTrack),
        });
        smf.tracks.push(t);
    }

    // Tracks 1-3: Melody, RH, LH
    for (hand, name_bytes, ch) in [
        (Hand::Melody, &melody_name, u4::new(0)),
        (Hand::Right, &rh_name, u4::new(1)),
        (Hand::Left, &lh_name, u4::new(2)),
    ] {
        let mut t = Vec::new();

        // Track name + program change
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Meta(MetaMessage::TrackName(name_bytes)),
        });
        t.push(TrackEvent {
            delta: u28::new(0),
            kind: TrackEventKind::Midi {
                channel: ch,
                message: MidiMessage::ProgramChange { program: HARP_PROGRAM },
            },
        });

        let mut last_rh: Vec<u8> = Vec::new();
        let mut last_lh: Vec<u8> = Vec::new();
        let mut pending: u32 = 0;

        for beat in &hymn.beats {
            match beat.beat_type.as_str() {
                "bar" => {}
                "rest" => {
                    if let Some(dur) = beat.duration {
                        pending += duration_to_ticks(dur);
                    }
                }
                "note" => {
                    let ticks = beat.duration.map(duration_to_ticks).unwrap_or(TICKS_PER_QUARTER as u32);

                    if let Some(ref rh) = beat.rh_midi { last_rh = rh.clone(); }
                    if let Some(ref lh) = beat.lh_midi { last_lh = lh.clone(); }

                    let src = match hand {
                        Hand::Melody => beat.midi.map(|m| vec![m]).unwrap_or_default(),
                        Hand::Right => last_rh.clone(),
                        Hand::Left => last_lh.clone(),
                    };

                    let notes = if mode_offset == 0 {
                        src
                    } else {
                        transpose_for_mode(&src, hymn.key_root, mode_offset, ref_str)
                    };

                    if notes.is_empty() {
                        pending += ticks;
                        continue;
                    }

                    // Note-on
                    for (i, &midi) in notes.iter().enumerate() {
                        t.push(TrackEvent {
                            delta: u28::new(if i == 0 { pending } else { 0 }),
                            kind: TrackEventKind::Midi {
                                channel: ch,
                                message: MidiMessage::NoteOn {
                                    key: u7::new(midi.min(127)),
                                    vel,
                                },
                            },
                        });
                    }

                    // Note-off
                    for (i, &midi) in notes.iter().enumerate() {
                        t.push(TrackEvent {
                            delta: u28::new(if i == 0 { ticks } else { 0 }),
                            kind: TrackEventKind::Midi {
                                channel: ch,
                                message: MidiMessage::NoteOff {
                                    key: u7::new(midi.min(127)),
                                    vel: u7::new(0),
                                },
                            },
                        });
                    }

                    pending = 0;
                }
                _ => {}
            }
        }

        t.push(TrackEvent {
            delta: u28::new(pending),
            kind: TrackEventKind::Meta(MetaMessage::EndOfTrack),
        });

        smf.tracks.push(t);
    }

    let mut buf = Vec::new();
    smf.write(&mut buf).expect("Failed to write MIDI");
    buf
}
