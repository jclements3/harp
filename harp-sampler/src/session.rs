use serde::{Deserialize, Serialize};
use crate::audio;

// ── Lever harp layout ──
// Typical 34-string lever harp: C1 to A6
// Each string is a diatonic note; levers raise by a half step.
// We sample each string natural + sharped (lever engaged).

/// What kind of sampling session to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SessionType {
    /// Every string, natural then sharped — full chromatic catalog
    Chromatic,
    /// Common chords in several keys
    Chords,
    /// Free recording with manual labeling
    Free,
}

/// A single note/chord to be sampled.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleTarget {
    /// Human-readable label (e.g., "C4", "C4# (lever)", "C major triad")
    pub label: String,
    /// Expected MIDI note(s)
    pub midi_notes: Vec<u8>,
    /// Whether lever is engaged (for single-note targets)
    pub lever: bool,
    /// Category for organization
    pub category: String,
}

/// Metadata saved alongside each WAV file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleMetadata {
    pub target: SampleTarget,
    pub wav_filename: String,
    pub sample_rate: u32,
    pub duration_secs: f32,
    pub rms_volume: f32,
    pub spectral_centroid_hz: f32,
    pub peak_frequencies_hz: Vec<f32>,
    pub peak_magnitudes_db: Vec<f32>,
    pub timestamp: String,
}

/// State of the sampling session.
pub struct Session {
    pub session_type: SessionType,
    pub targets: Vec<SampleTarget>,
    pub current_idx: usize,
    pub completed: Vec<SampleMetadata>,
    pub output_dir: std::path::PathBuf,
    pub session_name: String,
}

impl Session {
    pub fn new(session_type: SessionType, output_dir: std::path::PathBuf) -> Self {
        let session_name = format!(
            "{}_{:?}",
            chrono::Local::now().format("%Y%m%d_%H%M%S"),
            session_type
        );
        let dir = output_dir.join(&session_name);
        let targets = match session_type {
            SessionType::Chromatic => build_chromatic_targets(),
            SessionType::Chords => build_chord_targets(),
            SessionType::Free => Vec::new(),
        };
        Self {
            session_type,
            targets,
            current_idx: 0,
            completed: Vec::new(),
            output_dir: dir,
            session_name,
        }
    }

    pub fn current_target(&self) -> Option<&SampleTarget> {
        self.targets.get(self.current_idx)
    }

    pub fn is_complete(&self) -> bool {
        self.current_idx >= self.targets.len() && self.session_type != SessionType::Free
    }

    pub fn progress(&self) -> (usize, usize) {
        (self.completed.len(), self.targets.len())
    }

    /// Save a recorded sample, compute analysis, advance to next target.
    pub fn save_sample(
        &mut self,
        samples: &[f32],
        sample_rate: u32,
        rms: f32,
        spectrum: &[f32],
        label_override: Option<String>,
    ) -> Result<SampleMetadata, String> {
        std::fs::create_dir_all(&self.output_dir)
            .map_err(|e| format!("Can't create output dir: {e}"))?;

        let target = if self.session_type == SessionType::Free {
            SampleTarget {
                label: label_override.unwrap_or_else(|| {
                    format!("free_{:03}", self.completed.len() + 1)
                }),
                midi_notes: Vec::new(),
                lever: false,
                category: "free".into(),
            }
        } else {
            self.current_target()
                .cloned()
                .ok_or("No more targets")?
        };

        let filename = format!(
            "{:03}_{}.wav",
            self.completed.len() + 1,
            sanitize_filename(&target.label)
        );
        let wav_path = self.output_dir.join(&filename);

        audio::write_wav(&wav_path, samples, sample_rate)
            .map_err(|e| format!("WAV write failed: {e}"))?;

        let centroid = audio::spectral_centroid(spectrum, sample_rate);
        let peaks = audio::find_peaks(spectrum, sample_rate, 8);

        let meta = SampleMetadata {
            target,
            wav_filename: filename,
            sample_rate,
            duration_secs: samples.len() as f32 / sample_rate as f32,
            rms_volume: rms,
            spectral_centroid_hz: centroid,
            peak_frequencies_hz: peaks.iter().map(|p| p.0).collect(),
            peak_magnitudes_db: peaks.iter().map(|p| p.1).collect(),
            timestamp: chrono::Local::now().to_rfc3339(),
        };

        self.completed.push(meta.clone());

        // Save metadata JSON after each sample (incremental)
        let json_path = self.output_dir.join("samples.json");
        let json = serde_json::to_string_pretty(&self.completed).unwrap();
        std::fs::write(&json_path, json).ok();

        if self.session_type != SessionType::Free {
            self.current_idx += 1;
        }

        Ok(meta)
    }

    /// Skip the current target without recording.
    pub fn skip(&mut self) {
        if self.current_idx < self.targets.len() {
            self.current_idx += 1;
        }
    }
}

/// Build targets for every lever harp string, natural + sharped.
/// 34-string harp: C1, D1, E1, F1, G1, A1, B1, ... up to A6.
fn build_chromatic_targets() -> Vec<SampleTarget> {
    let mut targets = Vec::new();

    // Diatonic notes on a 34-string lever harp (C1 up to A6)
    // String numbering: string 1 = highest (A6), string 34 = lowest (C1)
    // MIDI: C1=24, D1=26, E1=28, F1=29, G1=31, A1=33, B1=35, C2=36, ...
    let diatonic_semitones: [u8; 7] = [0, 2, 4, 5, 7, 9, 11]; // C D E F G A B
    let diatonic_names: [&str; 7] = ["C", "D", "E", "F", "G", "A", "B"];

    // Build all strings: 5 full octaves (C1-B5) + partial octave (C6-A6)
    let mut strings: Vec<(u8, &str, u8)> = Vec::new(); // (midi, name, octave)

    for octave in 1..=5u8 {
        let base_midi = (octave + 1) * 12; // C1=24 when octave=1
        for (i, &semi) in diatonic_semitones.iter().enumerate() {
            strings.push((base_midi + semi, diatonic_names[i], octave));
        }
    }
    // Octave 6: C6, D6, E6, F6, G6, A6
    let base6 = 7 * 12; // 84 = C6
    for i in 0..6 {
        strings.push((base6 + diatonic_semitones[i], diatonic_names[i], 6));
    }

    // For each string: natural, then sharped (lever up = +1 semitone)
    for (midi, name, octave) in &strings {
        // Natural
        targets.push(SampleTarget {
            label: format!("{}{}", name, octave),
            midi_notes: vec![*midi],
            lever: false,
            category: "single_natural".into(),
        });
        // Sharped (lever engaged)
        let sharped_midi = midi + 1;
        let sharped_name = audio::midi_to_name(sharped_midi);
        targets.push(SampleTarget {
            label: format!("{} (lever)", sharped_name),
            midi_notes: vec![sharped_midi],
            lever: true,
            category: "single_lever".into(),
        });
    }

    targets
}

/// Build targets for common chords across several keys.
fn build_chord_targets() -> Vec<SampleTarget> {
    let mut targets = Vec::new();

    // Major, minor, and dim triads in keys commonly used on lever harp
    let keys: [(u8, &str); 8] = [
        (60, "C"), (62, "D"), (64, "E"), (65, "F"),
        (67, "G"), (69, "A"), (71, "B"), (63, "Eb"),
    ];

    let chord_types: [(&str, &[u8]); 4] = [
        ("major", &[0, 4, 7]),
        ("minor", &[0, 3, 7]),
        ("dom7",  &[0, 4, 7, 10]),
        ("dim",   &[0, 3, 6]),
    ];

    for (root_midi, root_name) in &keys {
        for (chord_name, intervals) in &chord_types {
            let midi_notes: Vec<u8> = intervals.iter().map(|&i| root_midi + i).collect();
            let note_names: Vec<String> = midi_notes.iter().map(|&m| audio::midi_to_name(m)).collect();
            targets.push(SampleTarget {
                label: format!("{} {} ({})", root_name, chord_name, note_names.join("-")),
                midi_notes,
                lever: false,
                category: format!("chord_{}", chord_name),
            });
        }
    }

    targets
}

fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            'a'..='z' | 'A'..='Z' | '0'..='9' | '-' | '_' => c,
            ' ' | '(' | ')' => '_',
            '#' => 's',
            _ => '_',
        })
        .collect()
}
