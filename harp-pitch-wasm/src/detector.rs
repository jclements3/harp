use rustfft::{num_complex::Complex, FftPlanner};
use std::f32::consts::PI;

/// Harmonic weights for scoring: fundamental, 2nd, 3rd, 4th partial
const HARMONIC_WEIGHTS: [f32; 4] = [1.0, 0.7, 0.4, 0.25];
/// Number of harmonics to zero out during subtraction
const SUBTRACTION_HARMONICS: usize = 6;
/// Maximum number of simultaneous notes to detect
const MAX_NOTES: usize = 4;

pub struct Detector {
    fft_size: usize,
    sample_rate: f32,
    hann_window: Vec<f32>,
    fft_scratch: Vec<Complex<f32>>,
    fft: std::sync::Arc<dyn rustfft::Fft<f32>>,
    sensitivity: f32,
}

impl Detector {
    pub fn new(fft_size: usize, sample_rate: u32) -> Self {
        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(fft_size);
        let scratch_len = fft.get_inplace_scratch_len();

        let hann_window: Vec<f32> = (0..fft_size)
            .map(|i| 0.5 * (1.0 - (2.0 * PI * i as f32 / fft_size as f32).cos()))
            .collect();

        Detector {
            fft_size,
            sample_rate: sample_rate as f32,
            hann_window,
            fft_scratch: vec![Complex::new(0.0, 0.0); scratch_len],
            fft,
            sensitivity: 1.0,
        }
    }

    pub fn set_sensitivity(&mut self, sensitivity: f32) {
        self.sensitivity = sensitivity.clamp(0.1, 5.0);
    }

    pub fn detect_notes(
        &mut self,
        samples: &[f32],
        available_pitch_classes: &[u8],
        low_midi: u8,
        high_midi: u8,
    ) -> Vec<u8> {
        let n = samples.len().min(self.fft_size);

        // 1. Noise gate — RMS check
        let rms = (samples[..n].iter().map(|&s| s * s).sum::<f32>() / n as f32).sqrt();
        if rms < 0.005 * self.sensitivity {
            return Vec::new();
        }

        // 2. Apply Hann window and convert to complex
        let mut buffer: Vec<Complex<f32>> = (0..self.fft_size)
            .map(|i| {
                if i < n {
                    Complex::new(samples[i] * self.hann_window[i], 0.0)
                } else {
                    Complex::new(0.0, 0.0)
                }
            })
            .collect();

        // 3. FFT
        self.fft
            .process_with_scratch(&mut buffer, &mut self.fft_scratch);

        // 4. Compute magnitude spectrum (only need first half)
        let half = self.fft_size / 2;
        let mut magnitudes: Vec<f32> = buffer[..half].iter().map(|c| c.norm()).collect();

        // 5. Adaptive noise floor — median magnitude × threshold
        let noise_floor = median_magnitude(&magnitudes) * 3.0 / self.sensitivity;

        // 6. Build candidate list from pedal state
        let candidates = build_candidates(available_pitch_classes, low_midi, high_midi);
        if candidates.is_empty() {
            return Vec::new();
        }

        // 7. Greedy harmonic subtraction
        let mut detected: Vec<u8> = Vec::new();
        let bin_freq = self.sample_rate / self.fft_size as f32;
        let absolute_threshold = noise_floor * 2.0;
        let mut first_score = 0.0f32;

        for iteration in 0..MAX_NOTES {
            // Score each remaining candidate against current (residual) spectrum
            let mut best_midi = 0u8;
            let mut best_score = 0.0f32;

            for &midi in &candidates {
                if detected.contains(&midi) {
                    continue;
                }
                let f0 = midi_to_freq(midi);
                let score =
                    score_candidate(f0, &magnitudes, bin_freq, half);
                if score > best_score {
                    best_score = score;
                    best_midi = midi;
                }
            }

            // Absolute threshold: must exceed noise floor
            if best_score <= absolute_threshold {
                break;
            }

            // Relative threshold: subsequent notes must be at least 10% of the
            // strongest note's score to avoid detecting harmonic residue artifacts
            if iteration == 0 {
                first_score = best_score;
            } else if best_score < first_score * 0.10 {
                break;
            }

            detected.push(best_midi);

            // Subtract harmonics of detected note from residual spectrum
            let f0 = midi_to_freq(best_midi);
            subtract_harmonics(f0, &mut magnitudes, bin_freq, half);
        }

        detected.sort();
        detected
    }
}

/// Score a candidate fundamental frequency against the magnitude spectrum.
/// Uses parabolic interpolation around each harmonic bin for better accuracy.
fn score_candidate(f0: f32, magnitudes: &[f32], bin_freq: f32, half: usize) -> f32 {
    let mut score = 0.0f32;

    for (h, &weight) in HARMONIC_WEIGHTS.iter().enumerate() {
        let harmonic_freq = f0 * (h + 1) as f32;
        let bin = harmonic_freq / bin_freq;
        let bin_idx = bin.round() as usize;

        if bin_idx == 0 || bin_idx + 1 >= half {
            continue;
        }

        // Parabolic interpolation for better peak estimation
        let alpha = magnitudes[bin_idx - 1];
        let beta = magnitudes[bin_idx];
        let gamma = magnitudes[bin_idx + 1];

        // Peak magnitude from parabolic fit
        let peak = beta - 0.25 * (alpha - gamma) * (alpha - gamma)
            / (alpha - 2.0 * beta + gamma + 1e-10);
        let peak = peak.max(beta); // at least as good as the center bin

        score += peak * weight;
    }

    score
}

/// Attenuate bins around each harmonic of the given fundamental.
/// Uses partial subtraction (70% attenuation) rather than zeroing, so that
/// energy from other notes sharing harmonic bins is preserved.
fn subtract_harmonics(f0: f32, magnitudes: &mut [f32], bin_freq: f32, half: usize) {
    for h in 0..SUBTRACTION_HARMONICS {
        let harmonic_freq = f0 * (h + 1) as f32;
        let center_bin = (harmonic_freq / bin_freq).round() as i32;
        // Attenuate a window of bins around the harmonic
        let width = 2 + h as i32;
        for b in (center_bin - width)..=(center_bin + width) {
            if b >= 0 && (b as usize) < half {
                let dist = (b - center_bin).unsigned_abs() as f32;
                // Stronger attenuation at center, weaker at edges
                let factor = if dist < 1.0 { 0.15 } else { 0.3 };
                magnitudes[b as usize] *= factor;
            }
        }
    }
}

/// Convert MIDI note number to frequency in Hz.
fn midi_to_freq(midi: u8) -> f32 {
    440.0 * 2.0f32.powf((midi as f32 - 69.0) / 12.0)
}

/// Build the set of candidate MIDI notes from available pitch classes within range.
fn build_candidates(available_pitch_classes: &[u8], low_midi: u8, high_midi: u8) -> Vec<u8> {
    let mut candidates = Vec::new();
    for midi in low_midi..=high_midi {
        let pc = midi % 12;
        if available_pitch_classes.iter().any(|&apc| apc == pc) {
            candidates.push(midi);
        }
    }
    candidates
}

/// Compute the median of a magnitude slice for adaptive noise floor estimation.
fn median_magnitude(mags: &[f32]) -> f32 {
    if mags.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<f32> = mags.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    const SAMPLE_RATE: u32 = 44100;
    const FFT_SIZE: usize = 8192;

    /// Generate a harp-like tone: fundamental + decaying harmonics
    fn harp_tone(freq: f32, amplitude: f32, samples: usize, sr: f32) -> Vec<f32> {
        let mut buf = vec![0.0f32; samples];
        // Harp-like profile: strong fundamental, decaying harmonics
        let harmonic_amps = [1.0, 0.5, 0.25, 0.12, 0.06, 0.03];
        for (h, &amp) in harmonic_amps.iter().enumerate() {
            let f = freq * (h + 1) as f32;
            for (i, sample) in buf.iter_mut().enumerate() {
                *sample += amplitude * amp * (2.0 * PI * f * i as f32 / sr).sin();
            }
        }
        buf
    }

    /// Mix multiple signals together
    fn mix(signals: &[Vec<f32>]) -> Vec<f32> {
        let len = signals.iter().map(|s| s.len()).max().unwrap_or(0);
        let mut result = vec![0.0f32; len];
        for sig in signals {
            for (i, &s) in sig.iter().enumerate() {
                result[i] += s;
            }
        }
        result
    }

    fn midi_freq(midi: u8) -> f32 {
        440.0 * 2.0f32.powf((midi as f32 - 69.0) / 12.0)
    }

    /// All 7 natural pitch classes (C major pedal setting)
    fn c_major_pcs() -> Vec<u8> {
        vec![0, 2, 4, 5, 7, 9, 11]
    }

    #[test]
    fn test_single_note_c4() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        let signal = harp_tone(midi_freq(60), 0.5, FFT_SIZE, SAMPLE_RATE as f32);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.contains(&60), "Should detect C4 (MIDI 60), got {:?}", result);
    }

    #[test]
    fn test_two_notes_c4_e4() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        let c4 = harp_tone(midi_freq(60), 0.4, FFT_SIZE, SAMPLE_RATE as f32);
        let e4 = harp_tone(midi_freq(64), 0.4, FFT_SIZE, SAMPLE_RATE as f32);
        let signal = mix(&[c4, e4]);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.contains(&60), "Should detect C4, got {:?}", result);
        assert!(result.contains(&64), "Should detect E4, got {:?}", result);
    }

    #[test]
    fn test_three_note_chord_c_major() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        let c4 = harp_tone(midi_freq(60), 0.35, FFT_SIZE, SAMPLE_RATE as f32);
        let e4 = harp_tone(midi_freq(64), 0.35, FFT_SIZE, SAMPLE_RATE as f32);
        let g4 = harp_tone(midi_freq(67), 0.35, FFT_SIZE, SAMPLE_RATE as f32);
        let signal = mix(&[c4, e4, g4]);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.contains(&60), "Should detect C4, got {:?}", result);
        assert!(result.contains(&64), "Should detect E4, got {:?}", result);
        assert!(result.contains(&67), "Should detect G4, got {:?}", result);
    }

    #[test]
    fn test_four_note_hand_position() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        let c4 = harp_tone(midi_freq(60), 0.3, FFT_SIZE, SAMPLE_RATE as f32);
        let e4 = harp_tone(midi_freq(64), 0.3, FFT_SIZE, SAMPLE_RATE as f32);
        let g4 = harp_tone(midi_freq(67), 0.3, FFT_SIZE, SAMPLE_RATE as f32);
        let c5 = harp_tone(midi_freq(72), 0.3, FFT_SIZE, SAMPLE_RATE as f32);
        let signal = mix(&[c4, e4, g4, c5]);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.contains(&60), "Should detect C4, got {:?}", result);
        assert!(result.contains(&64), "Should detect E4, got {:?}", result);
        assert!(result.contains(&67), "Should detect G4, got {:?}", result);
        assert!(result.contains(&72), "Should detect C5, got {:?}", result);
    }

    #[test]
    fn test_silence_returns_empty() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        let signal = vec![0.0f32; FFT_SIZE];
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.is_empty(), "Silence should return empty, got {:?}", result);
    }

    #[test]
    fn test_pedal_filtering() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        // Play F#4 (MIDI 66) but only allow C major pitch classes (no F# = pc 6)
        let signal = harp_tone(midi_freq(66), 0.5, FFT_SIZE, SAMPLE_RATE as f32);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        // Should NOT detect 66 since pitch class 6 is not in C major
        assert!(!result.contains(&66), "Should not detect F#4 with C major pedals, got {:?}", result);
    }

    #[test]
    fn test_harmonic_disambiguation_c3_not_c4() {
        let mut det = Detector::new(FFT_SIZE, SAMPLE_RATE);
        // Play only C3 (MIDI 48) — should not falsely detect C4 (60) from the 2nd harmonic
        let signal = harp_tone(midi_freq(48), 0.5, FFT_SIZE, SAMPLE_RATE as f32);
        let result = det.detect_notes(&signal, &c_major_pcs(), 48, 79);
        assert!(result.contains(&48), "Should detect C3 (48), got {:?}", result);
        // C4 should not appear — the subtraction of C3's harmonics should prevent it
        assert!(!result.contains(&60), "Should NOT detect C4 (60) from C3's harmonics, got {:?}", result);
    }
}
