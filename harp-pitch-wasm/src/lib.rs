mod detector;

use detector::Detector;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct HarpPitchDetector {
    inner: Detector,
}

#[wasm_bindgen]
impl HarpPitchDetector {
    #[wasm_bindgen(constructor)]
    pub fn new(fft_size: usize, sample_rate: u32) -> Self {
        HarpPitchDetector {
            inner: Detector::new(fft_size, sample_rate),
        }
    }

    /// Detect polyphonic notes from audio samples.
    ///
    /// - `samples`: Float32Array of time-domain audio (length = fft_size)
    /// - `available_pitch_classes`: Uint8Array of 7 pitch class values (0-11)
    ///   representing the current harp pedal state
    /// - `low_midi` / `high_midi`: MIDI range to search (typically 48-79)
    ///
    /// Returns: Uint8Array of detected MIDI notes (0-4 notes), sorted ascending.
    pub fn detect_notes(
        &mut self,
        samples: &[f32],
        available_pitch_classes: &[u8],
        low_midi: u8,
        high_midi: u8,
    ) -> Vec<u8> {
        self.inner
            .detect_notes(samples, available_pitch_classes, low_midi, high_midi)
    }

    pub fn set_sensitivity(&mut self, sensitivity: f32) {
        self.inner.set_sensitivity(sensitivity);
    }
}
