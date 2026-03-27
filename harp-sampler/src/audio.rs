use std::sync::{Arc, Mutex};

use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::Stream;
use rustfft::{FftPlanner, num_complex::Complex};

// ── Shared state between audio thread and UI ──

pub struct AudioState {
    /// Ring buffer of recent raw samples for waveform display
    pub waveform: Vec<f32>,
    /// Most recent magnitude spectrum (linear, 0..nyquist)
    pub spectrum: Vec<f32>,
    /// RMS volume level (0.0-1.0)
    pub volume: f32,
    /// Whether the audio stream is active
    pub running: bool,
    /// Device sample rate
    pub sample_rate: u32,
    /// Accumulator for recording (when armed)
    pub recording: bool,
    pub recorded_samples: Vec<f32>,
}

impl Default for AudioState {
    fn default() -> Self {
        Self {
            waveform: Vec::new(),
            spectrum: vec![0.0; 1024],
            volume: 0.0,
            running: false,
            sample_rate: 44100,
            recording: false,
            recorded_samples: Vec::new(),
        }
    }
}

/// Compute magnitude spectrum from a buffer of f32 samples.
/// Returns magnitudes for bins 0..n/2 (DC to Nyquist), dB scale.
pub fn compute_spectrum(samples: &[f32], fft_size: usize) -> Vec<f32> {
    let n = fft_size.min(samples.len());
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);

    // Apply Hann window
    let mut buffer: Vec<Complex<f32>> = samples[..n]
        .iter()
        .enumerate()
        .map(|(i, &s)| {
            let w = 0.5 * (1.0 - (2.0 * std::f32::consts::PI * i as f32 / n as f32).cos());
            Complex::new(s * w, 0.0)
        })
        .collect();

    fft.process(&mut buffer);

    let half = n / 2;
    buffer[..half]
        .iter()
        .map(|c| {
            let mag = c.norm() / half as f32;
            // dB scale, floor at -80 dB
            (20.0 * mag.max(1e-8).log10()).max(-80.0)
        })
        .collect()
}

/// Start capturing audio from the default input device.
pub fn start_capture(state: Arc<Mutex<AudioState>>) -> Result<Stream, Box<dyn std::error::Error>> {
    let host = cpal::default_host();
    let device = host
        .default_input_device()
        .ok_or("no input device available")?;

    let config = device.default_input_config()?;
    let sample_rate = config.sample_rate().0;

    {
        let mut s = state.lock().unwrap();
        s.sample_rate = sample_rate;
        s.running = true;
    }

    let channels = config.channels() as usize;
    let waveform_len = 2048; // ~46ms at 44.1kHz
    let fft_size = 4096;
    let mut accumulator: Vec<f32> = Vec::with_capacity(fft_size);

    let stream = device.build_input_stream(
        &config.into(),
        move |data: &[f32], _: &cpal::InputCallbackInfo| {
            // Downmix to mono
            let mono: Vec<f32> = data
                .chunks(channels)
                .map(|frame| frame.iter().sum::<f32>() / channels as f32)
                .collect();

            let mut s = state.lock().unwrap();

            // RMS volume
            let rms = (mono.iter().map(|x| x * x).sum::<f32>() / mono.len() as f32).sqrt();
            s.volume = rms;

            // Waveform ring buffer
            s.waveform.extend_from_slice(&mono);
            if s.waveform.len() > waveform_len {
                let excess = s.waveform.len() - waveform_len;
                s.waveform.drain(..excess);
            }

            // Accumulate for FFT
            accumulator.extend_from_slice(&mono);
            if accumulator.len() >= fft_size {
                s.spectrum = compute_spectrum(&accumulator[..fft_size], fft_size);
                accumulator.clear();
            }

            // Recording accumulator
            if s.recording {
                s.recorded_samples.extend_from_slice(&mono);
            }
        },
        |err| {
            log::error!("audio stream error: {}", err);
        },
        None,
    )?;

    stream.play()?;
    Ok(stream)
}

// ── WAV writer ──

/// Write mono f32 samples as a 16-bit PCM WAV file.
pub fn write_wav(path: &std::path::Path, samples: &[f32], sample_rate: u32) -> std::io::Result<()> {
    use std::io::Write;

    let num_samples = samples.len() as u32;
    let bytes_per_sample = 2u16; // 16-bit
    let data_size = num_samples * bytes_per_sample as u32;
    let file_size = 36 + data_size;

    let mut f = std::fs::File::create(path)?;

    // RIFF header
    f.write_all(b"RIFF")?;
    f.write_all(&file_size.to_le_bytes())?;
    f.write_all(b"WAVE")?;

    // fmt chunk
    f.write_all(b"fmt ")?;
    f.write_all(&16u32.to_le_bytes())?;      // chunk size
    f.write_all(&1u16.to_le_bytes())?;        // PCM format
    f.write_all(&1u16.to_le_bytes())?;        // mono
    f.write_all(&sample_rate.to_le_bytes())?;  // sample rate
    f.write_all(&(sample_rate * bytes_per_sample as u32).to_le_bytes())?; // byte rate
    f.write_all(&bytes_per_sample.to_le_bytes())?; // block align
    f.write_all(&16u16.to_le_bytes())?;       // bits per sample

    // data chunk
    f.write_all(b"data")?;
    f.write_all(&data_size.to_le_bytes())?;

    for &s in samples {
        let clamped = s.clamp(-1.0, 1.0);
        let i16_val = (clamped * 32767.0) as i16;
        f.write_all(&i16_val.to_le_bytes())?;
    }

    Ok(())
}

/// Compute spectral centroid (brightness) in Hz from a magnitude spectrum.
pub fn spectral_centroid(spectrum: &[f32], sample_rate: u32) -> f32 {
    let bin_hz = sample_rate as f32 / (spectrum.len() as f32 * 2.0);
    let mut weighted_sum = 0.0f32;
    let mut total_mag = 0.0f32;
    for (i, &mag_db) in spectrum.iter().enumerate() {
        let mag_linear = 10.0f32.powf(mag_db / 20.0);
        let freq = i as f32 * bin_hz;
        weighted_sum += freq * mag_linear;
        total_mag += mag_linear;
    }
    if total_mag > 0.0 { weighted_sum / total_mag } else { 0.0 }
}

/// Find the top N peaks in a magnitude spectrum, returned as (frequency_hz, magnitude_db).
pub fn find_peaks(spectrum: &[f32], sample_rate: u32, n: usize) -> Vec<(f32, f32)> {
    let bin_hz = sample_rate as f32 / (spectrum.len() as f32 * 2.0);
    let mut peaks: Vec<(f32, f32)> = Vec::new();

    for i in 1..spectrum.len() - 1 {
        if spectrum[i] > spectrum[i - 1] && spectrum[i] > spectrum[i + 1] && spectrum[i] > -40.0 {
            let freq = i as f32 * bin_hz;
            peaks.push((freq, spectrum[i]));
        }
    }

    peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    peaks.truncate(n);
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    peaks
}

/// Convert frequency to MIDI note number (fractional).
pub fn freq_to_midi(freq: f32) -> f32 {
    69.0 + 12.0 * (freq / 440.0).log2()
}

/// Convert MIDI note number to note name string.
pub fn midi_to_name(midi: u8) -> String {
    const NAMES: [&str; 12] = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"];
    let octave = (midi / 12) as i8 - 1;
    format!("{}{}", NAMES[(midi % 12) as usize], octave)
}
