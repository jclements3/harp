/// Simple synthesizer for reference tone playback.
/// Triangle waves strummed bottom to top.

use std::sync::{Arc, Mutex};
use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::Stream;

const MAX_VOICES: usize = 16;
const DECAY_RATE: f32 = 0.6;
const NOTE_VOLUME: f32 = 0.15;

struct Voice {
    freq: f32,
    phase: f32,
    amplitude: f32,
    delay_samples: u32,
    samples_played: u32,
}

pub struct SynthState {
    voices: Vec<Voice>,
    sample_rate: f32,
    pub master_volume: f32,
}

impl SynthState {
    fn new(sample_rate: f32) -> Self {
        Self {
            voices: Vec::new(),
            sample_rate,
            master_volume: 0.5,
        }
    }

    pub fn play_chord(&mut self, midi_notes: &[u8]) {
        self.voices.clear();
        let mut sorted: Vec<u8> = midi_notes.to_vec();
        sorted.sort();
        let strum_gap = 0.03;
        for (i, &midi) in sorted.iter().take(MAX_VOICES).enumerate() {
            let freq = 440.0 * 2.0f32.powf((midi as f32 - 69.0) / 12.0);
            let delay = (i as f32 * strum_gap * self.sample_rate) as u32;
            self.voices.push(Voice {
                freq,
                phase: i as f32 * 0.17,
                amplitude: NOTE_VOLUME,
                delay_samples: delay,
                samples_played: 0,
            });
        }
    }

    pub fn silence(&mut self) {
        self.voices.clear();
    }

    fn next_sample(&mut self) -> f32 {
        let dt = 1.0 / self.sample_rate;
        let mut out = 0.0f32;
        for voice in &mut self.voices {
            voice.samples_played += 1;
            if voice.samples_played < voice.delay_samples {
                continue;
            }
            let t = voice.phase.fract();
            let tri = if t < 0.25 { t * 4.0 }
                      else if t < 0.75 { 2.0 - t * 4.0 }
                      else { t * 4.0 - 4.0 };
            out += voice.amplitude * tri;
            voice.phase += voice.freq * dt;
            voice.amplitude *= 1.0 - DECAY_RATE * dt;
        }
        self.voices.retain(|v| v.amplitude > 0.001);
        (out * self.master_volume).clamp(-1.0, 1.0)
    }
}

pub struct SynthPlayer {
    pub state: Arc<Mutex<SynthState>>,
    _stream: Stream,
}

impl SynthPlayer {
    pub fn new() -> Result<Self, Box<dyn std::error::Error>> {
        let host = cpal::default_host();
        let device = host.default_output_device()
            .ok_or("no output device")?;
        let config = device.default_output_config()?;
        let sample_rate = config.sample_rate().0 as f32;
        let channels = config.channels() as usize;

        let state = Arc::new(Mutex::new(SynthState::new(sample_rate)));
        let state_clone = Arc::clone(&state);

        let stream = device.build_output_stream(
            &config.into(),
            move |data: &mut [f32], _: &cpal::OutputCallbackInfo| {
                let mut synth = state_clone.lock().unwrap();
                for frame in data.chunks_mut(channels) {
                    let sample = synth.next_sample();
                    for s in frame.iter_mut() {
                        *s = sample;
                    }
                }
            },
            |err| log::error!("audio output error: {}", err),
            None,
        )?;

        stream.play()?;

        Ok(SynthPlayer {
            state,
            _stream: stream,
        })
    }

    pub fn play_chord(&self, midi_notes: &[u8]) {
        if let Ok(mut synth) = self.state.lock() {
            synth.play_chord(midi_notes);
        }
    }

    #[allow(dead_code)]
    pub fn silence(&self) {
        if let Ok(mut synth) = self.state.lock() {
            synth.silence();
        }
    }
}
