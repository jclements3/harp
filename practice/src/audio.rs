use std::sync::{Arc, Mutex};

use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::Stream;
use pitch_detection::detector::mcleod::McLeodDetector;
use pitch_detection::detector::PitchDetector;

use crate::music::{self, Note};

const DETECTION_BUFFER_SIZE: usize = 2048;
const POWER_THRESHOLD: f64 = 5.0;
const CLARITY_THRESHOLD: f64 = 0.6;

pub struct AudioState {
    pub current_note: Option<Note>,
    pub note_history: Vec<Note>,
    pub volume: f32,
    pub running: bool,
    pub sample_rate: u32,
}

impl Default for AudioState {
    fn default() -> Self {
        Self {
            current_note: None,
            note_history: Vec::new(),
            volume: 0.0,
            running: false,
            sample_rate: 44100,
        }
    }
}

pub fn start_capture(
    state: Arc<Mutex<AudioState>>,
) -> Result<Stream, Box<dyn std::error::Error>> {
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

    let mut sample_buf: Vec<f64> = Vec::with_capacity(DETECTION_BUFFER_SIZE);

    let state_clone = Arc::clone(&state);
    let stream = device.build_input_stream(
        &config.into(),
        move |data: &[f32], _: &cpal::InputCallbackInfo| {
            let rms = (data.iter().map(|&s| s * s).sum::<f32>() / data.len() as f32).sqrt();

            for &s in data {
                sample_buf.push(s as f64);
            }

            if sample_buf.len() >= DETECTION_BUFFER_SIZE {
                let mut detector = McLeodDetector::new(DETECTION_BUFFER_SIZE, DETECTION_BUFFER_SIZE / 2);
                let pitch = detector.get_pitch(
                    &sample_buf[..DETECTION_BUFFER_SIZE],
                    sample_rate as usize,
                    POWER_THRESHOLD,
                    CLARITY_THRESHOLD,
                );

                let note = pitch.and_then(|p| music::freq_to_note(p.frequency as f32));

                if let Ok(mut s) = state_clone.lock() {
                    s.volume = rms;
                    s.current_note = note.clone();
                    if let Some(n) = note {
                        s.note_history.push(n);
                        let max_len = (sample_rate as usize / DETECTION_BUFFER_SIZE) * 2;
                        if s.note_history.len() > max_len {
                            s.note_history.remove(0);
                        }
                    }
                }

                sample_buf.clear();
            }
        },
        |err| eprintln!("audio error: {err}"),
        None,
    )?;

    stream.play()?;
    Ok(stream)
}
