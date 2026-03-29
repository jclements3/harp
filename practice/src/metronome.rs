use std::sync::{Arc, Mutex};

use cpal::traits::{DeviceTrait, HostTrait, StreamTrait};
use cpal::Stream;

pub struct Metronome {
    _stream: Option<Stream>,
    trigger: Arc<Mutex<bool>>,
}

impl Metronome {
    pub fn new() -> Self {
        let trigger = Arc::new(Mutex::new(false));
        let stream = Self::start_output(Arc::clone(&trigger));
        Self {
            _stream: stream,
            trigger,
        }
    }

    pub fn click(&self) {
        if let Ok(mut t) = self.trigger.lock() {
            *t = true;
        }
    }

    fn start_output(trigger: Arc<Mutex<bool>>) -> Option<Stream> {
        let host = cpal::default_host();
        let device = host.default_output_device()?;
        let config = device.default_output_config().ok()?;
        let sample_rate = config.sample_rate().0 as f32;
        let channels = config.channels() as usize;

        let click_freq = 880.0f32;
        let click_samples = (sample_rate * 0.04) as usize;

        let phase = Arc::new(Mutex::new(0usize));
        let remaining = Arc::new(Mutex::new(0usize));
        let trigger_clone = trigger;

        let stream = device
            .build_output_stream(
                &config.into(),
                move |data: &mut [f32], _: &cpal::OutputCallbackInfo| {
                    if let Ok(mut t) = trigger_clone.lock() {
                        if *t {
                            *t = false;
                            if let Ok(mut r) = remaining.lock() {
                                *r = click_samples;
                            }
                            if let Ok(mut p) = phase.lock() {
                                *p = 0;
                            }
                        }
                    }

                    let mut rem = remaining.lock().unwrap();
                    let mut ph = phase.lock().unwrap();

                    for frame in data.chunks_mut(channels) {
                        let sample = if *rem > 0 {
                            let t = *ph as f32 / sample_rate;
                            let env = (-t * 80.0).exp();
                            let val = (2.0 * std::f32::consts::PI * click_freq * t).sin() * env * 0.3;
                            *ph += 1;
                            *rem -= 1;
                            val
                        } else {
                            0.0
                        };
                        for s in frame.iter_mut() {
                            *s = sample;
                        }
                    }
                },
                |err| eprintln!("metronome audio error: {err}"),
                None,
            )
            .ok()?;

        stream.play().ok()?;
        Some(stream)
    }
}
