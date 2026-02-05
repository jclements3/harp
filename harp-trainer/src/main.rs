mod audio;
mod music;

use std::sync::{Arc, Mutex};

use eframe::egui;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([480.0, 640.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Harp Trainer",
        options,
        Box::new(|_cc| Ok(Box::new(HarpApp::new()))),
    )
}

struct HarpApp {
    state: Arc<Mutex<audio::AudioState>>,
    _stream: Option<cpal::Stream>,
    error: Option<String>,
}

impl HarpApp {
    fn new() -> Self {
        let state = Arc::new(Mutex::new(audio::AudioState::default()));
        let (stream, error) = match audio::start_capture(Arc::clone(&state)) {
            Ok(s) => (Some(s), None),
            Err(e) => (None, Some(format!("Audio error: {e}"))),
        };
        Self {
            state,
            _stream: stream,
            error,
        }
    }
}

impl eframe::App for HarpApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Request continuous repaint for real-time updates
        ctx.request_repaint();

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Harp Trainer");
            ui.separator();

            if let Some(err) = &self.error {
                ui.colored_label(egui::Color32::RED, err);
                return;
            }

            let state = self.state.lock().unwrap();

            // -- Volume meter --
            ui.horizontal(|ui| {
                ui.label("Volume:");
                let bar = egui::ProgressBar::new(state.volume.clamp(0.0, 1.0));
                ui.add(bar);
            });

            ui.add_space(12.0);

            // -- Current detected note --
            match &state.current_note {
                Some(note) => {
                    ui.horizontal(|ui| {
                        ui.label(
                            egui::RichText::new(format!("{}{}", note.name, note.octave))
                                .size(64.0)
                                .strong(),
                        );

                        // Cents indicator
                        let cents_text = if note.cents_off >= 0.0 {
                            format!("+{:.0}\u{00a2}", note.cents_off)
                        } else {
                            format!("{:.0}\u{00a2}", note.cents_off)
                        };
                        let cents_color = if note.cents_off.abs() < 10.0 {
                            egui::Color32::from_rgb(40, 167, 69)
                        } else if note.cents_off.abs() < 25.0 {
                            egui::Color32::from_rgb(255, 193, 7)
                        } else {
                            egui::Color32::from_rgb(220, 53, 69)
                        };
                        ui.label(
                            egui::RichText::new(cents_text).size(24.0).color(cents_color),
                        );
                    });
                }
                None => {
                    ui.label(
                        egui::RichText::new("Play a note...")
                            .size(32.0)
                            .color(egui::Color32::GRAY),
                    );
                }
            }

            ui.add_space(12.0);
            ui.separator();

            // -- Recent note history --
            ui.label(egui::RichText::new("Recent Notes").size(18.0));
            ui.add_space(4.0);

            let history_text: String = state
                .note_history
                .iter()
                .rev()
                .take(20)
                .map(|n| format!("{}{}", n.name, n.octave))
                .collect::<Vec<_>>()
                .join("  ");

            ui.label(egui::RichText::new(history_text).size(14.0));

            ui.add_space(12.0);
            ui.separator();

            // -- Chord detection from recent notes --
            ui.label(egui::RichText::new("Detected Chords").size(18.0));
            ui.add_space(4.0);

            // Gather unique pitch classes from last few detections
            let recent: Vec<u8> = state
                .note_history
                .iter()
                .rev()
                .take(8)
                .map(|n| n.midi % 12)
                .collect::<std::collections::HashSet<u8>>()
                .into_iter()
                .collect();

            if recent.len() >= 2 {
                let chords = music::match_chords(&recent);
                if chords.is_empty() {
                    ui.label("(no chord match)");
                } else {
                    for (root, name) in &chords {
                        ui.label(
                            egui::RichText::new(format!("{root} {name}")).size(16.0),
                        );
                    }
                }
            } else {
                ui.label("(need 2+ distinct notes)");
            }
        });
    }
}
