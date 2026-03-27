pub mod audio;
pub mod session;

use std::sync::{Arc, Mutex};

use eframe::egui;
use cpal::Stream;
use session::{Session, SessionType};

// ── Colors ──
const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241);
const CARD_BG: egui::Color32 = egui::Color32::WHITE;
const ACCENT: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const RED: egui::Color32 = egui::Color32::from_rgb(220, 38, 38);
const GREEN: egui::Color32 = egui::Color32::from_rgb(40, 167, 69);
const TEXT_MUTED: egui::Color32 = egui::Color32::from_rgb(136, 136, 136);
const TEXT_PRIMARY: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const BORDER: egui::Color32 = egui::Color32::from_rgb(204, 204, 204);
const WAVEFORM_COLOR: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const SPECTRUM_COLOR: egui::Color32 = egui::Color32::from_rgb(139, 92, 246);
const PEAK_COLOR: egui::Color32 = egui::Color32::from_rgb(220, 38, 38);

pub fn create_native_options() -> eframe::NativeOptions {
    eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([800.0, 700.0]),
        ..Default::default()
    }
}

pub fn run_app(options: eframe::NativeOptions) -> eframe::Result {
    eframe::run_native(
        "Harp Sampler",
        options,
        Box::new(|cc| {
            apply_style(&cc.egui_ctx);
            Ok(Box::new(SamplerApp::new()))
        }),
    )
}

fn apply_style(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing = egui::Vec2::new(8.0, 6.0);
    style.spacing.button_padding = egui::Vec2::new(14.0, 8.0);
    style.spacing.window_margin = egui::Margin::same(16);
    let r8 = egui::CornerRadius::same(8);
    style.visuals.widgets.noninteractive.corner_radius = r8;
    style.visuals.widgets.inactive.corner_radius = r8;
    style.visuals.widgets.hovered.corner_radius = r8;
    style.visuals.widgets.active.corner_radius = r8;
    style.visuals.panel_fill = BG;
    style.visuals.window_fill = BG;
    style.visuals.widgets.inactive.bg_stroke = egui::Stroke::new(2.0, BORDER);
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.selection.stroke = egui::Stroke::new(2.0, ACCENT);
    ctx.set_style(style);
}

// ── Screens ──

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Screen {
    Setup,
    Capture,
    Review,
}

// ── App state ──

struct SamplerApp {
    screen: Screen,
    audio_state: Arc<Mutex<audio::AudioState>>,
    _stream: Option<Stream>,
    audio_error: Option<String>,

    // Setup options
    session_type: SessionType,
    output_dir: String,

    // Active session
    session: Option<Session>,

    // Capture state
    is_recording: bool,
    #[allow(dead_code)]
    countdown: Option<(std::time::Instant, u8)>, // (start_time, seconds) — future: auto-record timer
    last_save_msg: Option<String>,
    free_label: String,

    // For continuous repaint during capture
    #[allow(dead_code)]
    needs_repaint: bool, // future: throttle repaints
}

impl SamplerApp {
    fn new() -> Self {
        let default_dir = if cfg!(target_os = "android") {
            "/storage/emulated/0/HarpSampler".to_string()
        } else {
            dirs_or_cwd("harp_samples")
        };

        Self {
            screen: Screen::Setup,
            audio_state: Arc::new(Mutex::new(audio::AudioState::default())),
            _stream: None,
            audio_error: None,
            session_type: SessionType::Chromatic,
            output_dir: default_dir,
            session: None,
            is_recording: false,
            countdown: None,
            last_save_msg: None,
            free_label: String::new(),
            needs_repaint: false,
        }
    }

    fn start_mic(&mut self) {
        if self._stream.is_some() {
            return;
        }
        match audio::start_capture(Arc::clone(&self.audio_state)) {
            Ok(s) => {
                self._stream = Some(s);
                self.audio_error = None;
            }
            Err(e) => {
                self.audio_error = Some(format!("Mic error: {e}"));
            }
        }
    }

    fn stop_mic(&mut self) {
        {
            let mut s = self.audio_state.lock().unwrap();
            s.running = false;
        }
        self._stream = None;
    }

    fn start_recording(&mut self) {
        let mut s = self.audio_state.lock().unwrap();
        s.recorded_samples.clear();
        s.recording = true;
        drop(s);
        self.is_recording = true;
        self.last_save_msg = None;
    }

    fn stop_recording(&mut self) {
        let mut s = self.audio_state.lock().unwrap();
        s.recording = false;
        drop(s);
        self.is_recording = false;
    }

    fn save_current_recording(&mut self) {
        let s = self.audio_state.lock().unwrap();
        let samples = s.recorded_samples.clone();
        let sample_rate = s.sample_rate;
        let rms = s.volume;
        let spectrum = s.spectrum.clone();
        drop(s);

        if samples.is_empty() {
            self.last_save_msg = Some("No audio recorded!".into());
            return;
        }

        if let Some(ref mut session) = self.session {
            let label = if session.session_type == SessionType::Free {
                Some(self.free_label.clone())
            } else {
                None
            };
            match session.save_sample(&samples, sample_rate, rms, &spectrum, label) {
                Ok(meta) => {
                    self.last_save_msg = Some(format!(
                        "Saved {} ({:.1}s, centroid {:.0} Hz)",
                        meta.wav_filename, meta.duration_secs, meta.spectral_centroid_hz
                    ));
                    self.free_label.clear();
                }
                Err(e) => {
                    self.last_save_msg = Some(format!("Error: {e}"));
                }
            }
        }
    }
}

impl eframe::App for SamplerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Continuous repaint during capture screen
        if self.screen == Screen::Capture {
            ctx.request_repaint();
        }

        match self.screen {
            Screen::Setup => self.ui_setup(ctx),
            Screen::Capture => self.ui_capture(ctx),
            Screen::Review => self.ui_review(ctx),
        }
    }
}

// ── Setup screen ──

impl SamplerApp {
    fn ui_setup(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.add_space(20.0);
            ui.heading(egui::RichText::new("Harp Sampler").strong().size(28.0).color(TEXT_PRIMARY));
            ui.add_space(4.0);
            ui.label(egui::RichText::new(
                "Record audio from your lever harp to build a sample library for pitch detection optimization."
            ).color(TEXT_MUTED));
            ui.add_space(20.0);

            // Session type cards
            ui.label(egui::RichText::new("Session Type").strong().size(16.0));
            ui.add_space(4.0);

            let types = [
                (SessionType::Chromatic, "Chromatic Scale",
                 "Every string natural + sharped. ~82 samples for a 34-string harp.\nBest for building a complete fingerprint of your instrument."),
                (SessionType::Chords, "Common Chords",
                 "Major, minor, dom7, dim triads in 8 keys. 32 chord samples.\nCaptures how your harp's overtones interact in harmony."),
                (SessionType::Free, "Free Recording",
                 "Record anything with manual labels. Good for specific passages,\nglissandi, harmonics, or edge cases."),
            ];

            for (st, title, desc) in types {
                let selected = self.session_type == st;
                let frame = egui::Frame::group(ui.style())
                    .fill(if selected { egui::Color32::from_rgb(219, 234, 254) } else { CARD_BG })
                    .stroke(egui::Stroke::new(
                        if selected { 2.0 } else { 1.0 },
                        if selected { ACCENT } else { BORDER },
                    ))
                    .inner_margin(12.0)
                    .corner_radius(8.0);

                let response = frame.show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut self.session_type, st, "");
                        ui.vertical(|ui| {
                            ui.label(egui::RichText::new(title).strong().size(14.0));
                            ui.label(egui::RichText::new(desc).small().color(TEXT_MUTED));
                        });
                    });
                }).response;

                if response.clicked() {
                    self.session_type = st;
                }
                ui.add_space(4.0);
            }

            ui.add_space(12.0);

            // Output directory
            ui.label(egui::RichText::new("Output Directory").strong().size(16.0));
            ui.horizontal(|ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut self.output_dir)
                        .desired_width(400.0)
                        .hint_text("Path to save samples...")
                );
                #[cfg(not(target_os = "android"))]
                if ui.button("Browse...").clicked() {
                    if let Some(dir) = rfd::FileDialog::new().pick_folder() {
                        self.output_dir = dir.to_string_lossy().to_string();
                    }
                }
            });

            ui.add_space(20.0);

            // Start button
            let btn = egui::Button::new(
                egui::RichText::new("  Start Sampling Session  ").size(18.0).color(egui::Color32::WHITE)
            ).fill(ACCENT).corner_radius(8.0);

            if ui.add(btn).clicked() {
                let dir = std::path::PathBuf::from(&self.output_dir);
                self.session = Some(Session::new(self.session_type, dir));
                self.start_mic();
                self.screen = Screen::Capture;
            }

            if let Some(ref err) = self.audio_error {
                ui.colored_label(RED, err);
            }
        });
    }
}

// ── Capture screen ──

impl SamplerApp {
    fn ui_capture(&mut self, ctx: &egui::Context) {
        let state = self.audio_state.lock().unwrap();
        let waveform = state.waveform.clone();
        let spectrum = state.spectrum.clone();
        let volume = state.volume;
        let sample_rate = state.sample_rate;
        let rec_len = state.recorded_samples.len();
        drop(state);

        let rec_secs = rec_len as f32 / sample_rate as f32;

        // Top panel: progress + nav
        egui::TopBottomPanel::top("capture_top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("< Back").clicked() {
                    self.stop_recording();
                    self.stop_mic();
                    self.screen = Screen::Setup;
                }
                ui.separator();

                if let Some(ref session) = self.session {
                    let (done, total) = session.progress();
                    if session.session_type != SessionType::Free {
                        ui.label(egui::RichText::new(
                            format!("Progress: {}/{}", done, total)
                        ).strong());

                        // Progress bar
                        let frac = if total > 0 { done as f32 / total as f32 } else { 0.0 };
                        ui.add(egui::ProgressBar::new(frac).desired_width(200.0));
                    } else {
                        ui.label(egui::RichText::new(
                            format!("{} samples recorded", done)
                        ).strong());
                    }
                }

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("Review Samples >>").clicked() {
                        self.stop_recording();
                        self.screen = Screen::Review;
                    }
                });
            });
        });

        // Bottom: status
        egui::TopBottomPanel::bottom("capture_bottom").show(ctx, |ui| {
            if let Some(ref msg) = self.last_save_msg {
                let color = if msg.starts_with("Error") { RED } else { GREEN };
                ui.colored_label(color, msg);
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let session_complete = self.session.as_ref().map(|s| s.is_complete()).unwrap_or(false);

            if session_complete {
                ui.add_space(40.0);
                ui.heading(egui::RichText::new("Session Complete!").color(GREEN).strong());
                ui.add_space(12.0);
                if let Some(ref session) = self.session {
                    ui.label(format!("{} samples collected.", session.completed.len()));
                    ui.label(format!("Saved to: {}", session.output_dir.display()));
                }
                ui.add_space(16.0);
                if ui.button("Review Samples").clicked() {
                    self.screen = Screen::Review;
                }
                return;
            }

            // Target note display
            if let Some(ref session) = self.session {
                if let Some(target) = session.current_target() {
                    ui.add_space(8.0);
                    egui::Frame::group(ui.style())
                        .fill(CARD_BG)
                        .inner_margin(16.0)
                        .corner_radius(8.0)
                        .show(ui, |ui| {
                            ui.horizontal(|ui| {
                                ui.label(egui::RichText::new("Play:").color(TEXT_MUTED).size(16.0));
                                ui.label(egui::RichText::new(&target.label)
                                    .strong().size(32.0).color(ACCENT));
                            });
                            if target.lever {
                                ui.label(egui::RichText::new("(engage lever before plucking)")
                                    .color(RED).italics());
                            }
                            if !target.midi_notes.is_empty() {
                                let note_names: Vec<String> = target.midi_notes.iter()
                                    .map(|&m| audio::midi_to_name(m))
                                    .collect();
                                ui.label(egui::RichText::new(
                                    format!("Expected: {}", note_names.join(", "))
                                ).color(TEXT_MUTED).small());
                            }
                        });
                }
            }

            // Free mode label input
            if self.session.as_ref().map(|s| s.session_type == SessionType::Free).unwrap_or(false) {
                ui.horizontal(|ui| {
                    ui.label("Label:");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.free_label)
                            .hint_text("e.g., C4, glissando, ...")
                            .desired_width(200.0)
                    );
                });
            }

            ui.add_space(8.0);

            // Volume meter
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new("Level:").color(TEXT_MUTED));
                let vol_bar = volume.min(1.0);
                let bar_color = if vol_bar > 0.8 { RED } else if vol_bar > 0.3 { GREEN } else { TEXT_MUTED };
                ui.add(egui::ProgressBar::new(vol_bar)
                    .desired_width(300.0)
                    .fill(bar_color));
                ui.label(egui::RichText::new(format!("{:.1} dB", 20.0 * volume.max(1e-8).log10()))
                    .monospace().small());
            });

            // Detected peaks
            if !spectrum.is_empty() {
                let peaks = audio::find_peaks(&spectrum, sample_rate, 5);
                if !peaks.is_empty() {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("Peaks:").color(TEXT_MUTED));
                        for (freq, _mag) in &peaks {
                            let midi = audio::freq_to_midi(*freq);
                            let midi_round = midi.round() as u8;
                            let cents = (midi - midi_round as f32) * 100.0;
                            let name = audio::midi_to_name(midi_round);
                            ui.label(egui::RichText::new(
                                format!("{} ({:+.0}c) {:.0}Hz", name, cents, freq)
                            ).monospace().small().color(PEAK_COLOR));
                        }
                    });
                }
            }

            ui.add_space(8.0);

            // Waveform display
            let waveform_height = 80.0;
            let (wf_rect, _) = ui.allocate_exact_size(
                egui::Vec2::new(ui.available_width(), waveform_height),
                egui::Sense::hover(),
            );
            draw_waveform(ui, wf_rect, &waveform, self.is_recording);

            ui.add_space(4.0);

            // Spectrum display
            let spectrum_height = 120.0;
            let (sp_rect, _) = ui.allocate_exact_size(
                egui::Vec2::new(ui.available_width(), spectrum_height),
                egui::Sense::hover(),
            );
            draw_spectrum(ui, sp_rect, &spectrum, sample_rate);

            ui.add_space(12.0);

            // Record controls
            ui.horizontal(|ui| {
                if self.is_recording {
                    // Recording indicator
                    ui.label(egui::RichText::new("  REC").color(RED).strong().size(18.0));
                    ui.label(egui::RichText::new(format!("{:.1}s", rec_secs)).monospace().size(16.0));

                    let stop_btn = egui::Button::new(
                        egui::RichText::new("  Stop & Save  ").size(16.0).color(egui::Color32::WHITE)
                    ).fill(GREEN).corner_radius(8.0);
                    if ui.add(stop_btn).clicked() {
                        self.stop_recording();
                        self.save_current_recording();
                    }

                    let discard_btn = egui::Button::new(
                        egui::RichText::new("  Discard  ").size(14.0)
                    ).corner_radius(8.0);
                    if ui.add(discard_btn).clicked() {
                        self.stop_recording();
                        self.last_save_msg = Some("Discarded.".into());
                    }
                } else {
                    let rec_btn = egui::Button::new(
                        egui::RichText::new("  Record  ").size(18.0).color(egui::Color32::WHITE)
                    ).fill(RED).corner_radius(8.0);
                    if ui.add(rec_btn).clicked() {
                        self.start_recording();
                    }

                    if self.session.as_ref().map(|s| s.session_type != SessionType::Free).unwrap_or(false) {
                        if ui.button("Skip >>").clicked() {
                            if let Some(ref mut session) = self.session {
                                session.skip();
                            }
                        }
                    }
                }
            });
        });
    }
}

// ── Review screen ──

impl SamplerApp {
    fn ui_review(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("review_top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("< Continue Recording").clicked() {
                    self.screen = Screen::Capture;
                }
                ui.separator();
                ui.heading(egui::RichText::new("Sample Review").strong());
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            if let Some(ref session) = self.session {
                ui.label(egui::RichText::new(
                    format!("Session: {} | {} samples", session.session_name, session.completed.len())
                ).color(TEXT_MUTED));
                ui.label(egui::RichText::new(
                    format!("Output: {}", session.output_dir.display())
                ).color(TEXT_MUTED).small());
                ui.add_space(8.0);

                if session.completed.is_empty() {
                    ui.label("No samples recorded yet.");
                    return;
                }

                // Summary stats
                egui::Frame::group(ui.style())
                    .fill(CARD_BG)
                    .inner_margin(12.0)
                    .corner_radius(8.0)
                    .show(ui, |ui| {
                        let total_dur: f32 = session.completed.iter().map(|s| s.duration_secs).sum();
                        let avg_centroid: f32 = session.completed.iter()
                            .map(|s| s.spectral_centroid_hz).sum::<f32>() / session.completed.len() as f32;
                        ui.horizontal(|ui| {
                            ui.label(egui::RichText::new(
                                format!("{} samples", session.completed.len())
                            ).strong());
                            ui.separator();
                            ui.label(format!("Total: {:.1}s", total_dur));
                            ui.separator();
                            ui.label(format!("Avg centroid: {:.0} Hz", avg_centroid));
                        });
                    });

                ui.add_space(8.0);

                // Sample table
                egui::ScrollArea::vertical().show(ui, |ui| {
                    egui::Grid::new("sample_table")
                        .striped(true)
                        .min_col_width(40.0)
                        .show(ui, |ui| {
                            ui.label(egui::RichText::new("#").strong());
                            ui.label(egui::RichText::new("Label").strong());
                            ui.label(egui::RichText::new("Duration").strong());
                            ui.label(egui::RichText::new("RMS").strong());
                            ui.label(egui::RichText::new("Centroid").strong());
                            ui.label(egui::RichText::new("Peaks").strong());
                            ui.label(egui::RichText::new("File").strong());
                            ui.end_row();

                            for (i, meta) in session.completed.iter().enumerate() {
                                ui.label(format!("{}", i + 1));
                                ui.label(egui::RichText::new(&meta.target.label).color(ACCENT));
                                ui.label(format!("{:.1}s", meta.duration_secs));
                                ui.label(format!("{:.3}", meta.rms_volume));
                                ui.label(format!("{:.0} Hz", meta.spectral_centroid_hz));

                                // Peak frequencies as note names
                                let peak_names: Vec<String> = meta.peak_frequencies_hz.iter()
                                    .take(3)
                                    .map(|&f| {
                                        let m = audio::freq_to_midi(f).round() as u8;
                                        audio::midi_to_name(m)
                                    })
                                    .collect();
                                ui.label(egui::RichText::new(peak_names.join(", ")).small());

                                ui.label(egui::RichText::new(&meta.wav_filename).small().color(TEXT_MUTED));
                                ui.end_row();
                            }
                        });
                });
            } else {
                ui.label("No active session.");
            }
        });
    }
}

// ── Drawing helpers ──

fn draw_waveform(ui: &mut egui::Ui, rect: egui::Rect, waveform: &[f32], recording: bool) {
    let painter = ui.painter_at(rect);

    // Background
    painter.rect_filled(rect, 4.0, CARD_BG);
    painter.rect_stroke(rect, 4.0, egui::Stroke::new(1.0, BORDER), egui::StrokeKind::Outside);

    if recording {
        painter.rect_stroke(rect, 4.0, egui::Stroke::new(2.0, RED), egui::StrokeKind::Outside);
    }

    if waveform.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "waiting for audio...",
            egui::FontId::monospace(12.0),
            TEXT_MUTED,
        );
        return;
    }

    // Center line
    let cy = rect.center().y;
    painter.line_segment(
        [egui::Pos2::new(rect.left(), cy), egui::Pos2::new(rect.right(), cy)],
        egui::Stroke::new(0.5, BORDER),
    );

    // Waveform
    let n = waveform.len();
    let w = rect.width();
    let h = rect.height() * 0.45;
    let step = (n as f32 / w).max(1.0) as usize;

    let points: Vec<egui::Pos2> = (0..w as usize)
        .map(|px| {
            let idx = (px as f32 * n as f32 / w) as usize;
            let end = (idx + step).min(n);
            let max_val = waveform[idx..end]
                .iter()
                .map(|s| s.abs())
                .fold(0.0f32, f32::max);
            egui::Pos2::new(
                rect.left() + px as f32,
                cy - max_val * h,
            )
        })
        .collect();

    let points_neg: Vec<egui::Pos2> = points.iter()
        .map(|p| egui::Pos2::new(p.x, cy + (cy - p.y)))
        .collect();

    // Draw as filled area
    let color = if recording { RED } else { WAVEFORM_COLOR };
    for i in 0..points.len() {
        painter.line_segment([points[i], points_neg[i]], egui::Stroke::new(1.0, color));
    }

    // Label
    painter.text(
        egui::Pos2::new(rect.left() + 4.0, rect.top() + 2.0),
        egui::Align2::LEFT_TOP,
        "Waveform",
        egui::FontId::monospace(10.0),
        TEXT_MUTED,
    );
}

fn draw_spectrum(ui: &mut egui::Ui, rect: egui::Rect, spectrum: &[f32], sample_rate: u32) {
    let painter = ui.painter_at(rect);

    painter.rect_filled(rect, 4.0, CARD_BG);
    painter.rect_stroke(rect, 4.0, egui::Stroke::new(1.0, BORDER), egui::StrokeKind::Outside);

    if spectrum.is_empty() || spectrum.iter().all(|&v| v <= -80.0) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "spectrum",
            egui::FontId::monospace(12.0),
            TEXT_MUTED,
        );
        return;
    }

    // We display 0 to ~4kHz (relevant range for harp fundamentals + first harmonics)
    let nyquist = sample_rate as f32 / 2.0;
    let max_display_hz = 4000.0f32.min(nyquist);
    let max_bin = ((max_display_hz / nyquist) * spectrum.len() as f32) as usize;
    let display_bins = &spectrum[1..max_bin.min(spectrum.len())]; // skip DC

    let db_min = -60.0f32;
    let db_max = 0.0f32;

    let w = rect.width();
    let h = rect.height();
    let bar_w = w / display_bins.len() as f32;

    for (i, &db) in display_bins.iter().enumerate() {
        let norm = ((db - db_min) / (db_max - db_min)).clamp(0.0, 1.0);
        let x = rect.left() + i as f32 * bar_w;
        let bar_h = norm * (h - 16.0);
        let y = rect.bottom() - bar_h;

        let color = if norm > 0.7 { PEAK_COLOR } else { SPECTRUM_COLOR };
        painter.rect_filled(
            egui::Rect::from_min_size(egui::Pos2::new(x, y), egui::Vec2::new(bar_w.max(1.0), bar_h)),
            0.0,
            color.linear_multiply(0.7),
        );
    }

    // Frequency axis labels
    for khz in [0.5, 1.0, 2.0, 3.0, 4.0] {
        let freq = khz * 1000.0;
        if freq > max_display_hz { break; }
        let x_frac = freq / max_display_hz;
        let x = rect.left() + x_frac * w;
        painter.text(
            egui::Pos2::new(x, rect.bottom() - 2.0),
            egui::Align2::CENTER_BOTTOM,
            format!("{}k", khz),
            egui::FontId::monospace(9.0),
            TEXT_MUTED,
        );
    }

    painter.text(
        egui::Pos2::new(rect.left() + 4.0, rect.top() + 2.0),
        egui::Align2::LEFT_TOP,
        "Spectrum (0-4kHz)",
        egui::FontId::monospace(10.0),
        TEXT_MUTED,
    );
}

fn dirs_or_cwd(subdir: &str) -> String {
    if let Ok(home) = std::env::var("HOME") {
        format!("{}/{}", home, subdir)
    } else {
        subdir.to_string()
    }
}

// ── Android entry point ──
#[cfg(target_os = "android")]
#[unsafe(no_mangle)]
fn android_main(app: winit::platform::android::activity::AndroidApp) {
    android_logger::init_once(
        android_logger::Config::default().with_max_level(log::LevelFilter::Info),
    );

    let options = eframe::NativeOptions {
        android_app: Some(app),
        renderer: eframe::Renderer::Glow,
        ..Default::default()
    };

    run_app(options).unwrap();
}
