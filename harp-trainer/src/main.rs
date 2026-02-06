mod audio;
mod drill;
mod metronome;
mod music;
mod staff;

use std::sync::{Arc, Mutex};

use eframe::egui;

use drill::{Difficulty, DrillConfig, DrillSession, InputMode, NoteResult, DURATION_OPTIONS};
use metronome::Metronome;
use music::PedalState;

// ── Colors matching the HTML version ───────────────────────────────

const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241); // #f8f6f1
const CARD_BG: egui::Color32 = egui::Color32::WHITE;
const TEXT_PRIMARY: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const TEXT_SECONDARY: egui::Color32 = egui::Color32::from_rgb(74, 74, 74);
const TEXT_MUTED: egui::Color32 = egui::Color32::from_rgb(136, 136, 136);
const ACCENT: egui::Color32 = egui::Color32::from_rgb(37, 99, 235); // #2563eb
const RED: egui::Color32 = egui::Color32::from_rgb(220, 38, 38);
const GREEN: egui::Color32 = egui::Color32::from_rgb(40, 167, 69);
const YELLOW: egui::Color32 = egui::Color32::from_rgb(255, 193, 7);
const BORDER: egui::Color32 = egui::Color32::from_rgb(204, 204, 204);
const TRACK_BG: egui::Color32 = egui::Color32::from_rgb(232, 232, 232);

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([480.0, 740.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Harp Trainer",
        options,
        Box::new(|cc| {
            apply_style(&cc.egui_ctx);
            Ok(Box::new(HarpApp::new()))
        }),
    )
}

fn apply_style(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();

    // Spacing
    style.spacing.item_spacing = egui::Vec2::new(8.0, 6.0);
    style.spacing.button_padding = egui::Vec2::new(12.0, 6.0);
    style.spacing.window_margin = egui::Margin::same(16);

    // Rounding
    let r8 = egui::CornerRadius::same(8);
    style.visuals.widgets.noninteractive.corner_radius = r8;
    style.visuals.widgets.inactive.corner_radius = r8;
    style.visuals.widgets.hovered.corner_radius = r8;
    style.visuals.widgets.active.corner_radius = r8;

    // Colors — warm theme
    style.visuals.panel_fill = BG;
    style.visuals.window_fill = BG;
    style.visuals.extreme_bg_color = CARD_BG;

    // Widget colors
    style.visuals.widgets.noninteractive.bg_fill = CARD_BG;
    style.visuals.widgets.noninteractive.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);

    style.visuals.widgets.inactive.bg_fill = CARD_BG;
    style.visuals.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, TEXT_SECONDARY);
    style.visuals.widgets.inactive.weak_bg_fill = CARD_BG;
    style.visuals.widgets.inactive.bg_stroke = egui::Stroke::new(2.0, BORDER);

    style.visuals.widgets.hovered.bg_fill = egui::Color32::from_rgb(240, 240, 240);
    style.visuals.widgets.hovered.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);
    style.visuals.widgets.hovered.bg_stroke = egui::Stroke::new(2.0, egui::Color32::from_rgb(153, 153, 153));

    style.visuals.widgets.active.bg_fill = egui::Color32::from_rgb(219, 234, 254); // #dbeafe
    style.visuals.widgets.active.fg_stroke = egui::Stroke::new(1.0, ACCENT);
    style.visuals.widgets.active.bg_stroke = egui::Stroke::new(2.0, ACCENT);

    style.visuals.selection.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.selection.stroke = egui::Stroke::new(2.0, ACCENT);

    // Separator
    style.visuals.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(238, 238, 238));

    ctx.set_style(style);
}

// ── App state ──────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Screen {
    Setup,
    Drill,
    Results,
}

struct HarpApp {
    screen: Screen,
    audio_state: Arc<Mutex<audio::AudioState>>,
    _stream: Option<cpal::Stream>,
    audio_error: Option<String>,
    drill: Option<DrillSession>,
    config: DrillConfig,
    pedals: PedalState,
    metronome: Metronome,
    confirm_count: u8,
    last_detected: Option<u8>,
    range_preset: usize,
}

const RANGE_PRESETS: [(&str, u8, u8); 3] = [
    ("Bass (C2-B3)", 36, 59),
    ("Treble (C4-C6)", 60, 84),
    ("Full (C2-C6)", 36, 84),
];

impl HarpApp {
    fn new() -> Self {
        let audio_state = Arc::new(Mutex::new(audio::AudioState::default()));
        Self {
            screen: Screen::Setup,
            audio_state,
            _stream: None,
            audio_error: None,
            drill: None,
            config: DrillConfig::default(),
            pedals: PedalState::default(),
            metronome: Metronome::new(),
            confirm_count: 0,
            last_detected: None,
            range_preset: 2,
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
        self._stream = None;
    }

    fn start_drill(&mut self) {
        let (_, low, high) = RANGE_PRESETS[self.range_preset];
        let config = DrillConfig {
            duration_secs: self.config.duration_secs,
            fingers: self.config.fingers,
            low_midi: low,
            high_midi: high,
            difficulty: self.config.difficulty,
            input_mode: self.config.input_mode,
        };
        self.drill = Some(DrillSession::new(config, &self.pedals));
        self.confirm_count = 0;
        self.last_detected = None;
        if self.config.input_mode == InputMode::Mic {
            self.start_mic();
        }
        self.screen = Screen::Drill;
    }
}

// ── Custom widget helpers ──────────────────────────────────────────

/// Draw a pill-shaped toggle button, returns true if clicked.
fn toggle_btn(ui: &mut egui::Ui, label: &str, active: bool) -> bool {
    let text = egui::RichText::new(label).size(13.0);
    let text = if active { text.color(egui::Color32::WHITE) } else { text.color(TEXT_SECONDARY) };

    let btn = egui::Button::new(text)
        .fill(if active { ACCENT } else { CARD_BG })
        .stroke(if active {
            egui::Stroke::new(2.0, ACCENT)
        } else {
            egui::Stroke::new(2.0, BORDER)
        })
        .corner_radius(8.0)
        .min_size(egui::Vec2::new(0.0, 32.0));

    ui.add(btn).clicked()
}

/// Draw a big blue action button.
fn action_btn(ui: &mut egui::Ui, label: &str, size: f32) -> bool {
    let btn = egui::Button::new(
        egui::RichText::new(label).size(size).color(egui::Color32::WHITE),
    )
    .fill(ACCENT)
    .stroke(egui::Stroke::NONE)
    .corner_radius(10.0)
    .min_size(egui::Vec2::new(200.0, 48.0));

    ui.add(btn).clicked()
}

/// Draw a secondary (outline) button.
fn secondary_btn(ui: &mut egui::Ui, label: &str) -> bool {
    let btn = egui::Button::new(
        egui::RichText::new(label).size(15.0).color(TEXT_SECONDARY),
    )
    .fill(CARD_BG)
    .stroke(egui::Stroke::new(2.0, BORDER))
    .corner_radius(10.0)
    .min_size(egui::Vec2::new(180.0, 40.0));

    ui.add(btn).clicked()
}

/// Section heading label.
fn section_label(ui: &mut egui::Ui, text: &str) {
    ui.label(
        egui::RichText::new(text)
            .size(14.0)
            .strong()
            .color(TEXT_SECONDARY),
    );
}

// ── UI ─────────────────────────────────────────────────────────────

impl eframe::App for HarpApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.request_repaint();
        match self.screen {
            Screen::Setup => self.draw_setup(ctx),
            Screen::Drill => self.draw_drill(ctx),
            Screen::Results => self.draw_results(ctx),
        }
    }
}

impl HarpApp {
    // ── Setup screen ───────────────────────────────────────────────

    fn draw_setup(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add_space(8.0);
                    ui.label(
                        egui::RichText::new("HARP SIGHT-READING")
                            .size(14.0)
                            .strong()
                            .color(TEXT_SECONDARY),
                    );
                });

                ui.add_space(8.0);

                // ── Key signature card ──
                self.draw_card(ui, |ui, pedals| {
                    // Key name + flat/sharp buttons
                    ui.vertical_centered(|ui| {
                        ui.horizontal(|ui| {
                            let flat_btn = egui::Button::new(
                                egui::RichText::new("\u{266d}+").size(18.0).strong().color(ACCENT),
                            )
                            .fill(CARD_BG)
                            .stroke(egui::Stroke::new(2.0, egui::Color32::from_rgb(147, 180, 240)))
                            .corner_radius(24.0)
                            .min_size(egui::Vec2::new(44.0, 44.0));
                            if ui.add(flat_btn).clicked() {
                                pedals.add_flat();
                            }

                            ui.add_space(12.0);
                            let key_name = pedals.key_name();
                            ui.vertical(|ui| {
                                ui.label(egui::RichText::new(&key_name).size(16.0).strong().color(TEXT_PRIMARY));
                                let (sharps, flats, is_standard) = pedals.key_sig_info();
                                let detail = if !is_standard {
                                    "custom".to_string()
                                } else if sharps > 0 {
                                    format!("{sharps}\u{266f}")
                                } else if flats > 0 {
                                    format!("{flats}\u{266d}")
                                } else {
                                    "no sharps/flats".to_string()
                                };
                                ui.label(egui::RichText::new(detail).size(11.0).color(TEXT_MUTED));
                            });
                            ui.add_space(12.0);

                            let sharp_btn = egui::Button::new(
                                egui::RichText::new("\u{266f}+").size(18.0).strong().color(RED),
                            )
                            .fill(CARD_BG)
                            .stroke(egui::Stroke::new(2.0, egui::Color32::from_rgb(240, 160, 160)))
                            .corner_radius(24.0)
                            .min_size(egui::Vec2::new(44.0, 44.0));
                            if ui.add(sharp_btn).clicked() {
                                pedals.add_sharp();
                            }
                        });
                    });

                    ui.add_space(8.0);

                    // Pedal diagram with tracks
                    Self::draw_pedal_tracks(ui, pedals);
                });

                ui.add_space(10.0);

                // ── Options card ──
                self.draw_options_card(ui);

                ui.add_space(12.0);

                // ── Start button ──
                ui.vertical_centered(|ui| {
                    if action_btn(ui, "Start Drill", 18.0) {
                        self.start_drill();
                    }
                });

                if let Some(err) = &self.audio_error {
                    ui.add_space(8.0);
                    ui.colored_label(RED, err);
                }

                ui.add_space(16.0);
            });
        });
    }

    /// Draw a white rounded card with shadow.
    fn draw_card(&mut self, ui: &mut egui::Ui, content: impl FnOnce(&mut egui::Ui, &mut PedalState)) {
        let frame = egui::Frame::NONE
            .fill(CARD_BG)
            .corner_radius(12.0)
            .inner_margin(egui::Margin::same(14))
            .shadow(egui::epaint::Shadow {
                offset: [0, 2],
                blur: 12,
                spread: 0,
                color: egui::Color32::from_black_alpha(20),
            });
        frame.show(ui, |ui| {
            content(ui, &mut self.pedals);
        });
    }

    fn draw_options_card(&mut self, ui: &mut egui::Ui) {
        let frame = egui::Frame::NONE
            .fill(CARD_BG)
            .corner_radius(12.0)
            .inner_margin(egui::Margin::same(14))
            .shadow(egui::epaint::Shadow {
                offset: [0, 2],
                blur: 12,
                spread: 0,
                color: egui::Color32::from_black_alpha(20),
            });
        frame.show(ui, |ui| {
            // Session length
            section_label(ui, "SESSION LENGTH");
            ui.horizontal(|ui| {
                for &(label, secs) in &DURATION_OPTIONS {
                    if toggle_btn(ui, label, self.config.duration_secs == secs) {
                        self.config.duration_secs = secs;
                    }
                }
            });

            ui.add_space(8.0);

            // Fingers per position
            section_label(ui, "FINGERS PER POSITION");
            ui.horizontal(|ui| {
                for &f in &[3usize, 4] {
                    let label = format!("{f} fingers");
                    if toggle_btn(ui, &label, self.config.fingers == f) {
                        self.config.fingers = f;
                    }
                }
            });
            ui.label(egui::RichText::new("Consecutive strings within hand span").size(11.0).color(TEXT_MUTED));

            ui.add_space(8.0);

            // Note range
            section_label(ui, "NOTE RANGE");
            ui.horizontal(|ui| {
                for (i, (label, _, _)) in RANGE_PRESETS.iter().enumerate() {
                    if toggle_btn(ui, label, self.range_preset == i) {
                        self.range_preset = i;
                    }
                }
            });
            ui.label(egui::RichText::new("Starts narrow, expands as you improve").size(11.0).color(TEXT_MUTED));

            ui.add_space(8.0);

            // Starting speed
            section_label(ui, "STARTING SPEED");
            ui.horizontal(|ui| {
                for diff in &[Difficulty::Relaxed, Difficulty::Moderate, Difficulty::Brisk] {
                    let label = format!("{} ({:.1}s)", diff.label(), diff.start_secs());
                    if toggle_btn(ui, &label, self.config.difficulty == *diff) {
                        self.config.difficulty = *diff;
                    }
                }
            });
            ui.label(egui::RichText::new("Adapts as you play. 4 on-pace \u{2192} speed up (90%+ accuracy).").size(11.0).color(TEXT_MUTED));

            ui.add_space(8.0);

            // Input mode
            section_label(ui, "INPUT MODE");
            ui.horizontal(|ui| {
                if toggle_btn(ui, "Manual (space)", self.config.input_mode == InputMode::Manual) {
                    self.config.input_mode = InputMode::Manual;
                }
                if toggle_btn(ui, "Microphone", self.config.input_mode == InputMode::Mic) {
                    self.config.input_mode = InputMode::Mic;
                }
            });
            let hint = if self.config.input_mode == InputMode::Manual {
                "Press Space to advance. Speed-only scoring."
            } else {
                "Play the note to advance. Accuracy + speed."
            };
            ui.label(egui::RichText::new(hint).size(11.0).color(TEXT_MUTED));
        });
    }

    /// Draw pedal tracks matching the HTML pedal diagram.
    fn draw_pedal_tracks(ui: &mut egui::Ui, pedals: &mut PedalState) {
        let pedal_names = ["D", "C", "B", "E", "F", "G", "A"];
        let track_w = 24.0f32;
        let track_h = 48.0f32;
        let dot_r = 8.0f32;
        let spacing = 10.0f32;

        // Calculate total width
        let divider_w = 12.0;
        let pedal_w = track_w + spacing;
        let total_w = 3.0 * pedal_w + divider_w + 4.0 * pedal_w - spacing;

        let (rect, response) = ui.allocate_exact_size(
            egui::Vec2::new(total_w, track_h + 24.0),
            egui::Sense::click(),
        );

        let painter = ui.painter_at(rect);
        let start_x = rect.left();
        let label_y = rect.top() + 4.0;
        let track_top = rect.top() + 20.0;

        let mut x = start_x;
        for (i, name) in pedal_names.iter().enumerate() {
            if i == 3 {
                // Divider between left (D,C,B) and right (E,F,G,A) groups
                let div_x = x + divider_w / 2.0;
                painter.line_segment(
                    [
                        egui::Pos2::new(div_x, track_top),
                        egui::Pos2::new(div_x, track_top + track_h),
                    ],
                    egui::Stroke::new(2.0, BORDER),
                );
                x += divider_w;
            }

            let cx = x + track_w / 2.0;
            let pos = pedals.positions[i];

            let (name_color, dot_color) = match pos {
                music::PedalPos::Flat => (ACCENT, ACCENT),
                music::PedalPos::Natural => (egui::Color32::from_rgb(102, 102, 102), egui::Color32::from_rgb(102, 102, 102)),
                music::PedalPos::Sharp => (RED, RED),
            };

            // Pedal letter
            painter.text(
                egui::Pos2::new(cx, label_y),
                egui::Align2::CENTER_TOP,
                name,
                egui::FontId::monospace(12.0),
                name_color,
            );

            // Track background
            let track_rect = egui::Rect::from_min_size(
                egui::Pos2::new(x, track_top),
                egui::Vec2::new(track_w, track_h),
            );
            painter.rect(track_rect, 12.0, TRACK_BG, egui::Stroke::new(1.0, BORDER), egui::StrokeKind::Outside);

            // Dot position: flat=top, natural=middle, sharp=bottom
            let dot_y = match pos {
                music::PedalPos::Flat => track_top + dot_r + 2.0,
                music::PedalPos::Natural => track_top + track_h / 2.0,
                music::PedalPos::Sharp => track_top + track_h - dot_r - 2.0,
            };
            painter.circle_filled(egui::Pos2::new(cx, dot_y), dot_r, dot_color);

            // Click detection per pedal
            let hit_rect = egui::Rect::from_min_size(
                egui::Pos2::new(x, rect.top()),
                egui::Vec2::new(track_w, rect.height()),
            );
            if response.clicked() {
                if let Some(click_pos) = response.interact_pointer_pos() {
                    if hit_rect.contains(click_pos) {
                        pedals.cycle(i);
                    }
                }
            }

            x += track_w + spacing;
        }

        // Labels below
        let label_bot = rect.bottom() + 2.0;
        let labels = ["\u{2191} flat", "\u{2014} natural", "\u{2193} sharp"];
        let label_positions = [rect.left() + 30.0, rect.center().x, rect.right() - 30.0];
        for (text, lx) in labels.iter().zip(label_positions.iter()) {
            painter.text(
                egui::Pos2::new(*lx, label_bot),
                egui::Align2::CENTER_TOP,
                text,
                egui::FontId::monospace(9.0),
                TEXT_MUTED,
            );
        }
    }

    // ── Drill screen ───────────────────────────────────────────────

    fn draw_drill(&mut self, ctx: &egui::Context) {
        // Process metronome
        if let Some(ref mut drill) = self.drill {
            if !drill.is_finished() && drill.should_click() {
                self.metronome.click();
            }
        }

        // Process player input FIRST (so hitting the note before the beat counts)
        if self.config.input_mode == InputMode::Mic {
            self.process_mic_input();
        }
        if self.config.input_mode == InputMode::Manual {
            if ctx.input(|i| i.key_pressed(egui::Key::Space)) {
                if let Some(ref mut drill) = self.drill {
                    drill.advance_manual();
                }
            }
        }

        // Then auto-advance at the beat if the player didn't act in time
        if let Some(ref mut drill) = self.drill {
            drill.check_auto_advance();
        }

        // Check if drill finished
        if let Some(ref drill) = self.drill {
            if drill.is_finished() {
                self.screen = Screen::Results;
                self.stop_mic();
                return;
            }
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            let drill = self.drill.as_ref().unwrap();
            let key_root = self.pedals.key_root();

            // ── Header row: time remaining + stats ──
            ui.horizontal(|ui| {
                // Time remaining (large)
                ui.label(
                    egui::RichText::new(drill.progress_text())
                        .size(22.0)
                        .strong()
                        .color(TEXT_PRIMARY),
                );
                ui.label(
                    egui::RichText::new(format!("pos {}", drill.positions_completed))
                        .size(12.0)
                        .color(TEXT_MUTED),
                );

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if drill.level > 0 {
                        ui.label(egui::RichText::new(format!("Lv {}", drill.level)).size(13.0).color(ACCENT));
                    }
                    ui.label(egui::RichText::new(format!("{:.0} NPM", drill.notes_per_minute())).size(13.0).color(TEXT_MUTED));
                    if self.config.input_mode == InputMode::Mic {
                        ui.label(egui::RichText::new(format!("{:.0}%", drill.accuracy())).size(13.0).color(TEXT_MUTED));
                    }
                });
            });

            // ── Range indicator ──
            ui.horizontal(|ui| {
                ui.label(egui::RichText::new(format!("Range: {}", drill.range_display())).size(11.0).color(TEXT_MUTED));
                ui.label(egui::RichText::new(format!("\u{00b7} {}", self.pedals.key_name())).size(11.0).color(TEXT_MUTED));
            });

            // ── Timing bar ──
            let frac = drill.note_elapsed_fraction();
            let bar_color = if frac < 0.7 { GREEN } else if frac < 1.0 { YELLOW } else { RED };
            let bar_frac = frac.min(1.5) / 1.5;

            let (bar_rect, _) = ui.allocate_exact_size(
                egui::Vec2::new(ui.available_width(), 6.0),
                egui::Sense::hover(),
            );
            let painter = ui.painter_at(bar_rect);
            painter.rect_filled(bar_rect, 3.0, TRACK_BG);
            let filled = egui::Rect::from_min_size(
                bar_rect.min,
                egui::Vec2::new(bar_rect.width() * bar_frac, bar_rect.height()),
            );
            painter.rect_filled(filled, 3.0, bar_color);
            let target_x = bar_rect.left() + bar_rect.width() * (1.0 / 1.5);
            painter.line_segment(
                [
                    egui::Pos2::new(target_x, bar_rect.top() - 1.0),
                    egui::Pos2::new(target_x, bar_rect.bottom() + 1.0),
                ],
                egui::Stroke::new(2.0, egui::Color32::from_rgb(80, 80, 80)),
            );

            // ── Streak dots + pace ──
            ui.horizontal(|ui| {
                let streak = drill.streak;
                for i in 0..4u32 {
                    let (dot, color) = if i < streak {
                        ("\u{25cf}", GREEN)
                    } else {
                        ("\u{25cb}", egui::Color32::from_rgb(200, 200, 200))
                    };
                    ui.label(egui::RichText::new(dot).size(12.0).color(color));
                }
                ui.label(egui::RichText::new(format!("{:.1}s", drill.target_secs)).size(11.0).color(TEXT_MUTED));
                ui.label(egui::RichText::new(format!("({:.0} BPM)", drill.current_bpm())).size(11.0).color(TEXT_MUTED));
            });

            ui.add_space(2.0);

            // ── Staff card with hand position ──
            let staff_frame = egui::Frame::NONE
                .fill(CARD_BG)
                .corner_radius(12.0)
                .inner_margin(egui::Margin::same(8))
                .shadow(egui::epaint::Shadow {
                    offset: [0, 2],
                    blur: 12,
                    spread: 0,
                    color: egui::Color32::from_black_alpha(20),
                });

            staff_frame.show(ui, |ui| {
                let staff_height = 320.0;
                let (rect, _) = ui.allocate_exact_size(
                    egui::Vec2::new(ui.available_width(), staff_height),
                    egui::Sense::hover(),
                );

                let painter = ui.painter_at(rect);

                // Legend
                staff::draw_legend(&painter, rect, key_root);

                // Staff lines
                let staff_rect = egui::Rect::from_min_max(
                    egui::Pos2::new(rect.left(), rect.top() + 28.0),
                    rect.max,
                );
                staff::draw_staff(&painter, staff_rect);

                // Hand position (all notes with current highlighted)
                staff::draw_hand_position(
                    &painter,
                    staff_rect,
                    &drill.position.notes,
                    drill.finger_index,
                    key_root,
                );

                // Detected note (mic mode, dimmer overlay)
                if self.config.input_mode == InputMode::Mic {
                    if let Ok(state) = self.audio_state.lock() {
                        if let Some(ref note) = state.current_note {
                            staff::draw_note_on_staff(
                                &painter, staff_rect, note.midi, key_root,
                                egui::Color32::from_rgb(156, 163, 175), 0.4,
                            );
                        }
                    }
                }
            });

            // ── Note name + solfege below staff ──
            if let Some(target_midi) = drill.current_note() {
                let note_name = music::NOTE_NAMES[target_midi as usize % 12];
                let octave = (target_midi as i32 / 12) - 1;
                let degree = music::scale_degree(target_midi, key_root);
                let solfege = if degree > 0 { music::SOLFEGE[(degree - 1) as usize] } else { "?" };
                let finger = drill.finger_index + 1;

                ui.vertical_centered(|ui| {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(format!("{note_name}{octave}")).size(28.0).strong().color(TEXT_PRIMARY));
                        ui.label(egui::RichText::new(format!("({solfege})")).size(18.0).color(TEXT_MUTED));
                        ui.label(egui::RichText::new(format!("finger {finger}")).size(14.0).color(TEXT_MUTED));
                    });
                });
            }

            // Mic volume
            if self.config.input_mode == InputMode::Mic {
                if let Ok(state) = self.audio_state.lock() {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("Vol:").size(11.0).color(TEXT_MUTED));
                        let (bar_rect, _) = ui.allocate_exact_size(
                            egui::Vec2::new(100.0, 4.0), egui::Sense::hover(),
                        );
                        let p = ui.painter_at(bar_rect);
                        p.rect_filled(bar_rect, 2.0, TRACK_BG);
                        let filled = egui::Rect::from_min_size(
                            bar_rect.min,
                            egui::Vec2::new(bar_rect.width() * state.volume.clamp(0.0, 1.0), bar_rect.height()),
                        );
                        p.rect_filled(filled, 2.0, GREEN);
                    });
                }
                if let Some(err) = &self.audio_error {
                    ui.colored_label(RED, err);
                }
            }

            ui.add_space(6.0);

            // ── Controls ──
            ui.vertical_centered(|ui| {
                if self.config.input_mode == InputMode::Manual {
                    if action_btn(ui, "Next (Space)", 16.0) {
                        if let Some(ref mut d) = self.drill {
                            d.advance_manual();
                        }
                    }
                }
                ui.add_space(4.0);
                if secondary_btn(ui, "Abort") {
                    self.screen = Screen::Setup;
                    self.stop_mic();
                    self.drill = None;
                }
            });
        });
    }

    fn process_mic_input(&mut self) {
        let detected_midi = {
            let state = self.audio_state.lock().unwrap();
            state.current_note.as_ref().map(|n| n.midi)
        };

        if let Some(midi) = detected_midi {
            if self.last_detected == Some(midi) {
                self.confirm_count = self.confirm_count.saturating_add(1);
            } else {
                self.last_detected = Some(midi);
                self.confirm_count = 1;
            }
            if self.confirm_count >= 3 {
                if let Some(ref mut drill) = self.drill {
                    let result = drill.check_note(midi);
                    if result == NoteResult::Hit {
                        self.confirm_count = 0;
                        self.last_detected = None;
                    }
                }
            }
        } else {
            self.confirm_count = 0;
            self.last_detected = None;
        }
    }

    // ── Results screen ─────────────────────────────────────────────

    fn draw_results(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.add_space(40.0);
            ui.vertical_centered(|ui| {
                ui.label(
                    egui::RichText::new("DRILL COMPLETE")
                        .size(14.0)
                        .strong()
                        .color(TEXT_SECONDARY),
                );
                ui.add_space(16.0);

                if let Some(ref drill) = self.drill {
                    // Results card
                    let frame = egui::Frame::NONE
                        .fill(CARD_BG)
                        .corner_radius(12.0)
                        .inner_margin(egui::Margin::same(24))
                        .shadow(egui::epaint::Shadow {
                            offset: [0, 2],
                            blur: 12,
                            spread: 0,
                            color: egui::Color32::from_black_alpha(20),
                        });

                    frame.show(ui, |ui| {
                        ui.vertical_centered(|ui| {
                            ui.label(
                                egui::RichText::new(format!("{:.1}", drill.notes_per_minute()))
                                    .size(48.0)
                                    .strong()
                                    .color(TEXT_PRIMARY),
                            );
                            ui.label(egui::RichText::new("notes / min").size(14.0).color(TEXT_MUTED));

                            if drill.config.input_mode == InputMode::Mic {
                                ui.add_space(12.0);
                                ui.label(
                                    egui::RichText::new(format!("{:.0}%", drill.accuracy()))
                                        .size(32.0)
                                        .color(TEXT_PRIMARY),
                                );
                                ui.label(egui::RichText::new("accuracy").size(14.0).color(TEXT_MUTED));
                                ui.label(
                                    egui::RichText::new(format!("{} hits, {} misses", drill.hits, drill.misses))
                                        .size(13.0)
                                        .color(TEXT_MUTED),
                                );
                            }

                            ui.add_space(12.0);

                            // Positions + notes
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} positions \u{00b7} {} notes",
                                    drill.positions_completed, drill.hits,
                                ))
                                .size(14.0)
                                .color(TEXT_PRIMARY),
                            );

                            ui.add_space(8.0);

                            // Pace info
                            ui.label(
                                egui::RichText::new(format!(
                                    "Final pace: {:.1}s/note ({:.0} BPM)",
                                    drill.target_secs, drill.current_bpm()
                                ))
                                .size(13.0)
                                .color(TEXT_MUTED),
                            );
                            if drill.level > 0 {
                                ui.label(
                                    egui::RichText::new(format!(
                                        "Leveled up {} time{}",
                                        drill.level,
                                        if drill.level == 1 { "" } else { "s" }
                                    ))
                                    .size(13.0)
                                    .color(ACCENT),
                                );
                            }

                            // Range info
                            ui.label(
                                egui::RichText::new(format!(
                                    "Range: {}",
                                    drill.range_display(),
                                ))
                                .size(13.0)
                                .color(TEXT_MUTED),
                            );

                            ui.add_space(4.0);
                            ui.label(
                                egui::RichText::new(format!(
                                    "{} \u{00b7} {} fingers \u{00b7} Started {:.1}s",
                                    self.pedals.key_name(),
                                    drill.config.fingers,
                                    drill.config.difficulty.start_secs(),
                                ))
                                .size(12.0)
                                .color(TEXT_MUTED),
                            );
                        });
                    });
                }

                ui.add_space(24.0);

                if action_btn(ui, "Try Again", 16.0) {
                    self.start_drill();
                }
                ui.add_space(6.0);
                if secondary_btn(ui, "New Drill") {
                    self.drill = None;
                    self.screen = Screen::Setup;
                }
            });
        });
    }
}
