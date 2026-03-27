pub mod abc;
pub mod chord;
pub mod music;
pub mod notation;
pub mod voicing;

use eframe::egui;
use abc::{Score, ScoreEvent, load_embedded_scores};
use notation::{NotationLayout, render_score};

// ── Colors ──
const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241);
const CARD_BG: egui::Color32 = egui::Color32::WHITE;
const TEXT_PRIMARY: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const TEXT_MUTED: egui::Color32 = egui::Color32::from_rgb(136, 136, 136);
const ACCENT: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const BORDER: egui::Color32 = egui::Color32::from_rgb(204, 204, 204);

pub fn create_native_options() -> eframe::NativeOptions {
    eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1000.0, 600.0]),
        ..Default::default()
    }
}

pub fn run_app(options: eframe::NativeOptions) -> eframe::Result {
    eframe::run_native(
        "Harp Player",
        options,
        Box::new(|cc| {
            apply_style(&cc.egui_ctx);
            Ok(Box::new(PlayerApp::new()))
        }),
    )
}

fn apply_style(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing = egui::Vec2::new(6.0, 4.0);
    style.spacing.button_padding = egui::Vec2::new(10.0, 6.0);
    style.visuals.window_fill = BG;
    style.visuals.panel_fill = BG;
    style.visuals.extreme_bg_color = CARD_BG; // text edit background

    let r8 = egui::CornerRadius::same(8);
    style.visuals.widgets.inactive.corner_radius = r8;
    style.visuals.widgets.hovered.corner_radius = r8;
    style.visuals.widgets.active.corner_radius = r8;

    style.visuals.widgets.inactive.bg_fill = CARD_BG;
    style.visuals.widgets.inactive.weak_bg_fill = CARD_BG;
    style.visuals.widgets.inactive.bg_stroke = egui::Stroke::new(2.0, BORDER);
    style.visuals.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);

    style.visuals.widgets.hovered.bg_fill = egui::Color32::from_rgb(240, 240, 240);
    style.visuals.widgets.hovered.weak_bg_fill = egui::Color32::from_rgb(240, 240, 240);

    style.visuals.widgets.active.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.widgets.active.bg_stroke = egui::Stroke::new(2.0, ACCENT);

    style.visuals.selection.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.selection.stroke = egui::Stroke::new(2.0, ACCENT);

    ctx.set_style(style);
}

const ALL_KEYS: [&str; 12] = ["C","Db","D","Eb","E","F","F#","G","Ab","A","Bb","B"];
const MODE_NAMES: [&str; 7] = ["Ionian","Dorian","Phrygian","Lydian","Mixolydian","Aeolian","Locrian"];

struct PlayerApp {
    scores: Vec<Score>,
    current_score: usize,

    // Playback state
    playing: bool,
    current_beat: f32,
    tempo_bpm: f32,
    last_frame_time: Option<f64>,

    // Key/mode
    current_key: String,
    mode_offset: i32,

    // View
    scroll_offset: f32,
    playhead_fraction: f32,

    search_text: String,
    status: String,
}

impl PlayerApp {
    fn new() -> Self {
        let scores = load_embedded_scores();
        let status = format!("{} hymns loaded", scores.len());
        let tempo = scores.first().map(|s| s.tempo as f32).unwrap_or(120.0);
        let key = scores.first().map(|s| s.key.clone()).unwrap_or("C".into());

        Self {
            scores,
            current_score: 0,
            playing: false,
            current_beat: 0.0,
            tempo_bpm: tempo,
            last_frame_time: None,
            current_key: key,
            mode_offset: 0,
            scroll_offset: 0.0,
            playhead_fraction: 0.25,
            search_text: String::new(),
            status,
        }
    }

    fn current_events(&self) -> &[ScoreEvent] {
        self.scores.get(self.current_score)
            .map(|s| s.events.as_slice())
            .unwrap_or(&[])
    }

    fn total_beats(&self) -> f32 {
        let mut total = 0.0f32;
        for event in self.current_events() {
            match event {
                ScoreEvent::Note { beats, .. } => total += beats,
                ScoreEvent::Rest { beats, .. } => total += beats,
                ScoreEvent::Bar => {}
            }
        }
        total
    }

    fn reset_playback(&mut self) {
        self.current_beat = 0.0;
        self.scroll_offset = 0.0;
        self.last_frame_time = None;
        self.playing = false;
    }
}

impl eframe::App for PlayerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Advance playback
        if self.playing {
            let now = ctx.input(|i| i.time);
            if let Some(last) = self.last_frame_time {
                let dt = (now - last) as f32;
                let beats_per_sec = self.tempo_bpm / 60.0;
                self.current_beat += dt * beats_per_sec;

                let total = self.total_beats();
                if self.current_beat >= total {
                    self.current_beat = 0.0;
                    self.playing = false;
                }
            }
            self.last_frame_time = Some(now);
            ctx.request_repaint();
        }

        // ── Top controls (2 rows) ──
        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.add_space(2.0);

            // Row 1: Transport + BPM + Key/Mode
            ui.horizontal(|ui| {
                // Play/Pause
                let play_text = if self.playing { "||" } else { ">" };
                let play_btn = egui::Button::new(
                    egui::RichText::new(play_text).size(14.0)
                        .color(if self.playing { egui::Color32::WHITE } else { TEXT_PRIMARY })
                ).fill(if self.playing { ACCENT } else { CARD_BG });
                if ui.add(play_btn).clicked() {
                    self.playing = !self.playing;
                    if self.playing { self.last_frame_time = None; }
                }

                if ui.button(egui::RichText::new("[]").size(14.0)).clicked() {
                    self.reset_playback();
                }

                if ui.button(egui::RichText::new("-").size(14.0)).clicked() {
                    self.tempo_bpm = (self.tempo_bpm - 20.0).max(10.0);
                }
                ui.add(egui::DragValue::new(&mut self.tempo_bpm)
                    .range(10.0..=240.0).speed(1.0).suffix(" BPM"));
                if ui.button(egui::RichText::new("+").size(14.0)).clicked() {
                    self.tempo_bpm = (self.tempo_bpm + 20.0).min(240.0);
                }

                ui.separator();

                // Key selector
                let mut selected_key = self.current_key.clone();
                egui::ComboBox::from_id_salt("key_select")
                    .selected_text(&self.current_key)
                    .width(50.0)
                    .show_ui(ui, |ui| {
                        for k in ALL_KEYS {
                            ui.selectable_value(&mut selected_key, k.to_string(), k);
                        }
                    });
                if selected_key != self.current_key {
                    self.current_key = selected_key;
                }

                // Mode buttons
                for (i, name) in MODE_NAMES.iter().enumerate() {
                    let selected = self.mode_offset == i as i32;
                    let btn = if selected {
                        egui::Button::new(
                            egui::RichText::new(*name).size(10.0).color(egui::Color32::WHITE)
                        ).fill(ACCENT)
                    } else {
                        egui::Button::new(egui::RichText::new(*name).size(10.0))
                    };
                    if ui.add(btn).clicked() {
                        self.mode_offset = i as i32;
                    }
                }

                // Beat counter (right-aligned)
                let total = self.total_beats();
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(egui::RichText::new(
                        format!("{:.0}/{:.0}", self.current_beat, total)
                    ).monospace().size(11.0).color(TEXT_MUTED));
                });
            });

            // Row 2: Hymn selector + Search
            ui.horizontal(|ui| {
                if !self.scores.is_empty() {
                    let label = self.scores.get(self.current_score)
                        .map(|s| format!("{}. {}", s.number, s.title))
                        .unwrap_or_default();

                    let old_score = self.current_score;
                    egui::ComboBox::from_id_salt("score_select")
                        .selected_text(&label)
                        .width(ui.available_width() * 0.45)
                        .show_ui(ui, |ui| {
                            for (i, s) in self.scores.iter().enumerate() {
                                let text = format!("{}. {}", s.number, s.title);
                                ui.selectable_value(&mut self.current_score, i, &text);
                            }
                        });

                    if self.current_score != old_score {
                        self.reset_playback();
                        if let Some(s) = self.scores.get(self.current_score) {
                            self.tempo_bpm = s.tempo as f32;
                            self.current_key = s.key.clone();
                        }
                    }

                    // Search (matches hymn selector height)
                    let response = ui.add(
                        egui::TextEdit::singleline(&mut self.search_text)
                            .hint_text("Search...")
                            .desired_width(200.0)
                            .text_color(TEXT_PRIMARY)
                            .background_color(CARD_BG)
                            .font(egui::FontSelection::FontId(egui::FontId::proportional(16.0)))
                    );
                    if response.changed() && !self.search_text.is_empty() {
                        let q = self.search_text.to_lowercase();
                        if let Some(idx) = self.scores.iter().position(|s|
                            s.title.to_lowercase().contains(&q) || s.number.contains(&q)
                        ) {
                            self.current_score = idx;
                            self.reset_playback();
                            if let Some(s) = self.scores.get(self.current_score) {
                                self.tempo_bpm = s.tempo as f32;
                                self.current_key = s.key.clone();
                            }
                        }
                    }
                }
            });

            ui.add_space(2.0);
        });

        // ── Bottom status ──
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.colored_label(TEXT_MUTED, &self.status);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(s) = self.scores.get(self.current_score) {
                        ui.label(egui::RichText::new(
                            format!("Key: {}  Meter: {}/{}  Tempo: {} BPM",
                                s.key, s.meter_num, s.meter_den, s.tempo)
                        ).small().color(TEXT_MUTED));
                    }
                });
            });
        });

        // ── Central score view ──
        egui::CentralPanel::default().show(ctx, |ui| {
            let avail = ui.available_size();

            // Score card
            egui::Frame::NONE
                .fill(CARD_BG)
                .inner_margin(8.0)
                .corner_radius(12.0)
                .shadow(egui::epaint::Shadow {
                    offset: [0, 2],
                    blur: 12,
                    spread: 0,
                    color: egui::Color32::from_black_alpha(20),
                })
                .show(ui, |ui| {
                    let (rect, _response) = ui.allocate_exact_size(
                        egui::Vec2::new(avail.x - 16.0, avail.y - 24.0),
                        egui::Sense::click_and_drag(),
                    );

                    let painter = ui.painter_at(rect);
                    let layout = NotationLayout::new(rect.top() + 20.0, 40.0);

                    // Auto-scroll: keep playhead at playhead_fraction of view width
                    let view_width = rect.width();
                    let playhead_x_target = view_width * self.playhead_fraction;
                    self.scroll_offset = self.current_beat * 40.0 - playhead_x_target + layout.left_margin;
                    if self.scroll_offset < 0.0 { self.scroll_offset = 0.0; }

                    let events = self.current_events().to_vec();
                    let key_root = self.scores.get(self.current_score)
                        .map(|s| s.key_root).unwrap_or(0);

                    render_score(
                        &painter,
                        &layout,
                        &events,
                        self.scroll_offset,
                        self.current_beat,
                        view_width,
                        key_root,
                    );
                });
        });
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

#[cfg(test)]
mod tests {
    use super::abc::*;
    use super::music::*;
    
    #[test]
    fn check_joy_to_world() {
        let scores = load_embedded_scores();
        let joy = scores.iter().find(|s| s.title.contains("Joy to the World")).unwrap();
        eprintln!("Hymn: {} - {} (key: {})", joy.number, joy.title, joy.key);
        for (i, event) in joy.events.iter().take(20).enumerate() {
            match event {
                ScoreEvent::Note { melody_midi, rh_midi, lh_midi, beats, chord_name, .. } => {
                    let mel_name = midi_to_name(*melody_midi);
                    let rh_names: Vec<String> = rh_midi.iter().map(|&m| midi_to_name(m)).collect();
                    let lh_names: Vec<String> = lh_midi.iter().map(|&m| midi_to_name(m)).collect();
                    eprintln!("  [{:2}] melody={} ({}) beats={:.1} RH={:?} LH={:?} chord={:?}",
                        i, mel_name, melody_midi, beats, rh_names, lh_names, chord_name);
                }
                ScoreEvent::Rest { beats } => eprintln!("  [{:2}] REST beats={:.1}", i, beats),
                ScoreEvent::Bar => eprintln!("  [{:2}] BAR", i),
            }
        }
    }
}

#[cfg(test)]
mod tests2 {
    use super::abc::*;
    
    #[test]
    fn check_strings() {
        let scores = load_embedded_scores();
        let joy = scores.iter().find(|s| s.title.contains("Joy to the World")).unwrap();
        for (i, event) in joy.events.iter().take(8).enumerate() {
            if let ScoreEvent::Note { rh_strings, lh_strings, chord_name, .. } = event {
                if chord_name.is_some() {
                    eprintln!("[{}] chord={:?} rh_str={:?} lh_str={:?}", i, chord_name, rh_strings, lh_strings);
                }
            }
        }
    }
}
// Add to tests
#[cfg(test)]
mod collision_test {
    use super::abc::*;
    #[test]
    fn find_collisions() {
        let scores = load_embedded_scores();
        // Check first hymn (Blessed Jesus)
        let s = &scores[0];
        eprintln!("Checking: {} - {}", s.number, s.title);
        for (i, event) in s.events.iter().enumerate() {
            if let ScoreEvent::Note { rh_midi, lh_midi, rh_strings, lh_strings, chord_name, .. } = event {
                // Check MIDI collision
                for rm in rh_midi {
                    for lm in lh_midi {
                        if rm == lm {
                            eprintln!("  MIDI COLLISION at beat {}: RH and LH both have MIDI {} chord={:?}", i, rm, chord_name);
                        }
                    }
                }
                // Check string collision
                for rs in rh_strings {
                    for ls in lh_strings {
                        if rs == ls {
                            eprintln!("  STRING COLLISION at beat {}: RH and LH both have string {} chord={:?}", i, rs, chord_name);
                        }
                    }
                }
            }
        }
        // Also check Joy to the World
        if let Some(joy) = scores.iter().find(|s| s.title.contains("Joy")) {
            eprintln!("Checking: {} - {}", joy.number, joy.title);
            for (i, event) in joy.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_midi, lh_midi, chord_name, .. } = event {
                    for rm in rh_midi {
                        for lm in lh_midi {
                            if rm == lm {
                                eprintln!("  MIDI COLLISION at beat {}: both have MIDI {} chord={:?}", i, rm, chord_name);
                            }
                        }
                    }
                }
            }
        }
    }
}
