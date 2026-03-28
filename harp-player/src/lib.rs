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
    // Scale UI for tablet readability (notation sizes are set independently)
    ctx.set_pixels_per_point(1.5);

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
const MODE_NAMES: [&str; 7] = ["Ion","Dor","Phr","Lyd","Mix","Aeo","Loc"];

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

    // Voicing
    rh_fingers: i32,
    lh_fingers: i32,

    // View
    scroll_offset: f32,
    playhead_fraction: f32,

    search_text: String,
    status: String,
    recent_hymns: Vec<usize>,
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
            rh_fingers: 4,
            lh_fingers: 4,
            scroll_offset: 0.0,
            playhead_fraction: 0.25,
            search_text: String::new(),
            status,
            recent_hymns: Vec::new(),
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

    fn select_hymn(&mut self, idx: usize) {
        if self.current_score != idx {
            // Add to recent list (deduplicate, keep last 10)
            self.recent_hymns.retain(|&i| i != idx);
            self.recent_hymns.insert(0, idx);
            self.recent_hymns.truncate(10);
            self.current_score = idx;
            self.reset_playback();
            if let Some(s) = self.scores.get(idx) {
                self.tempo_bpm = s.tempo as f32;
                self.current_key = s.key.clone();
            }
        }
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

        // ── Top controls: 3 columns ──
        // ── Top controls ──
        // Row 1: [Hymn ▾]              [Search]
        // Row 2: [⏮ ▶ ⏹]              [🕒 Recent ▾]        ⊙ Circle
        // Row 3: [G ▾] 120 BPM         RH [4] LH [4]       of 5ths
        // Row 4: [Ion][Dor][Phr][Lyd][Mix][Aeo][Loc]
        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                // ── Col 1: Circle of fifths ──
                let selected_pc = music::key_to_pc(&self.current_key);
                let diameter = 100.0;
                let cof_size = diameter * 0.5;
                let (cof_rect, _) = ui.allocate_exact_size(
                    egui::Vec2::splat(diameter),
                    egui::Sense::hover(),
                );
                let cof_center = cof_rect.center();
                let cof_painter = ui.painter_at(cof_rect);
                const COF_ORDER: [i32; 12] = [0, 7, 2, 9, 4, 11, 6, 1, 8, 3, 10, 5];
                const COF_LABELS: [&str; 12] = ["C","G","D","A","E","B","F#","Db","Ab","Eb","Bb","F"];
                let outer_r = cof_size - 2.0;
                let inner_r = cof_size * 0.55;
                cof_painter.circle_stroke(cof_center, outer_r, egui::Stroke::new(1.0, BORDER));
                for (i, (&pc, &label)) in COF_ORDER.iter().zip(COF_LABELS.iter()).enumerate() {
                    let angle = std::f32::consts::TAU * i as f32 / 12.0 - std::f32::consts::FRAC_PI_2;
                    let label_r = (outer_r + inner_r) * 0.5;
                    let lx = cof_center.x + label_r * angle.cos();
                    let ly = cof_center.y + label_r * angle.sin();
                    let is_current = pc == selected_pc;
                    let color = if is_current { ACCENT } else { TEXT_MUTED };
                    let font_size = if is_current { 10.0 } else { 8.0 };
                    if is_current {
                        cof_painter.circle_filled(
                            egui::Pos2::new(lx, ly), 10.0,
                            egui::Color32::from_rgb(219, 234, 254),
                        );
                    }
                    cof_painter.text(
                        egui::Pos2::new(lx, ly),
                        egui::Align2::CENTER_CENTER,
                        label,
                        egui::FontId::proportional(font_size),
                        color,
                    );
                }

                // ── Col 2: controls ──
                ui.vertical(|ui| {
                    ui.spacing_mut().item_spacing.y = 2.0;

                    // Row 1: [Hymn ▾] [Search]
                    ui.horizontal(|ui| {
                        if !self.scores.is_empty() {
                            let label = self.scores.get(self.current_score)
                                .map(|s| format!("{}. {}", s.number, s.title))
                                .unwrap_or_default();
                            let old_score = self.current_score;
                            egui::ComboBox::from_id_salt("score_select")
                                .selected_text(&label)
                                .width(150.0)
                                .show_ui(ui, |ui| {
                                    for (i, s) in self.scores.iter().enumerate() {
                                        let text = format!("{}. {}", s.number, s.title);
                                        ui.selectable_value(&mut self.current_score, i, &text);
                                    }
                                });
                            if self.current_score != old_score {
                                self.select_hymn(self.current_score);
                            }
                            let response = ui.add(
                                egui::TextEdit::singleline(&mut self.search_text)
                                    .hint_text("Search...")
                                    .desired_width(ui.available_width() - 4.0)
                                    .min_size(egui::Vec2::new(0.0, 24.0))
                                    .text_color(TEXT_PRIMARY)
                                    .background_color(CARD_BG)
                                    .font(egui::FontSelection::FontId(egui::FontId::proportional(13.0)))
                            );
                            if response.changed() && !self.search_text.is_empty() {
                                let q = self.search_text.to_lowercase();
                                if let Some(idx) = self.scores.iter().position(|s|
                                    s.title.to_lowercase().contains(&q) || s.number.contains(&q)
                                ) {
                                    self.select_hymn(idx);
                                }
                            }
                        }
                    });

                    // Row 2: [⏮ ▶ ⏹] [G▾] 120BPM RH[4] LH[4] [🕒 Recent ▾]
                    ui.horizontal(|ui| {
                        if ui.button(egui::RichText::new("\u{23EE}").size(14.0)).clicked() {
                            self.reset_playback();
                        }
                        let play_icon = if self.playing { "\u{23F8}" } else { "\u{25B6}" };
                        let play_btn = egui::Button::new(
                            egui::RichText::new(play_icon).size(14.0)
                                .color(if self.playing { egui::Color32::WHITE } else { TEXT_PRIMARY })
                        ).fill(if self.playing { ACCENT } else { CARD_BG });
                        if ui.add(play_btn).clicked() {
                            self.playing = !self.playing;
                            if self.playing { self.last_frame_time = None; }
                        }
                        if ui.button(egui::RichText::new("\u{23F9}").size(14.0)).clicked() {
                            self.reset_playback();
                        }

                        ui.separator();

                        let mut selected_key = self.current_key.clone();
                        egui::ComboBox::from_id_salt("key_select")
                            .selected_text(egui::RichText::new(&self.current_key).size(11.0))
                            .width(40.0)
                            .show_ui(ui, |ui| {
                                for k in ALL_KEYS {
                                    ui.selectable_value(&mut selected_key, k.to_string(), k);
                                }
                            });
                        if selected_key != self.current_key {
                            self.current_key = selected_key;
                        }
                        ui.add(egui::DragValue::new(&mut self.tempo_bpm)
                            .range(10.0..=240.0).speed(1.0).suffix(" BPM"));

                        ui.separator();

                        ui.label(egui::RichText::new("RH").size(10.0).color(TEXT_MUTED));
                        ui.add(egui::DragValue::new(&mut self.rh_fingers).range(1..=4).speed(0.1));
                        ui.label(egui::RichText::new("LH").size(10.0).color(TEXT_MUTED));
                        ui.add(egui::DragValue::new(&mut self.lh_fingers).range(1..=4).speed(0.1));

                    });

                    // Row 3: Mode buttons + Recent
                    ui.horizontal(|ui| {
                        for (i, name) in MODE_NAMES.iter().enumerate() {
                            let selected = self.mode_offset == i as i32;
                            let btn = if selected {
                                egui::Button::new(
                                    egui::RichText::new(*name).size(11.0).color(egui::Color32::WHITE)
                                ).fill(ACCENT)
                            } else {
                                egui::Button::new(egui::RichText::new(*name).size(11.0))
                            };
                            if ui.add(btn).clicked() {
                                self.mode_offset = i as i32;
                            }
                        }

                        ui.separator();

                        let mut recent_pick = self.current_score;
                        let recent_label = if self.recent_hymns.is_empty() { "Recent" } else { "\u{1F552} Recent" };
                        egui::ComboBox::from_id_salt("recent_select")
                            .selected_text(recent_label)
                            .width(ui.available_width() - 4.0)
                            .show_ui(ui, |ui| {
                                if self.recent_hymns.is_empty() {
                                    ui.label(egui::RichText::new("No recent hymns").color(TEXT_MUTED));
                                }
                                for &idx in &self.recent_hymns {
                                    if let Some(s) = self.scores.get(idx) {
                                        let text = format!("{}. {}", s.number, s.title);
                                        ui.selectable_value(&mut recent_pick, idx, &text);
                                    }
                                }
                            });
                        if recent_pick != self.current_score {
                            self.select_hymn(recent_pick);
                        }
                    });
                });

            });
            ui.add_space(2.0);
        });

        // ── Bottom: progress bar + status ──
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            // Progress bar (full width, draggable)
            let total = self.total_beats().max(1.0);
            let mut progress = self.current_beat / total;
            let slider = egui::Slider::new(&mut progress, 0.0..=1.0)
                .show_value(false)
                .trailing_fill(true);
            if ui.add(slider).changed() {
                self.current_beat = progress * total;
            }

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
                    let layout = NotationLayout::new(rect.top() + 20.0, 50.0);

                    // Auto-scroll: keep playhead at playhead_fraction of view width
                    let view_width = rect.width();
                    let playhead_x_target = view_width * self.playhead_fraction;
                    self.scroll_offset = self.current_beat * 60.0 - playhead_x_target + layout.left_margin;
                    if self.scroll_offset < 0.0 { self.scroll_offset = 0.0; }

                    let events = self.current_events().to_vec();
                    // Compute key_root from the user-selected key and mode
                    let selected_pc = music::key_to_pc(&self.current_key);
                    let key_root = ((selected_pc - music::MAJOR_SCALE[self.mode_offset as usize]) % 12 + 12) % 12;

                    render_score(
                        &painter,
                        &layout,
                        &events,
                        self.scroll_offset,
                        self.current_beat,
                        view_width,
                        key_root,
                        self.rh_fingers as usize,
                        self.lh_fingers as usize,
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
#[cfg(test)]
mod debug_hymn22 {
    use super::abc::*;
    use super::music::*;

    #[test]
    fn debug_blessed_jesus() {
        let scores = load_embedded_scores();
        let s = &scores[0];
        eprintln!("Hymn: {} - {} key={} key_root={}", s.number, s.title, s.key, s.key_root);
        
        let str_offset = STRING_NUMBER_OFFSET;
        
        for (i, event) in s.events.iter().take(10).enumerate() {
            if let ScoreEvent::Note { rh_strings, lh_strings, rh_midi, lh_midi, chord_name, is_chord_change, .. } = event {
                if *is_chord_change {
                    eprintln!("[{i}] chord={chord_name:?}");
                    eprintln!("  RH strings: {rh_strings:?}  RH midi: {rh_midi:?}");
                    for &rs in rh_strings {
                        let abs_str = rs + str_offset - 1;
                        let display = harp_string_to_midi(abs_str, s.key_root) + 12;
                        let pc = display % 12;
                        let oct = display / 12 - 1;
                        let names = ["C","C#","D","D#","E","F","F#","G","G#","A","A#","B"];
                        eprintln!("    str {} -> abs {} -> display midi {} = {}{}", rs, abs_str, display, names[pc as usize], oct);
                    }
                    eprintln!("  LH strings: {lh_strings:?}  LH midi: {lh_midi:?}");
                    for &ls in lh_strings {
                        let abs_str = ls + str_offset - 1;
                        let display = harp_string_to_midi(abs_str, s.key_root) + 12;
                        let pc = display % 12;
                        let oct = display / 12 - 1;
                        let names = ["C","C#","D","D#","E","F","F#","G","G#","A","A#","B"];
                        eprintln!("    str {} -> abs {} -> display midi {} = {}{}", ls, abs_str, display, names[pc as usize], oct);
                    }
                }
            }
        }
    }
}
