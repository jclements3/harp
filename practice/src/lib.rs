pub mod audio;
pub mod drill;
pub mod music;
pub mod staff;

use std::sync::{Arc, Mutex};
use eframe::egui;
use drill::{DrillConfig, DrillMode, DrillSession, NoteResult, load_progress, save_progress};
use music::PedalState;

const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241);
const CARD: egui::Color32 = egui::Color32::WHITE;
const PRI: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const MUT: egui::Color32 = egui::Color32::from_rgb(136, 136, 136);
const ACC: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const GRN: egui::Color32 = egui::Color32::from_rgb(40, 167, 69);
const BDR: egui::Color32 = egui::Color32::from_rgb(204, 204, 204);

const KEYS: [&str; 12] = ["C","Db","D","Eb","E","F","F#","G","Ab","A","Bb","B"];
const MODES: [&str; 7] = ["Ionian","Dorian","Phrygian","Lydian","Mixolydian","Aeolian","Locrian"];
const MAJ: [i32; 7] = [0, 2, 4, 5, 7, 9, 11];

pub fn create_native_options() -> eframe::NativeOptions {
    let icon = load_icon();
    let mut vp = egui::ViewportBuilder::default().with_inner_size([1000.0, 600.0]);
    if let Some(i) = icon { vp = vp.with_icon(i); }
    eframe::NativeOptions { viewport: vp, ..Default::default() }
}
fn load_icon() -> Option<egui::IconData> {
    let img = image::load_from_memory(include_bytes!("../icon.png")).ok()?.into_rgba8();
    Some(egui::IconData { width: img.width(), height: img.height(), rgba: img.into_raw() })
}
pub fn run_app(opt: eframe::NativeOptions) -> eframe::Result {
    eframe::run_native("Practice", opt, Box::new(|cc| { apply_style(&cc.egui_ctx); Ok(Box::new(App::new())) }))
}
#[cfg(target_os = "android")]
#[unsafe(no_mangle)]
fn android_main(app: winit::platform::android::activity::AndroidApp) {
    android_logger::init_once(android_logger::Config::default().with_max_level(log::LevelFilter::Info));
    run_app(eframe::NativeOptions { android_app: Some(app), renderer: eframe::Renderer::Glow, ..Default::default() }).unwrap();
}
fn apply_style(ctx: &egui::Context) {
    ctx.set_pixels_per_point(1.5);
    let mut s = (*ctx.style()).clone();
    s.spacing.item_spacing = egui::Vec2::new(6.0, 4.0);
    s.spacing.button_padding = egui::Vec2::new(10.0, 6.0);
    s.visuals.window_fill = BG; s.visuals.panel_fill = BG; s.visuals.extreme_bg_color = CARD;
    let r = egui::CornerRadius::same(8);
    s.visuals.widgets.inactive.corner_radius = r; s.visuals.widgets.hovered.corner_radius = r; s.visuals.widgets.active.corner_radius = r;
    s.visuals.widgets.inactive.bg_fill = CARD; s.visuals.widgets.inactive.weak_bg_fill = CARD;
    s.visuals.widgets.inactive.bg_stroke = egui::Stroke::new(2.0, BDR);
    s.visuals.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, PRI);
    s.visuals.widgets.hovered.bg_fill = egui::Color32::from_rgb(240,240,240);
    s.visuals.widgets.active.bg_fill = egui::Color32::from_rgb(219,234,254);
    s.visuals.widgets.active.bg_stroke = egui::Stroke::new(2.0, ACC);
    s.visuals.selection.bg_fill = egui::Color32::from_rgb(219,234,254);
    ctx.set_style(s);
}
fn kpc(k: &str) -> i32 { match k {"C"=>0,"Db"|"C#"=>1,"D"=>2,"Eb"=>3,"E"=>4,"F"=>5,"F#"|"Gb"=>6,"G"=>7,"Ab"=>8,"A"=>9,"Bb"=>10,"B"|"Cb"=>11,_=>0} }

struct App {
    drill: Option<DrillSession>, config: DrillConfig, progress: drill::Progress,
    key: String, mode_off: i32, paused: bool,
    use_mic: bool, score_mis: bool, random_fingers: bool,
    audio_st: Arc<Mutex<audio::AudioState>>, _mic: Option<cpal::Stream>,
    cc: u8, last_det: Option<u8>,
}
impl App {
    fn new() -> Self {
        Self {
            drill: None, config: DrillConfig::default(), progress: load_progress(),
            key: "Eb".into(), mode_off: 0, paused: false,
            use_mic: false, score_mis: false, random_fingers: false,
            audio_st: Arc::new(Mutex::new(audio::AudioState::default())), _mic: None,
            cc: 0, last_det: None,
        }
    }
    fn pedals(&self) -> PedalState {
        let mut p = PedalState::default();
        for &(k,n) in &[("G",1),("D",2),("A",3),("E",4),("B",5),("F#",6)] { if k==self.key { for _ in 0..n { p.add_sharp(); } return p; } }
        for &(k,n) in &[("F",1),("Bb",2),("Eb",3),("Ab",4),("Db",5)] { if k==self.key { for _ in 0..n { p.add_flat(); } return p; } }
        p
    }
    fn start(&mut self) {
        self.drill = Some(DrillSession::new(self.config.clone(), &self.pedals(), &self.progress));
        self.paused = false; self.cc = 0; self.last_det = None;
        if self.use_mic && self._mic.is_none() { if let Ok(s) = audio::start_capture(Arc::clone(&self.audio_st)) { self._mic = Some(s); } }
    }
    fn stop(&mut self) {
        if let Some(ref d) = self.drill { self.progress = d.save_to_progress(); }
        self.drill = None; self._mic = None;
    }
    fn mic_tick(&mut self) {
        let midi = self.audio_st.lock().ok().and_then(|s| s.current_note.as_ref().map(|n| n.midi));
        let Some(midi) = midi else { self.cc = 0; self.last_det = None; return; };
        if self.last_det == Some(midi) { self.cc = self.cc.saturating_add(1); } else { self.last_det = Some(midi); self.cc = 1; }
        if self.cc >= 3 {
            if let Some(ref mut d) = self.drill {
                let r = d.check_note(midi, self.score_mis);
                if r == NoteResult::Hit { self.cc = 0; self.last_det = None; }
            }
        }
    }
    fn running(&self) -> bool { self.drill.as_ref().map_or(false, |d| !d.is_finished()) }
    fn adjust_pace(&mut self, delta: f32) {
        self.progress.target_secs = (self.progress.target_secs + delta).clamp(0.4, 12.0);
        save_progress(&self.progress);
        if let Some(ref mut d) = self.drill { d.target_secs = self.progress.target_secs; }
    }
}

fn pill(ui: &mut egui::Ui, lbl: &str, sel: bool) -> bool {
    let t = if sel { egui::RichText::new(lbl).size(11.0).color(egui::Color32::WHITE) } else { egui::RichText::new(lbl).size(11.0) };
    ui.add(egui::Button::new(t).fill(if sel { ACC } else { CARD }).stroke(if sel { egui::Stroke::new(2.0,ACC) } else { egui::Stroke::new(1.5,BDR) }).corner_radius(6.0)).clicked()
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.request_repaint();
        let run = self.running();

        // Input
        if run && !self.paused {
            if self.use_mic { self.mic_tick(); }
            if ctx.input(|i| i.key_pressed(egui::Key::Space)) { if let Some(ref mut d) = self.drill { d.advance_manual(); } }
            if let Some(ref mut d) = self.drill { d.check_auto_advance(); }
        }
        if run && ctx.input(|i| i.key_pressed(egui::Key::P)) { self.paused = !self.paused; }
        if let Some(ref d) = self.drill { if d.is_finished() && d.challenges_completed > 0 { self.progress = d.save_to_progress(); } }

        let kr = { let pc = kpc(&self.key); (((pc - MAJ[self.mode_off as usize]) % 12 + 12) % 12) as u8 };

        // ── Row 4 (very bottom): Score + safe margin ──
        egui::TopBottomPanel::bottom("score").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if let Some(ref d) = self.drill {
                    ui.label(egui::RichText::new(format!("Lv {}", d.level)).size(14.0).strong().color(ACC));
                    ui.separator();
                    ui.label(egui::RichText::new(format!("{:.0} NPM", d.notes_per_minute())).size(14.0).strong().color(PRI));
                    ui.separator();
                    for i in 0..4u32 {
                        let (dot, col) = if i < d.streak { ("\u{25cf}", GRN) } else { ("\u{25cb}", egui::Color32::from_rgb(200,200,200)) };
                        ui.label(egui::RichText::new(dot).size(16.0).color(col));
                    }
                } else {
                    ui.label(egui::RichText::new(format!("Lv {}", self.progress.level)).size(14.0).strong().color(ACC));
                    if self.progress.best_npm > 0.0 {
                        ui.separator();
                        ui.label(egui::RichText::new(format!("Best {:.0} NPM", self.progress.best_npm)).size(12.0).color(MUT));
                    }
                }
            });
            // Safe margin for Android nav bar
            ui.add_space(16.0);
        });

        // ── Row 3: Key + Church modes ──
        egui::TopBottomPanel::bottom("modes").show(ctx, |ui| {
            ui.horizontal(|ui| {
                let mut sk = self.key.clone();
                egui::ComboBox::from_id_salt("k").selected_text(egui::RichText::new(&self.key).size(13.0).strong().color(PRI)).width(42.0)
                    .show_ui(ui, |ui| { for k in KEYS { ui.selectable_value(&mut sk, k.to_string(), k); } });
                if sk != self.key { self.key = sk; }

                ui.separator();

                ui.spacing_mut().button_padding = egui::Vec2::new(4.0, 2.0);
                for (i, nm) in MODES.iter().enumerate() { if pill(ui, nm, self.mode_off == i as i32) { self.mode_off = i as i32; } }
                ui.spacing_mut().button_padding = egui::Vec2::new(10.0, 6.0);

                ui.separator();
                if ui.add(egui::Button::new(egui::RichText::new("Calibrate").size(10.0))).clicked() {
                    // TODO: mic calibration
                }
            });
        });

        // ── Row 2: Controls ──
        egui::TopBottomPanel::bottom("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                // Transport
                if run {
                    let (ic, fl, cl) = if self.paused { ("\u{25B6}", CARD, PRI) } else { ("\u{23F8}", ACC, egui::Color32::WHITE) };
                    if ui.add(egui::Button::new(egui::RichText::new(ic).size(16.0).color(cl)).fill(fl).min_size(egui::Vec2::new(34.0, 28.0))).clicked() { self.paused = !self.paused; }
                    if ui.add(egui::Button::new(egui::RichText::new("\u{23F9}").size(16.0)).min_size(egui::Vec2::new(34.0, 28.0))).clicked() { self.stop(); }
                } else {
                    if ui.add(egui::Button::new(egui::RichText::new("\u{25B6}").size(16.0).color(egui::Color32::WHITE)).fill(ACC).min_size(egui::Vec2::new(34.0, 28.0))).clicked() { self.start(); }
                }

                ui.separator();

                // Fingers
                ui.label(egui::RichText::new("Fingers").size(10.0).color(MUT));
                for f in 1..=4u32 {
                    let sel = self.config.chord_size == f as usize;
                    if pill(ui, &f.to_string(), sel) { self.config.chord_size = f as usize; }
                }

                ui.separator();

                // Random
                ui.checkbox(&mut self.random_fingers, egui::RichText::new("Random").size(11.0));

                ui.separator();

                // Mic + Score
                ui.checkbox(&mut self.use_mic, egui::RichText::new("Mic").size(11.0));
                if self.use_mic {
                    ui.checkbox(&mut self.score_mis, egui::RichText::new("Score").size(11.0));
                }

                ui.separator();

                // Intervals / Chords
                if pill(ui, "Intervals", self.config.mode == DrillMode::Intervals) { self.config.mode = DrillMode::Intervals; }
                if pill(ui, "Chords", self.config.mode == DrillMode::Chords) { self.config.mode = DrillMode::Chords; }
            });
        });

        // ── Row 1 (top): Staff fills everything above ──
        egui::CentralPanel::default().show(ctx, |ui| {
            let avail = ui.available_size();
            let mut pace_drag = 0.0f32;

            egui::Frame::NONE.fill(CARD).inner_margin(6.0).corner_radius(12.0)
                .shadow(egui::epaint::Shadow { offset: [0, 2], blur: 12, spread: 0, color: egui::Color32::from_black_alpha(20) })
                .show(ui, |ui| {
                    let (rect, response) = ui.allocate_exact_size(egui::Vec2::new(avail.x - 12.0, avail.y - 12.0), egui::Sense::click_and_drag());
                    let p = ui.painter_at(rect);

                    if response.dragged() { pace_drag = response.drag_delta().y; }

                    if let Some(ref drill) = self.drill {
                        if !drill.is_finished() {
                            let t = if self.paused { drill.current_challenge().map_or(0.0, |tc| tc.start_time) } else { drill.elapsed() };
                            staff::render_scrolling(&p, rect, &drill.timeline, drill.current_idx, t, kr);
                            if self.use_mic { if let Ok(st) = self.audio_st.lock() { if let Some(ref n) = st.current_note { staff::draw_detected_note(&p, rect, n.midi, kr); } } }
                        } else {
                            p.text(rect.center(), egui::Align2::CENTER_CENTER,
                                format!("{:.1} NPM \u{00b7} Lv {} \u{00b7} Press \u{25B6}", drill.notes_per_minute(), drill.level),
                                egui::FontId::proportional(18.0), ACC);
                        }
                    } else {
                        let ly = staff::Layout::new(&rect);
                        staff::draw_staff_lines(&p, &ly, rect.width());
                        p.text(rect.center(), egui::Align2::CENTER_CENTER,
                            "Press \u{25B6}", egui::FontId::proportional(18.0), MUT);
                    }
                });

            if pace_drag.abs() > 0.5 { self.adjust_pace(pace_drag * 0.01); }
        });
    }
}
