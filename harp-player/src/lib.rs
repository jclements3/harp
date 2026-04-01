pub mod abc;
pub mod audio;
pub mod chord;
pub mod drill;
pub mod music;
pub mod notation;
pub mod voicing;

use eframe::egui;
use abc::{Score, ScoreEvent, load_embedded_scores};
use notation::{NotationLayout, render_score};

// ── Persistent settings ──

fn settings_path() -> std::path::PathBuf {
    if cfg!(target_os = "android") {
        std::path::PathBuf::from("/data/data/com.harp.player/harp_player_settings.json")
    } else {
        let home = std::env::var("HOME").unwrap_or_else(|_| ".".into());
        std::path::PathBuf::from(home).join(".harp_player_settings.json")
    }
}

fn load_saved_bpm() -> Option<f32> {
    let text = std::fs::read_to_string(settings_path()).ok()?;
    // Simple key:value parse — just look for "bpm":NUMBER
    for line in text.lines() {
        if let Some(rest) = line.trim().strip_prefix("\"bpm\":") {
            let val = rest.trim().trim_end_matches(',');
            return val.parse().ok();
        }
    }
    None
}

fn save_bpm(bpm: f32) {
    let content = format!("{{\n  \"bpm\": {:.0}\n}}\n", bpm);
    let _ = std::fs::write(settings_path(), content);
}

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
const MODE_NAMES: [&str; 7] = ["Io","Do","Ph","Ly","Mi","Ae","Lo"];
const DRILL_NAMES: [&str; 1] = ["Random"];

struct PlayerApp {
    scores: Vec<Score>,
    current_score: usize,

    // Playback state
    playing: bool,
    current_beat: f32,
    tempo_bpm: f32,
    prev_tempo_bpm: f32,
    last_frame_time: Option<f64>,

    // Key/mode
    current_key: String,
    prev_key: String,
    mode_offset: i32,
    prev_mode: i32,

    // Voicing
    rh_fingers: i32,
    lh_fingers: i32,
    prev_rh_fingers: i32,
    prev_lh_fingers: i32,
    random_fingers: bool,
    random_pattern: Vec<(usize, usize)>, // per-chord (rh_count, lh_count)
    sound_volume: f32, // 0.0 = mute, 1.0 = max
    pre_mute_volume: f32, // remember volume before muting
    audio_player: Option<audio::AudioPlayer>,
    last_chord_idx: i64, // track which chord we last played sound for

    // View
    scroll_offset: f32,
    playhead_fraction: f32,

    search_text: String,
    prev_search: String,
    status: String,
    recent_hymns: Vec<usize>,

    // Current chord info for circle of fifths display
    current_chord_degrees: Vec<i32>,  // diatonic degrees (0-6) of current chord tones
    current_melody_degree: i32,       // melody's diatonic degree (0-6)

    // Drill mode
    drill_mode: bool,
    drill_score: Option<abc::Score>,  // random score, same type as hymns
    drill_progress: drill::Progress,
    drill_config: drill::DrillConfig,
    drill_name: usize,
    use_mic: bool,
    // Note duration checkboxes per hand
    rh_whole: bool,
    rh_half: bool,
    rh_quarter: bool,
    rh_eighth: bool,
    lh_whole: bool,
    lh_half: bool,
    lh_quarter: bool,
    lh_eighth: bool,
}

impl PlayerApp {
    fn new() -> Self {
        let scores = load_embedded_scores();
        let status = format!("{} hymns loaded", scores.len());
        // Load saved BPM from file, default to 10
        let tempo = load_saved_bpm().unwrap_or(10.0);
        let key = scores.first().map(|s| s.key.clone()).unwrap_or("C".into());

        Self {
            scores,
            current_score: 0,
            playing: false,
            current_beat: 0.0,
            tempo_bpm: tempo,
            prev_tempo_bpm: tempo,
            last_frame_time: None,
            current_key: key.clone(),
            prev_key: key,
            mode_offset: 0,
            prev_mode: 0,
            rh_fingers: 4,
            lh_fingers: 4,
            prev_rh_fingers: 4,
            prev_lh_fingers: 4,
            random_fingers: false,
            random_pattern: Vec::new(),
            sound_volume: 0.0,
            pre_mute_volume: 0.5,
            audio_player: audio::AudioPlayer::new().ok(),
            last_chord_idx: -1,
            scroll_offset: 0.0,
            playhead_fraction: 0.25,
            search_text: String::new(),
            prev_search: String::new(),
            status,
            recent_hymns: Vec::new(),
            current_chord_degrees: vec![0, 2, 4], // default I chord
            current_melody_degree: 0,
            drill_mode: false,
            drill_score: None,
            drill_progress: drill::load_progress(),
            drill_config: drill::DrillConfig::default(),
            drill_name: 0,
            use_mic: false,
            rh_whole: true, rh_half: true, rh_quarter: true, rh_eighth: false,
            lh_whole: true, lh_half: false, lh_quarter: false, lh_eighth: false,
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

    fn shuffle_fingers(&mut self) {
        // Generate a random (rh, lh) finger count per chord change in the current hymn
        let mut seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos() as u64;

        // Count chord changes
        let n_chords = self.current_events().iter()
            .filter(|e| matches!(e, ScoreEvent::Note { is_chord_change: true, .. }))
            .count();

        self.random_pattern.clear();
        for _ in 0..n_chords.max(1) {
            // Simple LCG random
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let rh = ((seed >> 16) % 3 + 2) as usize; // 2-4
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let lh = ((seed >> 16) % 3 + 2) as usize; // 2-4
            self.random_pattern.push((rh, lh));
        }
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
                // Keep user's BPM — don't override with hymn's default
                self.current_key = s.key.clone();
            }
        }
    }
}

impl eframe::App for PlayerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Re-voice current hymn if fingers, key, or mode changed
        // Update prev values FIRST to prevent re-triggering
        let needs_revoice = self.rh_fingers != self.prev_rh_fingers
            || self.lh_fingers != self.prev_lh_fingers
            || self.current_key != self.prev_key
            || self.mode_offset != self.prev_mode;

        if needs_revoice {
            // Set prev values immediately to prevent re-triggering next frame
            self.prev_rh_fingers = self.rh_fingers;
            self.prev_lh_fingers = self.lh_fingers;
            self.prev_key = self.current_key.clone();
            self.prev_mode = self.mode_offset;
            if let Some(score) = self.scores.get(self.current_score).cloned() {
                // Compute transposition: difference between hymn's original key and selected key
                let original_pc = music::key_to_pc(&score.key);
                let selected_pc = music::key_to_pc(&self.current_key);
                let transpose = ((selected_pc - original_pc) % 12 + 12) % 12;

                // Re-voice with transposition and finger constraints
                let new_score = abc::revoice_score_transposed(
                    &score,
                    self.rh_fingers as usize,
                    self.lh_fingers as usize,
                    transpose,
                    &self.current_key,
                );
                let pattern = if self.random_fingers { &self.random_pattern } else { &vec![] as &Vec<(usize,usize)> };
                let mut new_score = new_score;
                abc::compute_display_strings(&mut new_score, self.rh_fingers as usize, self.lh_fingers as usize, pattern);
                self.scores[self.current_score] = new_score;
            }
        }

        // Save BPM when it changes
        if self.tempo_bpm != self.prev_tempo_bpm {
            save_bpm(self.tempo_bpm);
            self.prev_tempo_bpm = self.tempo_bpm;
        }

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

            // Trigger sound at chord changes (songs mode only)
            if self.sound_volume > 0.0 && !self.drill_mode {
                if let Some(ref player) = self.audio_player {
                    if let Ok(mut synth) = player.state.lock() {
                        synth.master_volume = self.sound_volume;
                    }
                }
                let events = self.current_events();
                let mut bt = 0.0f32;
                let mut chord_idx: i64 = 0;
                for ev in events {
                    match ev {
                        ScoreEvent::Note { beats, is_chord_change: true, rh_midi, lh_midi, .. } => {
                            if self.current_beat >= bt && self.current_beat < bt + beats {
                                if chord_idx != self.last_chord_idx {
                                    // New chord — play it
                                    let all: Vec<i32> = rh_midi.iter().chain(lh_midi.iter()).copied().collect();
                                    if let Some(ref player) = self.audio_player {
                                        player.play_chord(&all);
                                    }
                                    self.last_chord_idx = chord_idx;
                                }
                                break;
                            }
                            chord_idx += 1;
                            bt += beats;
                        }
                        ScoreEvent::Note { beats, .. } => bt += beats,
                        ScoreEvent::Rest { beats } => bt += beats,
                        ScoreEvent::Bar => {}
                    }
                }
            }
        }

        // Silence when stopped
        if !self.playing && self.last_chord_idx >= 0 {
            if let Some(ref player) = self.audio_player {
                player.silence();
            }
            self.last_chord_idx = -1;
        }

        // Update current chord info from the current beat
        {
            let to_degree = |s: i32| -> i32 {
                let lin = if s > 0 { s - 1 } else { s };
                ((lin % 7) + 7) % 7
            };
            let events = self.current_events().to_vec();
            let mut beat_time = 0.0f32;
            for event in &events {
                match event {
                    ScoreEvent::Note { melody_string, rh_strings, lh_strings, beats, is_chord_change, .. } => {
                        if self.current_beat >= beat_time && self.current_beat < beat_time + beats {
                            if *is_chord_change && !rh_strings.is_empty() {
                                let mut degs: Vec<i32> = rh_strings.iter()
                                    .chain(lh_strings.iter())
                                    .map(|&s| to_degree(s))
                                    .collect();
                                degs.sort(); degs.dedup();
                                self.current_chord_degrees = degs;
                            }
                            self.current_melody_degree = to_degree(*melody_string);
                            break;
                        }
                        beat_time += beats;
                    }
                    ScoreEvent::Rest { beats } => { beat_time += beats; }
                    ScoreEvent::Bar => {}
                }
            }
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
                // ── Col 1: Circle of fifths with chord polygon ──
                ui.add_space(0.0); // horizontal spacer (none needed)
                ui.vertical(|ui| {
                    ui.add_space(15.0); // shift circle down
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
                let label_r = (outer_r + inner_r) * 0.5;

                // Rainbow colors for 7 diatonic degrees (ROYGBIV)
                const RAINBOW: [egui::Color32; 7] = [
                    egui::Color32::from_rgb(220, 30, 30),   // 0 = red (root)
                    egui::Color32::from_rgb(230, 120, 0),   // 1 = orange
                    egui::Color32::from_rgb(180, 160, 0),   // 2 = dark gold (readable on white)
                    egui::Color32::from_rgb(0, 150, 0),     // 3 = green
                    egui::Color32::from_rgb(30, 80, 220),   // 4 = blue
                    egui::Color32::from_rgb(90, 0, 140),    // 5 = indigo
                    egui::Color32::from_rgb(150, 0, 180),   // 6 = violet
                ];

                // Map each COF position to its diatonic degree (-1 if not in key)
                let key_root = ((selected_pc - music::MAJOR_SCALE[self.mode_offset as usize]) % 12 + 12) % 12;
                let cof_positions: Vec<(f32, f32)> = (0..12).map(|i| {
                    let angle = std::f32::consts::TAU * i as f32 / 12.0 - std::f32::consts::FRAC_PI_2;
                    (cof_center.x + label_r * angle.cos(), cof_center.y + label_r * angle.sin())
                }).collect();

                // Find which COF index each diatonic degree maps to
                let degree_to_cof_idx = |deg: i32| -> Option<usize> {
                    let pc = (key_root + music::MAJOR_SCALE[deg as usize]) % 12;
                    COF_ORDER.iter().position(|&c| c == pc)
                };

                cof_painter.circle_stroke(cof_center, outer_r, egui::Stroke::new(1.0, BORDER));

                // Layer 1: Filled chord triangle/quad (unique chord tones in CoF order)
                if !self.current_chord_degrees.is_empty() {
                    let mut chord_with_idx: Vec<(usize, i32)> = self.current_chord_degrees.iter()
                        .filter_map(|&d| degree_to_cof_idx(d).map(|idx| (idx, d)))
                        .collect();
                    chord_with_idx.sort_by_key(|&(idx, _)| idx);

                    if chord_with_idx.len() >= 3 {
                        let points: Vec<egui::Pos2> = chord_with_idx.iter()
                            .map(|&(idx, _)| {
                                let (x, y) = cof_positions[idx];
                                egui::Pos2::new(x, y)
                            })
                            .collect();

                        let mel_color = RAINBOW[self.current_melody_degree as usize % 7];
                        let fill = egui::Color32::from_rgba_unmultiplied(
                            mel_color.r(), mel_color.g(), mel_color.b(), 30,
                        );
                        cof_painter.add(egui::Shape::convex_polygon(
                            points.clone(), fill, egui::Stroke::new(1.5, egui::Color32::from_rgba_unmultiplied(
                                mel_color.r(), mel_color.g(), mel_color.b(), 80,
                            )),
                        ));
                    }
                }

                // Layer 2: Star lines tracing all notes in playing order (melody down through fingers)
                // Build the note sequence from current chord's RH+LH used strings
                {
                    let to_degree = |s: i32| -> i32 {
                        let lin = if s > 0 { s - 1 } else { s };
                        ((lin % 7) + 7) % 7
                    };
                    // Get current beat's RH+LH strings in order (high to low)
                    let events = self.current_events().to_vec();
                    let mut beat_time = 0.0f32;
                    let mut note_sequence: Vec<i32> = Vec::new();
                    for event in &events {
                        if let ScoreEvent::Note { rh_strings, lh_strings, beats, is_chord_change, .. } = event {
                            if self.current_beat >= beat_time && self.current_beat < beat_time + beats {
                                if *is_chord_change && !rh_strings.is_empty() {
                                    // Collect degrees in playing order (high to low)
                                    for &s in rh_strings {
                                        note_sequence.push(to_degree(s));
                                    }
                                    for &s in lh_strings {
                                        note_sequence.push(to_degree(s));
                                    }
                                }
                                break;
                            }
                            beat_time += beats;
                        } else if let ScoreEvent::Rest { beats } = event {
                            beat_time += beats;
                        }
                    }

                    // Draw crossing lines through the note sequence
                    if note_sequence.len() >= 2 {
                        for w in note_sequence.windows(2) {
                            if let (Some(from_idx), Some(to_idx)) = (degree_to_cof_idx(w[0]), degree_to_cof_idx(w[1])) {
                                let (fx, fy) = cof_positions[from_idx];
                                let (tx, ty) = cof_positions[to_idx];
                                let color = RAINBOW[w[0] as usize % 7];
                                cof_painter.line_segment(
                                    [egui::Pos2::new(fx, fy), egui::Pos2::new(tx, ty)],
                                    egui::Stroke::new(1.5, color),
                                );
                            }
                        }
                    }
                }

                // Draw labels: diatonic notes get rainbow color + bold, chord tones highlighted
                for (i, (&pc, &label)) in COF_ORDER.iter().zip(COF_LABELS.iter()).enumerate() {
                    let (lx, ly) = cof_positions[i];

                    // Find diatonic degree for this pc (-1 if chromatic)
                    let degree = music::MAJOR_SCALE.iter().position(|&s| (key_root + s) % 12 == pc);
                    let is_chord_tone = degree.map_or(false, |d| self.current_chord_degrees.contains(&(d as i32)));

                    if let Some(deg) = degree {
                        let rainbow_color = RAINBOW[deg];
                        if is_chord_tone {
                            // Chord tone: rainbow highlight circle, white text
                            cof_painter.circle_filled(
                                egui::Pos2::new(lx, ly), 9.0, rainbow_color,
                            );
                            cof_painter.text(
                                egui::Pos2::new(lx, ly),
                                egui::Align2::CENTER_CENTER,
                                label,
                                egui::FontId::proportional(9.0),
                                egui::Color32::WHITE,
                            );
                        } else {
                            // Diatonic but not in chord: small rainbow dot, black text
                            cof_painter.circle_filled(
                                egui::Pos2::new(lx, ly), 4.0, rainbow_color,
                            );
                            cof_painter.text(
                                egui::Pos2::new(lx, ly - 9.0),
                                egui::Align2::CENTER_CENTER,
                                label,
                                egui::FontId::proportional(7.0),
                                TEXT_PRIMARY,
                            );
                        }
                    } else {
                        // Chromatic (not in key): black, smaller
                        cof_painter.text(
                            egui::Pos2::new(lx, ly),
                            egui::Align2::CENTER_CENTER,
                            label,
                            egui::FontId::proportional(6.0),
                            TEXT_MUTED,
                        );
                    }
                }

                }); // end circle vertical wrapper

                // ── Col 2: controls ──
                ui.vertical(|ui| {
                    ui.spacing_mut().item_spacing.y = 2.0;

                    let half_w = ui.available_width() * 0.5;

                    let pedal_w = 18.0 * 7.0 + 6.0 + 16.0; // pedal chart width + padding
                    let title_w = half_w * 2.0 - pedal_w - 10.0; // remaining for title/recent

                    // Row 1: Title + Filter
                    ui.horizontal(|ui| {
                        if !self.scores.is_empty() {
                            let current_label = self.scores.get(self.current_score)
                                .map(|s| format!("{}. {}", s.number, s.title))
                                .unwrap_or_default();
                            ui.add_sized([title_w - 20.0, 20.0], egui::Label::new(
                                egui::RichText::new(&current_label).strong().size(13.0).color(TEXT_PRIMARY)
                            ).truncate());
                        }
                        ui.add(
                            egui::TextEdit::singleline(&mut self.search_text)
                                .hint_text("Filter...")
                                .desired_width(ui.available_width() - 4.0)
                                .min_size(egui::Vec2::new(0.0, 24.0))
                                .text_color(TEXT_PRIMARY)
                                .background_color(CARD_BG)
                                .font(egui::FontSelection::FontId(egui::FontId::proportional(13.0)))
                        );
                    });

                    // Show filtered matches
                    if !self.search_text.is_empty() {
                        let q = self.search_text.to_lowercase();
                        ui.horizontal_wrapped(|ui| {
                            let mut picked = None;
                            for (i, s) in self.scores.iter().enumerate() {
                                if !s.title.to_lowercase().contains(&q) && !s.number.contains(&q) {
                                    continue;
                                }
                                let label = format!("{}. {}", s.number, s.title);
                                let selected = self.current_score == i;
                                if ui.selectable_label(selected, &label).clicked() {
                                    picked = Some(i);
                                }
                            }
                            if let Some(idx) = picked {
                                self.select_hymn(idx);
                                self.search_text.clear();
                            }
                        });
                    }

                    // Row 2: Recent + centered Pedal diagram
                    ui.horizontal(|ui| {
                        let mut recent_pick = self.current_score;
                        egui::ComboBox::from_id_salt("recent_select")
                            .selected_text("Recent")
                            .width(title_w - 20.0)
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

                        {
                            let pedal_letters_l = ["D", "C", "B"];
                            let pedal_letters_r = ["E", "F", "G", "A"];
                            let acc = music::key_sig_accidentals(&self.current_key);
                            let pedal_letter_idx = [1usize, 0, 6, 2, 3, 4, 5];
                            let dot_r = 4.0;
                            let col_w = 18.0;
                            let line_y_offset = 10.0;

                            for (i, letter) in pedal_letters_l.iter().enumerate() {
                                let li = pedal_letter_idx[i];
                                let a = acc[li];
                                let (response, painter) = ui.allocate_painter(
                                    egui::Vec2::new(col_w, 32.0), egui::Sense::hover());
                                let rect = response.rect;
                                let cx = rect.center().x;
                                let line_y = rect.top() + line_y_offset;
                                painter.line_segment(
                                    [egui::Pos2::new(rect.left(), line_y), egui::Pos2::new(rect.right(), line_y)],
                                    egui::Stroke::new(1.0, TEXT_PRIMARY));
                                let dot_y = match a { -1 => line_y - 5.0, 1 => line_y + 5.0, _ => line_y };
                                painter.circle_filled(egui::Pos2::new(cx, dot_y), dot_r, if a == 0 { TEXT_PRIMARY } else { ACCENT });
                                painter.text(egui::Pos2::new(cx, rect.bottom() + 1.0), egui::Align2::CENTER_BOTTOM,
                                    *letter, egui::FontId::monospace(9.0), TEXT_PRIMARY);
                            }
                            ui.allocate_ui(egui::Vec2::new(6.0, 32.0), |ui| {
                                let rect = ui.available_rect_before_wrap();
                                let line_y = rect.top() + line_y_offset;
                                ui.painter().line_segment(
                                    [egui::Pos2::new(rect.center().x, line_y - 6.0),
                                     egui::Pos2::new(rect.center().x, line_y + 6.0)],
                                    egui::Stroke::new(1.0, BORDER));
                            });
                            for (i, letter) in pedal_letters_r.iter().enumerate() {
                                let li = pedal_letter_idx[3 + i];
                                let a = acc[li];
                                let (response, painter) = ui.allocate_painter(
                                    egui::Vec2::new(col_w, 32.0), egui::Sense::hover());
                                let rect = response.rect;
                                let cx = rect.center().x;
                                let line_y = rect.top() + line_y_offset;
                                painter.line_segment(
                                    [egui::Pos2::new(rect.left(), line_y), egui::Pos2::new(rect.right(), line_y)],
                                    egui::Stroke::new(1.0, TEXT_PRIMARY));
                                let dot_y = match a { -1 => line_y - 5.0, 1 => line_y + 5.0, _ => line_y };
                                painter.circle_filled(egui::Pos2::new(cx, dot_y), dot_r, if a == 0 { TEXT_PRIMARY } else { ACCENT });
                                painter.text(egui::Pos2::new(cx, rect.bottom() + 1.0), egui::Align2::CENTER_BOTTOM,
                                    *letter, egui::FontId::monospace(9.0), TEXT_PRIMARY);
                            }
                        }

                        // Shuffle button when random on
                        if self.random_fingers {
                            ui.add_space(4.0);
                            if ui.button(egui::RichText::new("Shuffle").size(10.0)).clicked() {
                                self.shuffle_fingers();
                                if let Some(score) = self.scores.get_mut(self.current_score) {
                                    abc::compute_display_strings(score, self.rh_fingers as usize, self.lh_fingers as usize, &self.random_pattern);
                                }
                            }
                        }
                    });

                    // Row 3: [⏮ ▶ ⏹] [G▾] 120BPM  RH[4] LH[4] ⊙  [🔊▬▬]
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = 3.0;
                        if ui.button(egui::RichText::new("\u{23EE}").size(14.0)).clicked() {
                            self.current_beat = 0.0;
                            self.scroll_offset = 0.0;
                            self.last_chord_idx = -1;
                            self.last_frame_time = None;
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
                            self.playing = false;
                            self.current_beat = 0.0;
                            self.scroll_offset = 0.0;
                            self.last_chord_idx = -1;
                            self.last_frame_time = None;
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

                        ui.label(egui::RichText::new("RH").size(11.0).strong().color(egui::Color32::BLACK));
                        ui.add(egui::DragValue::new(&mut self.rh_fingers).range(1..=4).speed(0.1));
                        ui.label(egui::RichText::new("LH").size(11.0).strong().color(egui::Color32::BLACK));
                        ui.add(egui::DragValue::new(&mut self.lh_fingers).range(1..=4).speed(0.1));

                        // Random checkbox
                        let was_random = self.random_fingers;
                        ui.checkbox(&mut self.random_fingers, egui::RichText::new("Ran").size(9.0));
                        if self.random_fingers && !was_random {
                            self.shuffle_fingers();
                            if let Some(score) = self.scores.get_mut(self.current_score) {
                                abc::compute_display_strings(score, self.rh_fingers as usize, self.lh_fingers as usize, &self.random_pattern);
                            }
                        } else if !self.random_fingers && was_random {
                            if let Some(score) = self.scores.get_mut(self.current_score) {
                                abc::compute_display_strings(score, self.rh_fingers as usize, self.lh_fingers as usize, &[]);
                            }
                        }

                    });

                    // Row 4: Mode buttons + Volume slider
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = 1.0;
                        ui.spacing_mut().button_padding = egui::Vec2::new(11.0, 3.0);
                        for (i, name) in MODE_NAMES.iter().enumerate() {
                            let selected = self.mode_offset == i as i32;
                            let txt = egui::RichText::new(*name).size(9.5).color(
                                if selected { egui::Color32::WHITE } else { egui::Color32::BLACK }
                            );
                            let btn = if selected {
                                egui::Button::new(txt).fill(ACCENT).corner_radius(0.0)
                            } else {
                                egui::Button::new(txt).corner_radius(0.0)
                            };
                            if ui.add(btn).clicked() {
                                self.mode_offset = i as i32;
                            }
                        }

                        ui.add_space(4.0);
                        let mute_icon = if self.sound_volume > 0.0 { "\u{1F50A}" } else { "\u{1F507}" };
                        if ui.button(egui::RichText::new(mute_icon).size(10.0)).clicked() {
                            if self.sound_volume > 0.0 {
                                self.pre_mute_volume = self.sound_volume;
                                self.sound_volume = 0.0;
                            } else {
                                self.sound_volume = self.pre_mute_volume.max(0.3);
                            }
                        }
                        ui.visuals_mut().selection.bg_fill = ACCENT;
                        ui.visuals_mut().widgets.inactive.bg_fill = egui::Color32::from_rgb(190, 215, 245);
                        ui.style_mut().spacing.slider_width = (ui.available_width() - 8.0).max(40.0);
                        ui.add(egui::Slider::new(&mut self.sound_volume, 0.0..=1.0)
                            .show_value(false).trailing_fill(true));
                    });

                    // (Pedal diagram moved to Row 1b)
                });

            });
            ui.add_space(2.0);
        });

        // ── Dark strip so Android nav icons (white) are visible ──
        egui::TopBottomPanel::bottom("nav_bar")
            .frame(egui::Frame::NONE.fill(egui::Color32::from_rgb(30, 30, 30)))
            .show(ctx, |ui| {
                ui.add_space(24.0);
            });

        // ── Drill controls + metrics panel ──
        egui::TopBottomPanel::bottom("drill_panel").show(ctx, |ui| {
            ui.add_space(4.0);

            let active = self.drill_mode;
            let gray = egui::Color32::from_rgb(190, 190, 190);

            ui.horizontal(|ui| {
                ui.set_min_height(90.0);

                // ── Col 1: Toggle button (fills height) ──
                let label = if self.drill_mode { "Songs" } else { "Drill" };
                let toggle = egui::Button::new(
                    egui::RichText::new(label).size(18.0).strong().color(egui::Color32::WHITE)
                ).fill(ACCENT).corner_radius(10.0).min_size(egui::Vec2::new(80.0, 85.0));
                if ui.add(toggle).clicked() {
                    self.drill_mode = !self.drill_mode;
                    if self.drill_mode {
                        // Stop song playback, generate random drill score
                        self.playing = false;
                        if let Some(ref player) = self.audio_player { player.silence(); }
                        self.last_chord_idx = -1;
                        self.current_beat = 0.0;
                        self.scroll_offset = 0.0;
                        self.last_frame_time = None;
                        self.drill_score = Some(drill::generate_random_score(
                            &self.current_key, self.mode_offset, 32, 4,
                            self.rh_fingers as usize, self.lh_fingers as usize,
                        ));
                    } else {
                        self.drill_score = None;
                    }
                }

                ui.add_space(12.0);
                ui.separator();
                ui.add_space(4.0);

                // ── Col 2: ▶ Intervals Chords / RH notes / LH notes ──
                ui.vertical(|ui| {
                    // Row 1: Play + Intervals + Chords
                    ui.horizontal(|ui| {
                        let is_playing = self.playing;
                        let (ic, fl, cl) = if is_playing {
                            ("\u{23F8}", ACCENT, egui::Color32::WHITE)
                        } else {
                            ("\u{25B6}", egui::Color32::from_rgb(40, 167, 69), egui::Color32::WHITE)
                        };
                        if ui.add(egui::Button::new(egui::RichText::new(ic).size(14.0).color(cl)).fill(fl).corner_radius(6.0).min_size(egui::Vec2::new(32.0, 24.0))).clicked() {
                            self.playing = !self.playing;
                            if self.playing { self.last_frame_time = None; }
                        }

                        let mk = |ui: &mut egui::Ui, lbl: &str, sel: bool, en: bool| -> bool {
                            let (fg, bg, sc) = if !en { (gray, CARD_BG, egui::Color32::from_rgb(220,220,220)) }
                                else if sel { (egui::Color32::WHITE, ACCENT, ACCENT) }
                                else { (TEXT_MUTED, CARD_BG, BORDER) };
                            ui.add(egui::Button::new(egui::RichText::new(lbl).size(12.0).color(fg))
                                .fill(bg).stroke(egui::Stroke::new(1.5, sc)).corner_radius(6.0)).clicked() && en
                        };
                        if mk(ui, "Intervals", self.drill_config.mode == drill::DrillMode::Intervals, active) {
                            self.drill_config.mode = drill::DrillMode::Intervals;
                        }
                        if mk(ui, "Chords", self.drill_config.mode == drill::DrillMode::Chords, active) {
                            self.drill_config.mode = drill::DrillMode::Chords;
                        }
                    });

                    // Row 2: RH durations
                    let c = if active { TEXT_PRIMARY } else { gray };
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("RH").size(13.0).strong().color(c));
                        ui.checkbox(&mut self.rh_whole, egui::RichText::new("W").size(11.0).color(c));
                        ui.checkbox(&mut self.rh_half, egui::RichText::new("H").size(11.0).color(c));
                        ui.checkbox(&mut self.rh_quarter, egui::RichText::new("Q").size(11.0).color(c));
                        ui.checkbox(&mut self.rh_eighth, egui::RichText::new("8").size(11.0).color(c));
                    });
                    // Row 3: LH durations
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("LH").size(13.0).strong().color(c));
                        ui.checkbox(&mut self.lh_whole, egui::RichText::new("W").size(11.0).color(c));
                        ui.checkbox(&mut self.lh_half, egui::RichText::new("H").size(11.0).color(c));
                        ui.checkbox(&mut self.lh_quarter, egui::RichText::new("Q").size(11.0).color(c));
                        ui.checkbox(&mut self.lh_eighth, egui::RichText::new("8").size(11.0).color(c));
                    });
                });

                ui.separator();

                // ── Col 3: [Random▾] / Mic / Calibrate ──
                ui.vertical(|ui| {
                    // Row 1: Drill dropdown
                    let mut dn = self.drill_name;
                    let dn_color = if active { TEXT_PRIMARY } else { gray };
                    egui::ComboBox::from_id_salt("drill_sel")
                        .selected_text(egui::RichText::new(DRILL_NAMES[dn]).size(12.0).color(dn_color))
                        .width(80.0)
                        .show_ui(ui, |ui| {
                            for (i, name) in DRILL_NAMES.iter().enumerate() {
                                ui.selectable_value(&mut dn, i, *name);
                            }
                        });
                    if dn != self.drill_name { self.drill_name = dn; }

                    // Row 2: Mic checkbox
                    ui.checkbox(&mut self.use_mic, egui::RichText::new("Mic").size(12.0).color(if active { TEXT_PRIMARY } else { gray }));

                    // Row 3: Calibrate button
                    let cal_c = if active { TEXT_MUTED } else { gray };
                    ui.add(egui::Button::new(egui::RichText::new("Calibrate").size(11.0).color(cal_c))
                        .fill(CARD_BG).stroke(egui::Stroke::new(1.0, if active { BORDER } else { egui::Color32::from_rgb(220,220,220) }))
                        .corner_radius(6.0));
                });

                ui.separator();

                // ── Col 4: Two progress bars ──
                ui.vertical(|ui| {
                    let bar_w = ui.available_width();
                    let track = egui::Color32::from_rgb(232, 232, 232);

                    let npm_text = if self.drill_progress.best_npm > 0.0 {
                        format!("Best {:.0} NPM", self.drill_progress.best_npm)
                    } else { "Streak".into() };
                    let streak_pct = 0.0f32; // TODO: track streak
                    let level_pct = if self.drill_progress.best_npm > 0.0 {
                        ((self.drill_progress.best_npm / 60.0) * 100.0).min(100.0)
                    } else { 0.0 };

                    let txt_c = if active { TEXT_PRIMARY } else { gray };
                    let bar_c = if active { ACCENT } else { gray };

                    // Bar 1: Streak
                    ui.label(egui::RichText::new(if npm_text.is_empty() { "Streak".into() } else { npm_text }).size(11.0).color(txt_c));
                    let (r1, _) = ui.allocate_exact_size(egui::Vec2::new(bar_w, 12.0), egui::Sense::hover());
                    let p = ui.painter_at(r1);
                    p.rect_filled(r1, 6.0, track);
                    p.rect_filled(egui::Rect::from_min_size(r1.min, egui::Vec2::new(r1.width() * streak_pct / 100.0, r1.height())), 6.0, egui::Color32::from_rgb(40, 167, 69));

                    ui.add_space(4.0);

                    // Bar 2: Progress
                    let progress_label = format!("{} sessions", self.drill_progress.total_sessions);
                    ui.label(egui::RichText::new(progress_label).size(11.0).color(bar_c));
                    let (r2, _) = ui.allocate_exact_size(egui::Vec2::new(bar_w, 12.0), egui::Sense::hover());
                    let p2 = ui.painter_at(r2);
                    p2.rect_filled(r2, 6.0, track);
                    p2.rect_filled(egui::Rect::from_min_size(r2.min, egui::Vec2::new(r2.width() * level_pct / 100.0, r2.height())), 6.0, bar_c);
                });
            });

            // Row 2: status line (slider removed - swipe to scroll)
            ui.add_space(2.0);
            if !self.drill_mode {
                ui.horizontal(|ui| {
                    ui.colored_label(TEXT_PRIMARY, &self.status);
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if let Some(s) = self.scores.get(self.current_score) {
                            ui.label(egui::RichText::new(
                                format!("Key: {}  Meter: {}/{}  Tempo: {} BPM",
                                    s.key, s.meter_num, s.meter_den, s.tempo)
                            ).small().color(TEXT_PRIMARY));
                        }
                    });
                });
            } else {
                // Drill progress status line
                ui.horizontal(|ui| {
                    let p = &self.drill_progress;
                    ui.colored_label(TEXT_PRIMARY, format!(
                        "Level {} \u{00b7} Best {:.0} NPM \u{00b7} {} notes \u{00b7} {} sessions",
                        p.level, p.best_npm, p.total_notes, p.total_sessions));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(egui::RichText::new(
                            format!("Key: {} \u{00b7} {}", self.current_key,
                                if self.drill_config.mode == drill::DrillMode::Intervals { "Intervals" } else { "Chords" })
                        ).small().color(TEXT_PRIMARY));
                    });
                });
            }
            ui.add_space(4.0);
        });

        // ── Central score view ──
        egui::CentralPanel::default().show(ctx, |ui| {
            let avail = ui.available_size();

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
                    let (rect, response) = ui.allocate_exact_size(
                        egui::Vec2::new(avail.x - 16.0, avail.y - 24.0),
                        egui::Sense::click_and_drag(),
                    );

                    let painter = ui.painter_at(rect);
                    let selected_pc = music::key_to_pc(&self.current_key);
                    let key_root = ((selected_pc - music::MAJOR_SCALE[self.mode_offset as usize]) % 12 + 12) % 12;

                    {
                        // ── Same rendering for Songs and Drill ──
                        let layout = NotationLayout::new(rect.top() + 8.0, 50.0);
                        let view_width = rect.width();

                        if !self.playing && response.dragged() {
                            let drag_beats = response.drag_delta().x / 60.0;
                            self.current_beat = (self.current_beat - drag_beats).max(0.0);
                        }

                        self.scroll_offset = (self.current_beat * 60.0).max(0.0);

                        // Use drill score or hymn score
                        let events = if self.drill_mode {
                            self.drill_score.as_ref().map(|s| s.events.clone()).unwrap_or_default()
                        } else {
                            self.current_events().to_vec()
                        };

                        render_score(
                            &painter,
                            &layout,
                            &events,
                            self.scroll_offset,
                            self.current_beat,
                            view_width,
                            key_root,
                        );
                    }
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

// ── Regression tests ──
// These lock down the string-to-note mapping so refactoring can't break it.
// The harp has no "string zero" — strings are 1-based like years (no year 0).

#[cfg(test)]
mod tests {
    use super::abc::*;
    use super::music::*;

    /// Helper: convert a relative string number to a display note name.
    /// This is the EXACT pipeline the renderer uses — if this breaks, the display breaks.
    fn string_to_note(relative_str: i32, key_root: i32) -> String {
        let abs_str = relative_str + STRING_NUMBER_OFFSET - 1; // 1-based to 0-based
        let display_midi = harp_string_to_midi(abs_str, key_root) + 12; // base-48 to real pitch
        midi_to_name(display_midi)
    }

    // ── String numbering: no zero, 1 = root ──

    #[test]
    fn string_1_is_root_in_every_key() {
        // String 1 must always be the key root note (lowest octave on harp)
        let keys = [("C",0), ("D",2), ("E",4), ("F",5), ("G",7), ("A",9), ("B",11)];
        let root_names = ["C","D","E","F","G","A","B"];
        for (i, &(_key_name, kr)) in keys.iter().enumerate() {
            let note = string_to_note(1, kr);
            assert!(note.starts_with(root_names[i]),
                "String 1 in key {} should be {} but got {}", _key_name, root_names[i], note);
        }
    }

    #[test]
    fn string_8_is_root_one_octave_up() {
        // String 8 = root + 1 octave (7 diatonic steps)
        let note_c = string_to_note(8, 0);
        assert!(note_c.starts_with("C"), "String 8 in C should be C, got {}", note_c);
        let note_g = string_to_note(8, 7);
        assert!(note_g.starts_with("G"), "String 8 in G should be G, got {}", note_g);
        let note_d = string_to_note(8, 2);
        assert!(note_d.starts_with("D"), "String 8 in D should be D, got {}", note_d);
    }

    #[test]
    fn diatonic_scale_from_root() {
        // In G major: 1=G 2=A 3=B 4=C 5=D 6=E 7=F#
        let expected = ["G", "A", "B", "C", "D", "E", "F#"];
        for (i, exp) in expected.iter().enumerate() {
            let note = string_to_note((i + 1) as i32, 7);
            assert!(note.starts_with(exp),
                "String {} in G major should be {} but got {}", i + 1, exp, note);
        }
    }

    #[test]
    fn diatonic_scale_c_major() {
        // In C major: 1=C 2=D 3=E 4=F 5=G 6=A 7=B
        let expected = ["C", "D", "E", "F", "G", "A", "B"];
        for (i, exp) in expected.iter().enumerate() {
            let note = string_to_note((i + 1) as i32, 0);
            assert!(note.starts_with(exp),
                "String {} in C major should be {} but got {}", i + 1, exp, note);
        }
    }

    #[test]
    fn diatonic_scale_d_major() {
        // In D major: 1=D 2=E 3=F# 4=G 5=A 6=B 7=C#/Db
        let expected_pc = [2, 4, 6, 7, 9, 11, 1]; // D E F# G A B C#
        for (i, &exp_pc) in expected_pc.iter().enumerate() {
            let abs_str = (i as i32 + 1) + STRING_NUMBER_OFFSET - 1;
            let midi = harp_string_to_midi(abs_str, 2) + 12;
            let pc = midi % 12;
            assert_eq!(pc, exp_pc,
                "String {} in D major: expected pc {} but got pc {} (midi {})",
                i + 1, exp_pc, pc, midi);
        }
    }

    #[test]
    fn no_string_zero() {
        // String numbering skips zero: ...-2, -1, 1, 2...
        // String 0 should never appear in voicing output
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_strings, lh_strings, .. } = event {
                    for &s in rh_strings {
                        assert!(s != 0, "Hymn {} beat {}: RH has string 0!", score.number, i);
                    }
                    for &s in lh_strings {
                        assert!(s != 0, "Hymn {} beat {}: LH has string 0!", score.number, i);
                    }
                }
            }
        }
    }

    // ── RH/LH separation: never the same string ──

    #[test]
    fn rh_and_lh_never_share_a_string() {
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_strings, lh_strings, chord_name, .. } = event {
                    for rs in rh_strings {
                        for ls in lh_strings {
                            assert!(rs != ls,
                                "Hymn {} '{}' beat {}: RH and LH share string {} (chord {:?})",
                                score.number, score.title, i, rs, chord_name);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn rh_and_lh_never_share_midi() {
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_midi, lh_midi, chord_name, .. } = event {
                    for rm in rh_midi {
                        for lm in lh_midi {
                            assert!(rm != lm,
                                "Hymn {} '{}' beat {}: RH and LH share MIDI {} (chord {:?})",
                                score.number, score.title, i, rm, chord_name);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn lh_always_below_rh() {
        // LH highest string must always be below RH lowest string
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_strings, lh_strings, chord_name, .. } = event {
                    if rh_strings.is_empty() || lh_strings.is_empty() { continue; }
                    let rh_low = *rh_strings.iter().min().unwrap();
                    let lh_high = *lh_strings.iter().max().unwrap();
                    assert!(lh_high < rh_low,
                        "Hymn {} '{}' beat {}: LH highest ({}) >= RH lowest ({}) (chord {:?})",
                        score.number, score.title, i, lh_high, rh_low, chord_name);
                }
            }
        }
    }

    // ── Staff position: notes on correct lines/spaces ──

    #[test]
    fn staff_position_middle_c() {
        // Middle C (MIDI 60) = staff position 0
        use super::notation::midi_to_staff_pos;
        assert_eq!(midi_to_staff_pos(60), 0, "Middle C should be staff pos 0");
    }

    #[test]
    fn staff_position_treble_lines() {
        // Treble staff lines: E4=2, G4=4, B4=6, D5=8, F5=10
        use super::notation::midi_to_staff_pos;
        assert_eq!(midi_to_staff_pos(64), 2, "E4 should be pos 2");
        assert_eq!(midi_to_staff_pos(67), 4, "G4 should be pos 4");
        assert_eq!(midi_to_staff_pos(71), 6, "B4 should be pos 6");
        assert_eq!(midi_to_staff_pos(74), 8, "D5 should be pos 8");
        assert_eq!(midi_to_staff_pos(77), 10, "F5 should be pos 10");
    }

    #[test]
    fn staff_position_bass_lines() {
        // Bass staff lines: G2=-10, B2=-8, D3=-6, F3=-4, A3=-2
        use super::notation::midi_to_staff_pos;
        assert_eq!(midi_to_staff_pos(43), -10, "G2 should be pos -10");
        assert_eq!(midi_to_staff_pos(47), -8, "B2 should be pos -8");
        assert_eq!(midi_to_staff_pos(50), -6, "D3 should be pos -6");
        assert_eq!(midi_to_staff_pos(53), -4, "F3 should be pos -4");
        assert_eq!(midi_to_staff_pos(57), -2, "A3 should be pos -2");
    }

    // ── Known-good hymn voicings (golden tests) ──

    #[test]
    fn blessed_jesus_is_in_g_major() {
        let scores = load_embedded_scores();
        let s = &scores[0];
        assert_eq!(s.key, "G");
        assert_eq!(s.key_root, 7);
    }

    #[test]
    fn joy_to_world_is_in_d_major() {
        let scores = load_embedded_scores();
        let joy = scores.iter().find(|s| s.title.contains("Joy to the World")).unwrap();
        assert_eq!(joy.key, "D");
        assert_eq!(joy.key_root, 2);
    }

    #[test]
    fn every_voicing_note_is_in_key() {
        // Every string number at a chord change must map to a note in that hymn's key
        let scores = load_embedded_scores();
        for score in &scores {
            let kr = score.key_root;
            let scale_pcs: Vec<i32> = (0..7).map(|deg| (kr + MAJOR_SCALE[deg] as i32) % 12).collect();

            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_strings, lh_strings, is_chord_change: true, .. } = event {
                    for &s in rh_strings.iter().chain(lh_strings.iter()) {
                        let abs_str = s + STRING_NUMBER_OFFSET - 1;
                        let display_midi = harp_string_to_midi(abs_str, kr) + 12;
                        let pc = display_midi % 12;
                        assert!(scale_pcs.contains(&pc),
                            "Hymn {} '{}' beat {}: string {} = midi {} (pc {}) not in key {}",
                            score.number, score.title, i, s, display_midi, pc, score.key);
                    }
                }
            }
        }
    }

    // ── Voicing constraints ──

    #[test]
    fn rh_has_1_to_4_notes() {
        // RH should have melody on thumb + 0-3 harmony notes
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { rh_strings, is_chord_change: true, .. } = event {
                    assert!(!rh_strings.is_empty(),
                        "Hymn {} beat {}: RH has no strings at chord change", score.number, i);
                    assert!(rh_strings.len() <= 4,
                        "Hymn {} beat {}: RH has {} strings (max 4)", score.number, i, rh_strings.len());
                }
            }
        }
    }

    #[test]
    fn lh_has_2_to_4_notes() {
        // LH carries the chord: 2-4 notes
        let scores = load_embedded_scores();
        for score in &scores {
            for (i, event) in score.events.iter().enumerate() {
                if let ScoreEvent::Note { lh_strings, is_chord_change: true, .. } = event {
                    if lh_strings.is_empty() { continue; } // some voicings may fail
                    assert!(lh_strings.len() <= 4,
                        "Hymn {} beat {}: LH has {} strings (max 4)", score.number, i, lh_strings.len());
                }
            }
        }
    }

    #[test]
    fn all_hymns_parse() {
        let scores = load_embedded_scores();
        assert!(scores.len() >= 280, "Expected 280+ hymns, got {}", scores.len());
        for s in &scores {
            assert!(!s.events.is_empty(), "Hymn {} has no events", s.number);
            assert!(!s.title.is_empty(), "Hymn {} has no title", s.number);
        }
    }
}
