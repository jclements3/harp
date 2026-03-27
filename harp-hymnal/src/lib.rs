pub mod abc;
pub mod chord;
pub mod music;
pub mod voicing;

use eframe::egui;
use abc::{Beat, Hymn, parse_abc};
use music::*;

// ── Embedded ABC data (for Android and fallback) ──
const EMBEDDED_ABC_BYTES: &[u8] = include_bytes!("../../data/OpenHymnal.abc");

// ── Colors matching the HTML version ──
const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241);       // #f8f6f1
const CARD_BG: egui::Color32 = egui::Color32::WHITE;
const TEXT_PRIMARY: egui::Color32 = egui::Color32::from_rgb(42, 42, 42); // #2a2a2a
const TEXT_SECONDARY: egui::Color32 = egui::Color32::from_rgb(74, 74, 74); // #4a4a4a
const TEXT_MUTED: egui::Color32 = egui::Color32::from_rgb(153, 153, 153); // #999
const CHORD_LABEL: egui::Color32 = egui::Color32::from_rgb(170, 170, 170); // #aaa
const ACCENT: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);     // #2563eb
const RH_FINGER: egui::Color32 = egui::Color32::from_rgb(51, 51, 51);   // #333
const LH_FINGER: egui::Color32 = egui::Color32::from_rgb(102, 102, 102); // #666
const BORDER: egui::Color32 = egui::Color32::from_rgb(204, 204, 204);   // #ccc
const BAR_COLOR: egui::Color32 = egui::Color32::from_rgb(51, 51, 51);   // #333
const KEY_BADGE_BG: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);

const ALL_KEYS: [&str; 12] = ["C","Db","D","Eb","E","F","F#","G","Ab","A","Bb","B"];
const MODE_NAMES: [&str; 7] = ["Ionian","Dorian","Phrygian","Lydian","Mixolydian","Aeolian","Locrian"];

// Font sizes matching the HTML: cells 0.85rem (~12px), chord row 0.75rem (~10.5px), labels 0.65rem (~9px)
const CELL_SIZE: f32 = 12.0;
const CHORD_SIZE: f32 = 11.0;
const LABEL_SIZE: f32 = 9.5;
const BAR_W: f32 = 10.0;   // barline column width
const LABEL_W: f32 = 32.0; // row label column width
const ROW_H: f32 = 15.0;   // row height

pub fn apply_style(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing = egui::Vec2::new(4.0, 2.0);
    style.spacing.button_padding = egui::Vec2::new(8.0, 4.0);
    style.visuals.window_fill = BG;
    style.visuals.panel_fill = BG;
    style.visuals.extreme_bg_color = CARD_BG;

    // Rounded corners
    let r8 = egui::CornerRadius::same(8);
    style.visuals.widgets.noninteractive.corner_radius = r8;
    style.visuals.widgets.inactive.corner_radius = r8;
    style.visuals.widgets.hovered.corner_radius = r8;
    style.visuals.widgets.active.corner_radius = r8;

    // Inactive buttons: white bg, gray border, dark text
    style.visuals.widgets.inactive.bg_fill = CARD_BG;
    style.visuals.widgets.inactive.weak_bg_fill = CARD_BG;
    style.visuals.widgets.inactive.bg_stroke = egui::Stroke::new(2.0, BORDER);
    style.visuals.widgets.inactive.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);

    // Hovered: light gray bg
    style.visuals.widgets.hovered.bg_fill = egui::Color32::from_rgb(240, 240, 240);
    style.visuals.widgets.hovered.weak_bg_fill = egui::Color32::from_rgb(240, 240, 240);
    style.visuals.widgets.hovered.bg_stroke = egui::Stroke::new(2.0, egui::Color32::from_rgb(153, 153, 153));
    style.visuals.widgets.hovered.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);

    // Active/pressed: light blue
    style.visuals.widgets.active.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.widgets.active.weak_bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.widgets.active.bg_stroke = egui::Stroke::new(2.0, ACCENT);
    style.visuals.widgets.active.fg_stroke = egui::Stroke::new(1.0, ACCENT);

    // Non-interactive (labels, panels)
    style.visuals.widgets.noninteractive.bg_fill = CARD_BG;
    style.visuals.widgets.noninteractive.fg_stroke = egui::Stroke::new(1.0, TEXT_PRIMARY);
    style.visuals.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(238, 238, 238));

    // Selection
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(219, 234, 254);
    style.visuals.selection.stroke = egui::Stroke::new(2.0, ACCENT);

    // Separator
    style.visuals.widgets.noninteractive.bg_stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(238, 238, 238));

    ctx.set_style(style);
}

pub fn create_native_options() -> eframe::NativeOptions {
    eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1000.0, 700.0]),
        ..Default::default()
    }
}

pub fn run_app(options: eframe::NativeOptions) -> eframe::Result {
    eframe::run_native(
        "Harp Hymnal",
        options,
        Box::new(|cc| {
            apply_style(&cc.egui_ctx);
            Ok(Box::new(HarpHymnal::new()))
        }),
    )
}

pub struct HarpHymnal {
    hymns: Vec<Hymn>,
    current: usize,
    mode_offset: i32,
    current_key: String,
    search_text: String,
    status: String,
}

impl HarpHymnal {
    pub fn new() -> Self {
        let mut app = Self {
            hymns: Vec::new(),
            current: 0,
            mode_offset: 0,
            current_key: "C".into(),
            search_text: String::new(),
            status: String::new(),
        };

        #[allow(unused_mut)]
        let mut loaded = false;

        #[cfg(not(target_os = "android"))]
        {
            for path in &["data/OpenHymnal.abc", "data/harp_leadsheets.abc", "OpenHymnal.abc", "harp_leadsheets.abc"] {
                if let Ok(text) = std::fs::read_to_string(path) {
                    app.load_abc(&text);
                    loaded = true;
                    break;
                }
                if let Ok(exe) = std::env::current_exe() {
                    if let Some(dir) = exe.parent() {
                        let p = dir.join(path);
                        if let Ok(text) = std::fs::read_to_string(&p) {
                            app.load_abc(&text);
                            loaded = true;
                            break;
                        }
                    }
                }
            }
        }

        if !loaded {
            let text = String::from_utf8_lossy(EMBEDDED_ABC_BYTES);
            app.load_abc(&text);
        }

        app
    }

    pub fn load_abc(&mut self, text: &str) {
        self.hymns = parse_abc(text);
        self.current = 0;
        if let Some(h) = self.hymns.first() {
            self.current_key = h.key.clone();
        }
        self.status = format!("{} hymns loaded.", self.hymns.len());
    }

    fn current_hymn(&self) -> Option<&Hymn> {
        self.hymns.get(self.current)
    }

    fn key_display(&self) -> String {
        let base_pc = key_to_pc(&self.current_key);
        let modal_root = (base_pc + MAJOR_SCALE[self.mode_offset as usize]) % 12;
        format!("{} {}", PC_NAMES[modal_root as usize], MODE_NAMES[self.mode_offset as usize])
    }
}

impl eframe::App for HarpHymnal {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // ── Top controls panel ──
        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.add_space(4.0);

            // Row 1: Key badge + key selector + mode buttons (centered)
            ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                ui.horizontal(|ui| {
                    // Key badge (dark bg, white text, like HTML)
                    let key_text = self.key_display();
                    let badge = egui::Frame::NONE
                        .fill(KEY_BADGE_BG)
                        .inner_margin(egui::Margin::symmetric(6, 3))
                        .corner_radius(6.0);
                    badge.show(ui, |ui| {
                        ui.label(egui::RichText::new(&key_text)
                            .monospace().size(13.0).strong()
                            .color(egui::Color32::WHITE));
                    });

                    ui.add_space(4.0);

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

                    ui.add_space(4.0);

                    // Mode buttons
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

                    // File menu (desktop only)
                    #[cfg(not(target_os = "android"))]
                    {
                        ui.add_space(12.0);
                        if ui.button("Open ABC...").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("ABC files", &["abc"])
                                .pick_file()
                            {
                                if let Ok(text) = std::fs::read_to_string(&path) {
                                    self.load_abc(&text);
                                }
                            }
                        }
                    }
                });
            });

            // Row 2: Hymn selector + search
            ui.horizontal(|ui| {
                if !self.hymns.is_empty() {
                    let label = self.current_hymn()
                        .map(|h| format!("{}. {}", h.number, h.title))
                        .unwrap_or_default();

                    egui::ComboBox::from_id_salt("hymn_select")
                        .selected_text(&label)
                        .width(ui.available_width() * 0.55)
                        .show_ui(ui, |ui| {
                            for (i, h) in self.hymns.iter().enumerate() {
                                let text = format!("{}. {}", h.number, h.title);
                                if ui.selectable_value(&mut self.current, i, &text).clicked() {
                                    self.current_key = h.key.clone();
                                }
                            }
                        });

                    ui.add_space(8.0);

                    let response = ui.add(
                        egui::TextEdit::singleline(&mut self.search_text)
                            .hint_text("Search...")
                            .desired_width(180.0)
                    );
                    if response.changed() && !self.search_text.is_empty() {
                        let q = self.search_text.to_lowercase();
                        if let Some(idx) = self.hymns.iter().position(|h|
                            h.title.to_lowercase().contains(&q) || h.number.contains(&q)
                        ) {
                            self.current = idx;
                            self.current_key = self.hymns[idx].key.clone();
                        }
                    }
                }
            });

            ui.add_space(2.0);
        });

        // ── Bottom status ──
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.colored_label(TEXT_MUTED, &self.status);
        });

        // ── Central content ──
        egui::CentralPanel::default().show(ctx, |ui| {
            if self.hymns.is_empty() {
                ui.centered_and_justified(|ui| {
                    ui.label("Load an ABC file to begin.");
                });
                return;
            }

            let hymn = self.hymns[self.current].clone();
            let offset = self.mode_offset;
            let avail_w = ui.available_width();

            egui::ScrollArea::vertical().show(ui, |ui| {
                // Title
                ui.add_space(8.0);
                ui.label(egui::RichText::new(format!("{}. {}", hymn.number, hymn.title))
                    .strong().size(20.0).color(TEXT_SECONDARY));
                ui.add_space(8.0);

                // Leadsheet in a card
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
                        self.render_leadsheet(ui, &hymn, offset, avail_w - 24.0);
                    });

                ui.add_space(16.0);

                // Chord table in a card
                egui::Frame::NONE
                    .fill(CARD_BG)
                    .inner_margin(12.0)
                    .corner_radius(10.0)
                    .shadow(egui::epaint::Shadow {
                        offset: [0, 1],
                        blur: 6,
                        spread: 0,
                        color: egui::Color32::from_black_alpha(15),
                    })
                    .show(ui, |ui| {
                        self.render_chord_table(ui, &hymn);
                    });

                ui.add_space(16.0);
            });
        });
    }
}

// ── Leadsheet rendering with measure-fit wrapping ──

/// A single measure's worth of beat columns.
struct Measure {
    beats: Vec<BeatCol>,
}

struct BeatCol {
    is_chord_change: bool,
    melody_string: Option<i32>,
    rh: Vec<i32>,
    lh: Vec<i32>,
    full_chord: String,
    rh_chord: String,
    lh_chord: String,
}

impl HarpHymnal {
    fn render_leadsheet(&self, ui: &mut egui::Ui, hymn: &Hymn, offset: i32, avail_w: f32) {
        // 1. Split beats into measures
        let measures = split_into_measures(hymn);
        if measures.is_empty() {
            return;
        }

        // 2. Estimate column width: monospace char at CELL_SIZE is ~7.5px.
        //    Max cell content is ~5 chars (e.g. "V^1|3"), so ~38px + 6px padding = 44px.
        let char_w = avail_w / 140.0; // scale relative to screen; ~7.5 at 1050px
        let cell_w = (char_w * 5.0 + 6.0).max(28.0);

        // 3. Calculate how many beats fit per line, then convert to measures
        let beats_per_line = ((avail_w - LABEL_W - 16.0) / (cell_w + 1.0)) as usize;

        // Greedily pack measures into lines by beat count
        let mut lines: Vec<(usize, usize)> = Vec::new();
        let mut line_start = 0;
        let mut line_beats = 0;

        for (i, m) in measures.iter().enumerate() {
            let m_beats = m.beats.len() + 1; // +1 for barline
            if line_beats + m_beats > beats_per_line && i > line_start {
                lines.push((line_start, i));
                line_start = i;
                line_beats = m_beats;
            } else {
                line_beats += m_beats;
            }
        }
        lines.push((line_start, measures.len()));

        // 4. Render each line as its own grid
        let row_labels = ["", "RH", "LH", "T", "2", "3", "4", "T", "2", "3", "4"];
        let num_rows = row_labels.len();

        for (line_idx, &(start, end)) in lines.iter().enumerate() {
            // Horizontal rule above each line group
            {
                let rect = ui.available_rect_before_wrap();
                let y = rect.top();
                ui.painter().line_segment(
                    [egui::Pos2::new(rect.left(), y), egui::Pos2::new(rect.left() + avail_w, y)],
                    egui::Stroke::new(1.0, BORDER),
                );
                ui.add_space(2.0);
            }

            // Count total columns for this line to compute even cell width
            let total_beats: usize = (start..end).map(|i| measures[i].beats.len()).sum();
            let total_bars = end - start;
            let usable = avail_w - LABEL_W - (total_bars as f32 * BAR_W) - 8.0;
            let col_w = (usable / total_beats as f32).min(52.0).max(20.0);

            let id = format!("ls_line_{}", line_idx);

            egui::Grid::new(id)
                .min_col_width(0.0)
                .spacing(egui::Vec2::new(0.0, 0.0))
                .show(ui, |ui| {
                    for row in 0..num_rows {
                        // Row label
                        ui.allocate_ui(egui::Vec2::new(LABEL_W, ROW_H), |ui| {
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                ui.add_space(2.0);
                                ui.label(egui::RichText::new(row_labels[row])
                                    .monospace().size(LABEL_SIZE).color(TEXT_MUTED));
                            });
                        });

                        for m_idx in start..end {
                            let measure = &measures[m_idx];

                            // Barline
                            ui.allocate_ui(egui::Vec2::new(BAR_W, ROW_H), |ui| {
                                if row > 0 {
                                    ui.centered_and_justified(|ui| {
                                        ui.label(egui::RichText::new("|")
                                            .monospace().size(CELL_SIZE).color(BAR_COLOR));
                                    });
                                }
                            });

                            // Beat columns
                            for beat in &measure.beats {
                                ui.allocate_ui(egui::Vec2::new(col_w, ROW_H), |ui| {
                                    let (text, color, size, bold) = beat_cell(beat, row, offset);
                                    if !text.is_empty() {
                                        let mut rt = egui::RichText::new(&text)
                                            .monospace().size(size).color(color);
                                        if bold { rt = rt.strong(); }
                                        ui.centered_and_justified(|ui| {
                                            ui.label(rt);
                                        });
                                    }
                                });
                            }
                        }

                        ui.end_row();
                    }
                });

            if line_idx < lines.len() - 1 {
                ui.add_space(2.0);
            }
        }
    }

    fn render_chord_table(&self, ui: &mut egui::Ui, hymn: &Hymn) {
        ui.label(egui::RichText::new("Chords in this hymn").strong().size(14.0).color(TEXT_SECONDARY));
        ui.add_space(4.0);

        // Header
        egui::Frame::NONE
            .fill(KEY_BADGE_BG)
            .inner_margin(egui::Margin::symmetric(8, 4))
            .corner_radius(egui::CornerRadius { nw: 6, ne: 6, sw: 0, se: 0 })
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    let hw = 90.0;
                    ui.allocate_ui(egui::Vec2::new(hw, 16.0), |ui| {
                        ui.label(egui::RichText::new("Harp").strong().color(egui::Color32::WHITE).size(11.0));
                    });
                    ui.allocate_ui(egui::Vec2::new(140.0, 16.0), |ui| {
                        ui.label(egui::RichText::new("LH / RH Strings").strong().color(egui::Color32::WHITE).size(11.0));
                    });
                    ui.allocate_ui(egui::Vec2::new(120.0, 16.0), |ui| {
                        ui.label(egui::RichText::new("Full Chord").strong().color(egui::Color32::WHITE).size(11.0));
                    });
                    ui.label(egui::RichText::new("Notes").strong().color(egui::Color32::WHITE).size(11.0));
                });
            });

        let mut seen = std::collections::HashSet::new();
        let mut row_idx = 0;
        for beat in &hymn.beats {
            if let Beat::Note { chord: Some(chord), rh_strings, lh_strings, rh_midi, lh_midi, full_chord, .. } = beat {
                if seen.insert(chord.clone()) {
                    let bg = if row_idx % 2 == 0 { CARD_BG } else { egui::Color32::from_rgb(245, 245, 242) };
                    egui::Frame::NONE
                        .fill(bg)
                        .inner_margin(egui::Margin::symmetric(8, 3))
                        .show(ui, |ui| {
                            ui.horizontal(|ui| {
                                let hw = 90.0;
                                ui.allocate_ui(egui::Vec2::new(hw, 16.0), |ui| {
                                    ui.label(egui::RichText::new(chord).monospace().color(ACCENT).strong().size(12.0));
                                });

                                let mut lh_s = lh_strings.clone(); lh_s.sort();
                                let mut rh_s = rh_strings.clone(); rh_s.sort();
                                let strings = format!("[{}] [{}]",
                                    lh_s.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(" "),
                                    rh_s.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(" "),
                                );
                                ui.allocate_ui(egui::Vec2::new(140.0, 16.0), |ui| {
                                    ui.label(egui::RichText::new(&strings).monospace().strong().size(11.0));
                                });

                                ui.allocate_ui(egui::Vec2::new(120.0, 16.0), |ui| {
                                    ui.label(egui::RichText::new(full_chord.as_deref().unwrap_or(""))
                                        .size(11.0).color(TEXT_PRIMARY));
                                });

                                let mut lh_n: Vec<i32> = lh_midi.clone(); lh_n.sort();
                                let mut rh_n: Vec<i32> = rh_midi.clone(); rh_n.sort();
                                let notes = format!("[{}] [{}]",
                                    lh_n.iter().map(|&m| midi_to_name(m)).collect::<Vec<_>>().join(" "),
                                    rh_n.iter().map(|&m| midi_to_name(m)).collect::<Vec<_>>().join(" "),
                                );
                                ui.label(egui::RichText::new(&notes).monospace().size(10.0).color(LH_FINGER));
                            });
                        });
                    row_idx += 1;
                }
            }
        }
    }
}

fn split_into_measures(hymn: &Hymn) -> Vec<Measure> {
    let mut measures: Vec<Measure> = Vec::new();
    let mut current_beats: Vec<BeatCol> = Vec::new();

    for beat in &hymn.beats {
        match beat {
            Beat::Bar => {
                if !current_beats.is_empty() {
                    measures.push(Measure { beats: std::mem::take(&mut current_beats) });
                }
            }
            Beat::Note { midi, chord, full_chord, rh_chord, lh_chord, rh_strings, lh_strings, .. } => {
                let melody_string = midi_to_harp_string(*midi, hymn.key_root)
                    .map(|s| to_relative_string(s, hymn.key_root));

                current_beats.push(BeatCol {
                    is_chord_change: chord.is_some(),
                    melody_string,
                    rh: rh_strings.clone(),
                    lh: lh_strings.clone(),
                    full_chord: full_chord.clone().unwrap_or_default(),
                    rh_chord: rh_chord.clone().unwrap_or_default(),
                    lh_chord: lh_chord.clone().unwrap_or_default(),
                });
            }
        }
    }
    // Don't forget trailing beats after last bar
    if !current_beats.is_empty() {
        measures.push(Measure { beats: current_beats });
    }
    measures
}

/// Get display text, color, size, and boldness for a beat cell at a given row.
fn beat_cell(beat: &BeatCol, row: usize, offset: i32) -> (String, egui::Color32, f32, bool) {
    match row {
        // Row 0: full chord name (gray, smaller)
        0 => {
            if beat.is_chord_change {
                (beat.full_chord.clone(), CHORD_LABEL, CHORD_SIZE, true)
            } else {
                (String::new(), CHORD_LABEL, CHORD_SIZE, false)
            }
        }
        // Row 1: RH chord (black, bold)
        1 => {
            if beat.is_chord_change {
                (beat.rh_chord.clone(), TEXT_PRIMARY, CELL_SIZE, true)
            } else {
                (String::new(), TEXT_PRIMARY, CELL_SIZE, false)
            }
        }
        // Row 2: LH chord (black, bold)
        2 => {
            if beat.is_chord_change {
                (beat.lh_chord.clone(), TEXT_PRIMARY, CELL_SIZE, true)
            } else {
                (String::new(), TEXT_PRIMARY, CELL_SIZE, false)
            }
        }
        // Rows 3-6: RH fingers (thumb=3 blue, rest dark)
        3..=6 => {
            let finger = row - 3;
            let color = if finger == 0 { ACCENT } else { RH_FINGER };
            if beat.is_chord_change && finger < beat.rh.len() {
                (format!("{}", beat.rh[finger] + offset), color, CELL_SIZE, finger == 0)
            } else if finger == 0 {
                let text = beat.melody_string
                    .map(|s| format!("{}", s + offset))
                    .unwrap_or_default();
                (text, ACCENT, CELL_SIZE, true)
            } else {
                (String::new(), color, CELL_SIZE, false)
            }
        }
        // Rows 7-10: LH fingers (thumb=7 blue, rest gray)
        7..=10 => {
            let finger = row - 7;
            let color = if finger == 0 { ACCENT } else { LH_FINGER };
            if beat.is_chord_change && finger < beat.lh.len() {
                (format!("{}", beat.lh[finger] + offset), color, CELL_SIZE, finger == 0)
            } else {
                (String::new(), color, CELL_SIZE, false)
            }
        }
        _ => (String::new(), TEXT_PRIMARY, CELL_SIZE, false),
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
