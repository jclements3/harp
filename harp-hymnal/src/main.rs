mod abc;
mod chord;
mod music;
mod voicing;

use eframe::egui;
use abc::{Beat, Hymn, parse_abc};
use music::*;

// ── Colors ──
const BG: egui::Color32 = egui::Color32::from_rgb(248, 246, 241);
const ACCENT: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const TEXT_MUTED: egui::Color32 = egui::Color32::from_rgb(153, 153, 153);

const ALL_KEYS: [&str; 12] = ["C","Db","D","Eb","E","F","F#","G","Ab","A","Bb","B"];
const MODE_NAMES: [&str; 7] = ["Ionian","Dorian","Phrygian","Lydian","Mixolydian","Aeolian","Locrian"];

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([900.0, 700.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Harp Hymnal",
        options,
        Box::new(|cc| {
            apply_style(&cc.egui_ctx);
            Ok(Box::new(HarpHymnal::new()))
        }),
    )
}

fn apply_style(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.spacing.item_spacing = egui::Vec2::new(8.0, 4.0);
    style.spacing.button_padding = egui::Vec2::new(10.0, 4.0);
    style.visuals.window_fill = BG;
    style.visuals.panel_fill = BG;
    ctx.set_style(style);
}

struct HarpHymnal {
    hymns: Vec<Hymn>,
    current: usize,
    mode_offset: i32,
    current_key: String,
    search_text: String,
    status: String,
}

impl HarpHymnal {
    fn new() -> Self {
        let mut app = Self {
            hymns: Vec::new(),
            current: 0,
            mode_offset: 0,
            current_key: "C".into(),
            search_text: String::new(),
            status: "No file loaded. Use File > Open to load an ABC file.".into(),
        };

        // Try to auto-load OpenHymnal.abc from the app directory
        for path in &["data/OpenHymnal.abc", "data/harp_leadsheets.abc", "OpenHymnal.abc", "harp_leadsheets.abc"] {
            if let Ok(text) = std::fs::read_to_string(path) {
                app.load_abc(&text);
                break;
            }
            // Also try relative to executable
            if let Ok(exe) = std::env::current_exe() {
                if let Some(dir) = exe.parent() {
                    let p = dir.join(path);
                    if let Ok(text) = std::fs::read_to_string(&p) {
                        app.load_abc(&text);
                        break;
                    }
                }
            }
        }

        app
    }

    fn load_abc(&mut self, text: &str) {
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
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open ABC...").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("ABC files", &["abc"])
                            .pick_file()
                        {
                            if let Ok(text) = std::fs::read_to_string(&path) {
                                self.load_abc(&text);
                            }
                        }
                        ui.close_menu();
                    }
                });
            });
        });

        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.colored_label(TEXT_MUTED, &self.status);
        });

        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                // Key badge
                let key_text = self.key_display();
                ui.label(egui::RichText::new(&key_text).strong().monospace().size(16.0));

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

                ui.separator();

                // Mode buttons
                for (i, name) in MODE_NAMES.iter().enumerate() {
                    let selected = self.mode_offset == i as i32;
                    if ui.selectable_label(selected, *name).clicked() {
                        self.mode_offset = i as i32;
                    }
                }
            });

            ui.horizontal(|ui| {
                // Hymn selector
                if !self.hymns.is_empty() {
                    let label = self.current_hymn()
                        .map(|h| format!("{}. {}", h.number, h.title))
                        .unwrap_or_default();

                    egui::ComboBox::from_id_salt("hymn_select")
                        .selected_text(&label)
                        .width(400.0)
                        .show_ui(ui, |ui| {
                            for (i, h) in self.hymns.iter().enumerate() {
                                let text = format!("{}. {}", h.number, h.title);
                                if ui.selectable_value(&mut self.current, i, &text).clicked() {
                                    self.current_key = h.key.clone();
                                }
                            }
                        });

                    ui.separator();

                    // Search
                    let response = ui.add(
                        egui::TextEdit::singleline(&mut self.search_text)
                            .hint_text("Search...")
                            .desired_width(200.0)
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
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            if self.hymns.is_empty() {
                ui.centered_and_justified(|ui| {
                    ui.label("Load an ABC file to begin.");
                });
                return;
            }

            let hymn = self.hymns[self.current].clone();
            let offset = self.mode_offset;

            egui::ScrollArea::both().show(ui, |ui| {
                // Title
                ui.heading(egui::RichText::new(format!("{}. {}", hymn.number, hymn.title)).strong());
                ui.add_space(8.0);

                // Leadsheet grid
                self.render_leadsheet(ui, &hymn, offset);

                ui.add_space(16.0);

                // Chord table
                self.render_chord_table(ui, &hymn);
            });
        });
    }
}

impl HarpHymnal {
    fn render_leadsheet(&self, ui: &mut egui::Ui, hymn: &Hymn, offset: i32) {
        // Collect columns
        struct Col {
            is_bar: bool,
            is_chord_change: bool,
            melody_string: Option<i32>,
            rh: Vec<i32>,
            lh: Vec<i32>,
            full_chord: String,
            rh_chord: String,
            lh_chord: String,
        }

        let mut cols: Vec<Col> = Vec::new();

        for beat in &hymn.beats {
            match beat {
                Beat::Bar => {
                    cols.push(Col {
                        is_bar: true, is_chord_change: false,
                        melody_string: None, rh: vec![], lh: vec![],
                        full_chord: String::new(), rh_chord: String::new(), lh_chord: String::new(),
                    });
                }
                Beat::Note { midi, chord, full_chord, rh_chord, lh_chord, rh_strings, lh_strings, .. } => {
                    let melody_string = midi_to_harp_string(*midi, hymn.key_root)
                        .map(|s| to_relative_string(s, hymn.key_root));

                    cols.push(Col {
                        is_bar: false,
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

        let row_labels = ["Chord", "RH", "LH", "thumb", "index", "middle", "ring", "thumb", "index", "middle", "ring"];

        egui::Grid::new("leadsheet")
            .min_col_width(4.0)
            .spacing(egui::Vec2::new(1.0, 0.0))
            .show(ui, |ui| {
                for (row, label) in row_labels.iter().enumerate() {
                    // Row label
                    ui.label(egui::RichText::new(*label).small().color(TEXT_MUTED));

                    for col in &cols {
                        if col.is_bar {
                            ui.label(egui::RichText::new("|").monospace().color(egui::Color32::DARK_GRAY));
                            continue;
                        }

                        let text = match row {
                            0 => { // full chord
                                if col.is_chord_change { col.full_chord.clone() } else { String::new() }
                            }
                            1 => { // RH chord
                                if col.is_chord_change { col.rh_chord.clone() } else { String::new() }
                            }
                            2 => { // LH chord
                                if col.is_chord_change { col.lh_chord.clone() } else { String::new() }
                            }
                            3..=6 => { // RH fingers (thumb=3, index=4, ...)
                                let finger = row - 3;
                                if col.is_chord_change && finger < col.rh.len() {
                                    format!("{}", col.rh[finger] + offset)
                                } else if finger == 0 {
                                    col.melody_string.map(|s| format!("{}", s + offset)).unwrap_or_default()
                                } else {
                                    String::new()
                                }
                            }
                            7..=10 => { // LH fingers
                                let finger = row - 7;
                                if col.is_chord_change && finger < col.lh.len() {
                                    format!("{}", col.lh[finger] + offset)
                                } else {
                                    String::new()
                                }
                            }
                            _ => String::new(),
                        };

                        let color = match row {
                            0 => TEXT_MUTED,
                            1 | 2 => egui::Color32::BLACK,
                            3 | 7 => ACCENT, // thumbs
                            4..=6 => egui::Color32::DARK_GRAY,
                            8..=10 => egui::Color32::GRAY,
                            _ => egui::Color32::BLACK,
                        };

                        let rt = egui::RichText::new(&text).monospace().size(12.0).color(color);
                        let rt = if matches!(row, 1 | 2) { rt.strong() } else { rt };
                        ui.label(rt);
                    }
                    ui.end_row();
                }
            });
    }

    fn render_chord_table(&self, ui: &mut egui::Ui, hymn: &Hymn) {
        ui.label(egui::RichText::new("Chords in this hymn").strong());

        egui::Grid::new("chord_table")
            .striped(true)
            .min_col_width(60.0)
            .show(ui, |ui| {
                ui.label(egui::RichText::new("Harp").strong());
                ui.label(egui::RichText::new("LH RH Strings").strong());
                ui.label(egui::RichText::new("Full Chord").strong());
                ui.label(egui::RichText::new("Notes").strong());
                ui.end_row();

                let mut seen = std::collections::HashSet::new();
                for beat in &hymn.beats {
                    if let Beat::Note { chord: Some(chord), rh_strings, lh_strings, rh_midi, lh_midi, full_chord, .. } = beat {
                        if seen.insert(chord.clone()) {
                            ui.label(egui::RichText::new(chord).monospace().color(ACCENT).strong());

                            let mut lh_s = lh_strings.clone(); lh_s.sort();
                            let mut rh_s = rh_strings.clone(); rh_s.sort();
                            let strings = format!("[{}] [{}]",
                                lh_s.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(" "),
                                rh_s.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(" "),
                            );
                            ui.label(egui::RichText::new(&strings).monospace());

                            ui.label(full_chord.as_deref().unwrap_or(""));

                            let mut lh_n: Vec<i32> = lh_midi.clone(); lh_n.sort();
                            let mut rh_n: Vec<i32> = rh_midi.clone(); rh_n.sort();
                            let notes = format!("[{}] [{}]",
                                lh_n.iter().map(|&m| midi_to_name(m)).collect::<Vec<_>>().join(" "),
                                rh_n.iter().map(|&m| midi_to_name(m)).collect::<Vec<_>>().join(" "),
                            );
                            ui.label(egui::RichText::new(&notes).monospace().small());

                            ui.end_row();
                        }
                    }
                }
            });
    }
}
