/// Music notation renderer using egui painter.
///
/// Draws treble+bass grand staff with realistic noteheads, stems, flags,
/// ledger lines, and barlines. Notes are positioned by beat time on the X axis
/// and staff position on the Y axis.

use eframe::egui;

use crate::abc::ScoreEvent;

// ── Layout constants ──
const STAFF_LINE_SPACING: f32 = 14.0;   // pixels between staff lines
const STAFF_LINES: usize = 5;
const GRAND_STAFF_GAP: f32 = 84.0;      // gap between treble and bass staves
const NOTEHEAD_RX: f32 = 8.0;           // notehead ellipse x-radius
const NOTEHEAD_RY: f32 = 5.5;           // notehead ellipse y-radius
const STEM_LENGTH: f32 = 42.0;
const FLAG_LENGTH: f32 = 18.0;
const LEDGER_EXTEND: f32 = 6.0;         // how far ledger lines extend past notehead
const BEAT_WIDTH: f32 = 60.0;           // pixels per beat (horizontal spacing)

// ── Colors ──
const STAFF_COLOR: egui::Color32 = egui::Color32::from_rgb(180, 180, 180);
const NOTE_COLOR: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const ACTIVE_COLOR: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const PLAYHEAD_COLOR: egui::Color32 = egui::Color32::from_rgba_premultiplied(37, 99, 235, 80);
const BAR_COLOR: egui::Color32 = egui::Color32::from_rgb(160, 160, 160);
const REST_COLOR: egui::Color32 = egui::Color32::from_rgb(120, 120, 120);
const CHORD_NAME_COLOR: egui::Color32 = egui::Color32::from_rgb(170, 170, 170);
const RH_LABEL_COLOR: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const LH_LABEL_COLOR: egui::Color32 = egui::Color32::from_rgb(42, 42, 42);
const THUMB_COLOR: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const FINGER_COLOR: egui::Color32 = egui::Color32::from_rgb(100, 100, 100);

// Label row heights above the treble staff
const LABEL_ROW_H: f32 = 19.0;
const NUM_LABEL_ROWS: usize = 3; // Chord, RH, LH
const LABEL_BLOCK_H: f32 = NUM_LABEL_ROWS as f32 * LABEL_ROW_H;

/// Convert MIDI note to staff position.
/// Middle C (MIDI 60) = 0. Each step = one diatonic position.
/// Positive = up from middle C, negative = down.
const DIATONIC_IN_OCTAVE: [i32; 12] = [0, 0, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6];

pub fn midi_to_staff_pos(midi: i32) -> i32 {
    let octave = midi / 12 - 5; // octave relative to middle C's octave
    let pc = (midi % 12) as usize;
    octave * 7 + DIATONIC_IN_OCTAVE[pc]
}

/// Y coordinate for a staff position on the treble staff.
/// Staff pos 0 = middle C = one ledger line below treble staff.
/// Treble staff bottom line = E4 (pos 2), top = F5 (pos 9).
pub fn treble_y(staff_pos: i32, treble_top_y: f32) -> f32 {
    // Top line of treble = F5 = pos 10
    // Each line is STAFF_LINE_SPACING apart
    let top_line_pos = 10; // F5
    treble_top_y + (top_line_pos - staff_pos) as f32 * (STAFF_LINE_SPACING / 2.0)
}

/// Y coordinate for a staff position on the bass staff.
/// Bass staff top line = A3 (pos -2), bottom = G2 (pos -10).
pub fn bass_y(staff_pos: i32, bass_top_y: f32) -> f32 {
    let top_line_pos = -2; // A3
    bass_top_y + (top_line_pos - staff_pos) as f32 * (STAFF_LINE_SPACING / 2.0)
}

/// Which staff a MIDI note belongs on. Treble for >= middle C, bass below.
fn is_treble(midi: i32) -> bool {
    midi >= 60
}

// Finger rows below the bass staff
const FINGER_ROW_H: f32 = 18.0;
const FINGER_GAP: f32 = 72.0; // gap between bass staff bottom and finger rows
const NUM_FINGER_ROWS: usize = 8; // RH: T 2 3 4, LH: T 2 3 4

pub struct NotationLayout {
    pub label_top_y: f32,  // top of the chord/RH/LH label rows
    pub treble_top_y: f32,
    pub bass_top_y: f32,
    pub finger_top_y: f32, // top of the finger number rows
    pub total_height: f32,
    pub left_margin: f32,
}

impl NotationLayout {
    pub fn new(top_y: f32, left_margin: f32) -> Self {
        let label_top_y = top_y;
        let treble_top_y = top_y + LABEL_BLOCK_H + 28.0; // gap between chord labels and staff
        let treble_height = (STAFF_LINES - 1) as f32 * STAFF_LINE_SPACING;
        let bass_top_y = treble_top_y + treble_height + GRAND_STAFF_GAP;
        let bass_height = (STAFF_LINES - 1) as f32 * STAFF_LINE_SPACING;
        let finger_top_y = bass_top_y + bass_height + FINGER_GAP;
        let finger_block = NUM_FINGER_ROWS as f32 * FINGER_ROW_H;
        let total_height = finger_top_y + finger_block - top_y + 10.0;

        Self { label_top_y, treble_top_y, bass_top_y, finger_top_y, total_height, left_margin }
    }

    /// Get Y position for a MIDI note.
    pub fn note_y(&self, midi: i32) -> f32 {
        let pos = midi_to_staff_pos(midi);
        if is_treble(midi) {
            treble_y(pos, self.treble_top_y)
        } else {
            bass_y(pos, self.bass_top_y)
        }
    }

    /// X position for a beat time.
    pub fn beat_x(&self, beat_time: f32, scroll_offset: f32) -> f32 {
        self.left_margin + beat_time * BEAT_WIDTH - scroll_offset
    }
}

/// Get the Y position for a MIDI note on the best staff (treble if >= middle C, bass otherwise).
fn note_y_on_best_staff(midi: i32, layout: &NotationLayout) -> f32 {
    let pos = midi_to_staff_pos(midi);
    if midi >= 60 {
        treble_y(pos, layout.treble_top_y)
    } else {
        bass_y(pos, layout.bass_top_y)
    }
}

/// Draw ledger lines for a note at the given staff position on the specified staff.
fn draw_ledger_lines(painter: &egui::Painter, x: f32, pos: i32, on_treble: bool, layout: &NotationLayout) {
    let lw = NOTEHEAD_RX + LEDGER_EXTEND;

    if on_treble {
        // Treble staff lines: pos 2 (bottom) to pos 10 (top)
        // Ledger below (middle C area)
        if pos <= 0 {
            let mut p = 0;
            while p >= pos {
                let y = treble_y(p, layout.treble_top_y);
                painter.line_segment(
                    [egui::Pos2::new(x - lw, y), egui::Pos2::new(x + lw, y)],
                    egui::Stroke::new(1.0, STAFF_COLOR),
                );
                p -= 2;
            }
        }
        // Ledger above
        if pos >= 12 {
            let mut p = 12;
            while p <= pos + 1 {
                let y = treble_y(p, layout.treble_top_y);
                painter.line_segment(
                    [egui::Pos2::new(x - lw, y), egui::Pos2::new(x + lw, y)],
                    egui::Stroke::new(1.0, STAFF_COLOR),
                );
                p += 2;
            }
        }
    } else {
        // Bass staff lines: pos -2 (top) to pos -10 (bottom)
        // Ledger above (middle C area)
        if pos >= 0 {
            let mut p = 0;
            while p <= pos + 1 {
                let y = bass_y(p, layout.bass_top_y);
                painter.line_segment(
                    [egui::Pos2::new(x - lw, y), egui::Pos2::new(x + lw, y)],
                    egui::Stroke::new(1.0, STAFF_COLOR),
                );
                p += 2;
            }
        }
        // Ledger below
        if pos <= -12 {
            let mut p = -12;
            while p >= pos {
                let y = bass_y(p, layout.bass_top_y);
                painter.line_segment(
                    [egui::Pos2::new(x - lw, y), egui::Pos2::new(x + lw, y)],
                    egui::Stroke::new(1.0, STAFF_COLOR),
                );
                p -= 2;
            }
        }
    }
}

/// Draw a notehead at (x, y).
fn draw_notehead(painter: &egui::Painter, x: f32, y: f32, filled: bool, color: egui::Color32) {
    let n_points = 16;
    let tilt = -0.3f32;
    let points: Vec<egui::Pos2> = (0..n_points).map(|i| {
        let angle = 2.0 * std::f32::consts::PI * i as f32 / n_points as f32;
        let ex = NOTEHEAD_RX * angle.cos();
        let ey = NOTEHEAD_RY * angle.sin();
        let rx = ex * tilt.cos() - ey * tilt.sin();
        let ry = ex * tilt.sin() + ey * tilt.cos();
        egui::Pos2::new(x + rx, y + ry)
    }).collect();

    if filled {
        painter.add(egui::Shape::convex_polygon(points, color, egui::Stroke::NONE));
    } else {
        painter.add(egui::Shape::convex_polygon(points, egui::Color32::TRANSPARENT, egui::Stroke::new(1.5, color)));
    }
}

/// Draw a note with stem UP (right side) — used for RH notes.
/// `pos` = staff position (0 = middle C), `on_treble` = which staff to draw on.
fn draw_note_stem_up(
    painter: &egui::Painter, layout: &NotationLayout,
    x: f32, pos: i32, on_treble: bool, beats: f32, is_active: bool,
) {
    let color = if is_active { ACTIVE_COLOR } else { NOTE_COLOR };
    let y = if on_treble { treble_y(pos, layout.treble_top_y) } else { bass_y(pos, layout.bass_top_y) };
    let filled = beats < 2.0;
    let has_stem = beats < 4.0;
    let flags = if beats <= 0.25 { 2 } else if beats <= 0.5 { 1 } else { 0 };

    draw_ledger_lines(painter, x, pos, on_treble, layout);
    draw_notehead(painter, x, y, filled, color);

    if has_stem {
        let stem_x = x + NOTEHEAD_RX - 0.5;
        let stem_y2 = y - STEM_LENGTH;
        painter.line_segment(
            [egui::Pos2::new(stem_x, y), egui::Pos2::new(stem_x, stem_y2)],
            egui::Stroke::new(1.2, color),
        );
        for f in 0..flags {
            let flag_y = stem_y2 + f as f32 * 6.0;
            let pts = [
                egui::Pos2::new(stem_x, flag_y),
                egui::Pos2::new(stem_x + 8.0, flag_y + FLAG_LENGTH * 0.5),
                egui::Pos2::new(stem_x + 4.0, flag_y + FLAG_LENGTH),
            ];
            painter.line_segment([pts[0], pts[1]], egui::Stroke::new(1.5, color));
            painter.line_segment([pts[1], pts[2]], egui::Stroke::new(1.5, color));
        }
    }
}

/// Draw a note with stem DOWN (left side) — used for LH notes.
/// `pos` = staff position (0 = middle C), `on_treble` = which staff to draw on.
fn draw_note_stem_down(
    painter: &egui::Painter, layout: &NotationLayout,
    x: f32, pos: i32, on_treble: bool, beats: f32, is_active: bool,
) {
    let color = if is_active { ACTIVE_COLOR } else { NOTE_COLOR };
    let y = if on_treble { treble_y(pos, layout.treble_top_y) } else { bass_y(pos, layout.bass_top_y) };
    let filled = beats < 2.0;
    let has_stem = beats < 4.0;
    let flags = if beats <= 0.25 { 2 } else if beats <= 0.5 { 1 } else { 0 };

    draw_ledger_lines(painter, x, pos, on_treble, layout);
    draw_notehead(painter, x, y, filled, color);

    if has_stem {
        let stem_x = x - NOTEHEAD_RX + 0.5;
        let stem_y2 = y + STEM_LENGTH;
        painter.line_segment(
            [egui::Pos2::new(stem_x, y), egui::Pos2::new(stem_x, stem_y2)],
            egui::Stroke::new(1.2, color),
        );
        for f in 0..flags {
            let flag_y = stem_y2 - f as f32 * 6.0;
            let pts = [
                egui::Pos2::new(stem_x, flag_y),
                egui::Pos2::new(stem_x - 8.0, flag_y - FLAG_LENGTH * 0.5),
                egui::Pos2::new(stem_x - 4.0, flag_y - FLAG_LENGTH),
            ];
            painter.line_segment([pts[0], pts[1]], egui::Stroke::new(1.5, color));
            painter.line_segment([pts[1], pts[2]], egui::Stroke::new(1.5, color));
        }
    }
}

/// Draw the grand staff (treble + bass, 5 lines each).
pub fn draw_staff(painter: &egui::Painter, layout: &NotationLayout, width: f32) {
    let x0 = 0.0;
    let x1 = width;

    // Treble staff
    for i in 0..STAFF_LINES {
        let y = layout.treble_top_y + i as f32 * STAFF_LINE_SPACING;
        painter.line_segment(
            [egui::Pos2::new(x0, y), egui::Pos2::new(x1, y)],
            egui::Stroke::new(1.0, STAFF_COLOR),
        );
    }

    // Bass staff
    for i in 0..STAFF_LINES {
        let y = layout.bass_top_y + i as f32 * STAFF_LINE_SPACING;
        painter.line_segment(
            [egui::Pos2::new(x0, y), egui::Pos2::new(x1, y)],
            egui::Stroke::new(1.0, STAFF_COLOR),
        );
    }

    // Brace/bracket connecting staves (simple vertical line)
    let top = layout.treble_top_y;
    let bottom = layout.bass_top_y + (STAFF_LINES - 1) as f32 * STAFF_LINE_SPACING;
    painter.line_segment(
        [egui::Pos2::new(x0, top), egui::Pos2::new(x0, bottom)],
        egui::Stroke::new(2.0, NOTE_COLOR),
    );
}

/// Draw a clef symbol (simplified text representation).
pub fn draw_clefs(painter: &egui::Painter, layout: &NotationLayout) {
    let treble_center_y = layout.treble_top_y + 2.0 * STAFF_LINE_SPACING;
    painter.text(
        egui::Pos2::new(8.0, treble_center_y),
        egui::Align2::LEFT_CENTER,
        "\u{1D11E}", // treble clef unicode (may not render, fallback below)
        egui::FontId::proportional(42.0),
        NOTE_COLOR,
    );

    let bass_center_y = layout.bass_top_y + 2.0 * STAFF_LINE_SPACING;
    painter.text(
        egui::Pos2::new(8.0, bass_center_y),
        egui::Align2::LEFT_CENTER,
        "\u{1D122}", // bass clef unicode
        egui::FontId::proportional(36.0),
        NOTE_COLOR,
    );
}


/// Draw a rest symbol at the given X position on the treble staff.
pub fn draw_rest(painter: &egui::Painter, x: f32, beats: f32, layout: &NotationLayout) {
    let center_y = layout.treble_top_y + 2.0 * STAFF_LINE_SPACING;
    let symbol = match beats as u8 {
        4 => "\u{1D13B}", // whole rest
        2 => "\u{1D13C}", // half rest
        1 => "\u{1D13D}", // quarter rest
        _ => "-",
    };
    painter.text(
        egui::Pos2::new(x, center_y),
        egui::Align2::CENTER_CENTER,
        symbol,
        egui::FontId::proportional(20.0),
        REST_COLOR,
    );
}

/// Draw a barline spanning both staves.
pub fn draw_barline(painter: &egui::Painter, x: f32, layout: &NotationLayout) {
    let top = layout.treble_top_y;
    let bottom = layout.bass_top_y + (STAFF_LINES - 1) as f32 * STAFF_LINE_SPACING;
    painter.line_segment(
        [egui::Pos2::new(x, top), egui::Pos2::new(x, bottom)],
        egui::Stroke::new(1.0, BAR_COLOR),
    );
}

/// Draw the playhead (vertical line at current position).
pub fn draw_playhead(painter: &egui::Painter, x: f32, layout: &NotationLayout) {
    let top = layout.treble_top_y - 10.0;
    let bottom = layout.bass_top_y + (STAFF_LINES - 1) as f32 * STAFF_LINE_SPACING + 10.0;

    // Glow/highlight band
    let rect = egui::Rect::from_min_max(
        egui::Pos2::new(x - 3.0, top),
        egui::Pos2::new(x + 3.0, bottom),
    );
    painter.rect_filled(rect, 0.0, PLAYHEAD_COLOR);

    // Sharp line
    painter.line_segment(
        [egui::Pos2::new(x, top), egui::Pos2::new(x, bottom)],
        egui::Stroke::new(2.0, ACTIVE_COLOR),
    );
}

/// Draw the static row labels on the left side.
fn draw_row_labels(painter: &egui::Painter, layout: &NotationLayout) {
    let labels = ["Chord", "RH", "LH"];
    let colors = [CHORD_NAME_COLOR, RH_LABEL_COLOR, LH_LABEL_COLOR];
    let x = 2.0;
    for (i, label) in labels.iter().enumerate() {
        let y = layout.label_top_y + i as f32 * LABEL_ROW_H + LABEL_ROW_H * 0.5;
        painter.text(
            egui::Pos2::new(x, y),
            egui::Align2::LEFT_CENTER,
            label,
            egui::FontId::monospace(14.0),
            colors[i],
        );
    }
}

/// Draw the static finger row labels below the bass staff.
fn draw_finger_labels(painter: &egui::Painter, layout: &NotationLayout) {
    let labels = ["T", "2", "3", "4", "T", "2", "3", "4"];
    let x = 2.0;
    for (i, label) in labels.iter().enumerate() {
        let y = layout.finger_top_y + i as f32 * FINGER_ROW_H + FINGER_ROW_H * 0.5;
        let color = match i {
            0 | 4 => THUMB_COLOR,
            _ => FINGER_COLOR,
        };
        // Add RH/LH group label on first row of each group
        if i == 0 {
            painter.text(
                egui::Pos2::new(x, y - FINGER_ROW_H * 0.1),
                egui::Align2::LEFT_CENTER,
                "RH",
                egui::FontId::monospace(10.0),
                RH_LABEL_COLOR,
            );
        } else if i == 4 {
            painter.text(
                egui::Pos2::new(x, y - FINGER_ROW_H * 0.1),
                egui::Align2::LEFT_CENTER,
                "LH",
                egui::FontId::monospace(10.0),
                LH_LABEL_COLOR,
            );
        }
        painter.text(
            egui::Pos2::new(x + 20.0, y),
            egui::Align2::LEFT_CENTER,
            label,
            egui::FontId::monospace(14.0),
            color,
        );
    }
}

/// Draw finger string numbers at a given X for a chord change.
fn draw_finger_values(
    painter: &egui::Painter,
    layout: &NotationLayout,
    x: f32,
    rh_strings: &[i32],
    lh_strings: &[i32],
    is_active: bool,
) {
    let rh_lh_gap = FINGER_ROW_H; // gap between RH and LH finger rows
    let rh_row_y = |row: usize| layout.finger_top_y + row as f32 * FINGER_ROW_H + FINGER_ROW_H * 0.5;
    let lh_row_y = |row: usize| layout.finger_top_y + (4.0 * FINGER_ROW_H) + rh_lh_gap + row as f32 * FINGER_ROW_H + FINGER_ROW_H * 0.5;

    // RH fingers: rows 0-3 (T, 2, 3, 4)
    for (i, &s) in rh_strings.iter().take(4).enumerate() {
        let color = if is_active { ACTIVE_COLOR }
                    else if i == 0 { THUMB_COLOR }
                    else { FINGER_COLOR };
        painter.text(
            egui::Pos2::new(x, rh_row_y(i)),
            egui::Align2::CENTER_CENTER,
            &s.to_string(),
            egui::FontId::monospace(15.0),
            color,
        );
    }

    // LH fingers (below RH with gap)
    for (i, &s) in lh_strings.iter().take(4).enumerate() {
        let color = if is_active { ACTIVE_COLOR }
                    else if i == 0 { THUMB_COLOR }
                    else { FINGER_COLOR };
        painter.text(
            egui::Pos2::new(x, lh_row_y(i)),
            egui::Align2::CENTER_CENTER,
            &s.to_string(),
            egui::FontId::monospace(15.0),
            color,
        );
    }
}

/// Draw label row values at a given X for a chord change.
fn draw_label_values(
    painter: &egui::Painter,
    layout: &NotationLayout,
    x: f32,
    chord_name: &str,
    rh_chord: &str,
    lh_chord: &str,
    _rh_strings: &[i32],
    _lh_strings: &[i32],
    is_active: bool,
) {
    let row_y = |row: usize| layout.label_top_y + row as f32 * LABEL_ROW_H + LABEL_ROW_H * 0.5;

    // Row 0: Chord name
    painter.text(
        egui::Pos2::new(x, row_y(0)),
        egui::Align2::CENTER_CENTER,
        chord_name,
        egui::FontId::monospace(14.0),
        if is_active { ACTIVE_COLOR } else { CHORD_NAME_COLOR },
    );

    // Row 1: RH chord
    painter.text(
        egui::Pos2::new(x, row_y(1)),
        egui::Align2::CENTER_CENTER,
        rh_chord,
        egui::FontId::monospace(15.0),
        if is_active { ACTIVE_COLOR } else { RH_LABEL_COLOR },
    );

    // Row 2: LH chord
    painter.text(
        egui::Pos2::new(x, row_y(2)),
        egui::Align2::CENTER_CENTER,
        lh_chord,
        egui::FontId::monospace(15.0),
        if is_active { ACTIVE_COLOR } else { LH_LABEL_COLOR },
    );
}

/// Render a full harp score.
/// Notes are positioned by string number → staff position.
/// RH → treble staff (stem up), LH → bass staff (stem down).
pub fn render_score(
    painter: &egui::Painter,
    layout: &NotationLayout,
    events: &[ScoreEvent],
    scroll_offset: f32,
    current_beat: f32,
    view_width: f32,
    key_root: i32,
) {
    // Convert relative string number (skips zero) to staff position (0 = middle C)
    let to_pos = |s: i32| -> i32 { crate::music::relative_string_to_staff_pos(s, key_root) };

    draw_staff(painter, layout, view_width);

    let playhead_x = layout.beat_x(current_beat, scroll_offset);
    if playhead_x > 0.0 && playhead_x < view_width {
        draw_playhead(painter, playhead_x, layout);
    }

    let mut beat_time: f32 = 0.0;

    for event in events {
        match event {
            ScoreEvent::Note { melody_string: _, display_rh, display_lh,
                                beats, chord_name, rh_chord, lh_chord, is_chord_change, .. } => {
                let x = layout.beat_x(beat_time + beats / 2.0, scroll_offset);
                let is_active = current_beat >= beat_time && current_beat < beat_time + beats;

                if x > -50.0 && x < view_width + 50.0 {
                    // Just draw what was pre-computed — no note selection here
                    let rh_used = display_rh;
                    let lh_used = display_lh;

                    if *is_chord_change && !rh_used.is_empty() {
                        for &s in rh_used {
                            let pos = to_pos(s);
                            draw_note_stem_up(painter, layout, x, pos, true, *beats, is_active);
                        }
                        for &s in lh_used {
                            let pos = to_pos(s);
                            draw_note_stem_down(painter, layout, x, pos, false, *beats, is_active);
                        }
                    } else if !rh_used.is_empty() {
                        // Melody only
                        let pos = to_pos(rh_used[0]);
                        draw_note_stem_up(painter, layout, x, pos, true, *beats, is_active);
                    }

                    if *is_chord_change {
                        draw_label_values(
                            painter, layout, x,
                            chord_name.as_deref().unwrap_or(""),
                            rh_chord.as_deref().unwrap_or(""),
                            lh_chord.as_deref().unwrap_or(""),
                            &rh_used,
                            &lh_used,
                            is_active,
                        );
                        draw_finger_values(
                            painter, layout, x,
                            &rh_used,
                            &lh_used,
                            is_active,
                        );
                    }
                }
                beat_time += beats;
            }
            ScoreEvent::Rest { beats } => {
                let x = layout.beat_x(beat_time + beats / 2.0, scroll_offset);
                if x > -50.0 && x < view_width + 50.0 {
                    draw_rest(painter, x, *beats, layout);
                }
                beat_time += beats;
            }
            ScoreEvent::Bar => {
                let x = layout.beat_x(beat_time, scroll_offset);
                if x > 0.0 && x < view_width {
                    draw_barline(painter, x, layout);
                }
            }
        }
    }
}
