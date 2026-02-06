use eframe::egui::{self, Color32, Pos2, Rect, Stroke};

use crate::music;

// ── Layout constants (matching the HTML version) ───────────────────

/// Half-space height on the staff
const STEP: f32 = 9.0;
/// Horizontal radius of a notehead
const NOTE_RX: f32 = STEP * 1.1;
/// Vertical radius of a notehead
const NOTE_RY: f32 = STEP * 0.85;

/// Treble staff lines as diatonic positions: E4=30, G4=32, B4=34, D5=36, F5=38
const TREBLE_LINES: [i32; 5] = [30, 32, 34, 36, 38];
/// Bass staff lines as diatonic positions: G2=18, B2=20, D3=22, F3=24, A3=26
const BASS_LINES: [i32; 5] = [18, 20, 22, 24, 26];

/// Middle C staff position = 28
const MIDDLE_C_POS: i32 = 28;

const STAFF_COLOR: Color32 = Color32::from_rgb(176, 171, 160);
const LABEL_COLOR: Color32 = Color32::from_rgb(153, 153, 153);

/// Convert a diatonic staff position to a Y coordinate.
/// `center_y` is the Y pixel for middle C (pos=28).
fn pos_to_y(pos: i32, center_y: f32) -> f32 {
    center_y - (pos - MIDDLE_C_POS) as f32 * STEP
}

// ── Public drawing API ─────────────────────────────────────────────

/// Draw the grand staff (5 treble lines + 5 bass lines + clef labels).
/// `rect` is the available area; the staff is centered vertically in it.
pub fn draw_staff(painter: &egui::Painter, rect: Rect) {
    let center_y = rect.center().y;
    let left = rect.left() + 10.0;
    let right = rect.right() - 10.0;
    let stroke = Stroke::new(1.0, STAFF_COLOR);

    for &pos in &TREBLE_LINES {
        let y = pos_to_y(pos, center_y);
        painter.line_segment([Pos2::new(left, y), Pos2::new(right, y)], stroke);
    }
    for &pos in &BASS_LINES {
        let y = pos_to_y(pos, center_y);
        painter.line_segment([Pos2::new(left, y), Pos2::new(right, y)], stroke);
    }

    // Clef labels
    let treble_center_y = pos_to_y(34, center_y);
    let bass_center_y = pos_to_y(22, center_y);
    let label_x = left - 2.0;

    painter.text(
        Pos2::new(label_x, treble_center_y),
        egui::Align2::RIGHT_CENTER,
        "\u{1D11E}",
        egui::FontId::proportional(28.0),
        LABEL_COLOR,
    );
    painter.text(
        Pos2::new(label_x, bass_center_y),
        egui::Align2::RIGHT_CENTER,
        "\u{1D122}",
        egui::FontId::proportional(28.0),
        LABEL_COLOR,
    );
}

/// Draw ledger lines for a note at the given staff position.
pub fn draw_ledger_lines(painter: &egui::Painter, center_x: f32, center_y_staff: f32, staff_pos: i32) {
    let stroke = Stroke::new(1.0, STAFF_COLOR);
    let hw = NOTE_RX + 4.0;

    let mut lines = Vec::new();

    // Above treble staff
    let mut p = 40;
    while p <= staff_pos {
        lines.push(p);
        p += 2;
    }
    // Below bass staff
    p = 16;
    while p >= staff_pos {
        lines.push(p);
        p -= 2;
    }
    // Middle C
    if staff_pos == MIDDLE_C_POS {
        lines.push(MIDDLE_C_POS);
    }
    // Between staves (above bass, below treble)
    if staff_pos > 26 && staff_pos < 30 && staff_pos % 2 == 0 {
        lines.push(staff_pos);
    }

    for lp in lines {
        let y = pos_to_y(lp, center_y_staff);
        painter.line_segment(
            [Pos2::new(center_x - hw, y), Pos2::new(center_x + hw, y)],
            stroke,
        );
    }
}

/// Draw a shaped notehead at `center` with the given scale degree (1-7, 0=chromatic).
pub fn draw_shaped_note(painter: &egui::Painter, center: Pos2, rx: f32, ry: f32, degree: u8, color: Color32) {
    match degree {
        1 => {
            // Do → triangle (point up)
            let pts = vec![
                Pos2::new(center.x, center.y - ry),
                Pos2::new(center.x + rx, center.y + ry),
                Pos2::new(center.x - rx, center.y + ry),
            ];
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        2 => {
            // Re → half-moon (right half of ellipse)
            let n = 20;
            let mut pts = Vec::with_capacity(n + 1);
            for i in 0..=n {
                let angle = -std::f32::consts::FRAC_PI_2
                    + std::f32::consts::PI * (i as f32 / n as f32);
                pts.push(Pos2::new(
                    center.x + rx * angle.cos(),
                    center.y + ry * angle.sin(),
                ));
            }
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        3 => {
            // Mi → diamond
            let pts = vec![
                Pos2::new(center.x, center.y - ry),
                Pos2::new(center.x + rx, center.y),
                Pos2::new(center.x, center.y + ry),
                Pos2::new(center.x - rx, center.y),
            ];
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        4 => {
            // Fa → right-triangle (flag-like)
            let pts = vec![
                Pos2::new(center.x - rx, center.y - ry),
                Pos2::new(center.x + rx, center.y + ry),
                Pos2::new(center.x - rx, center.y + ry),
            ];
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        5 => {
            // Sol → ellipse/circle
            let n = 24;
            let mut pts = Vec::with_capacity(n);
            for i in 0..n {
                let angle = 2.0 * std::f32::consts::PI * (i as f32 / n as f32);
                pts.push(Pos2::new(
                    center.x + rx * angle.cos(),
                    center.y + ry * angle.sin(),
                ));
            }
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        6 => {
            // La → square/rect
            let pts = vec![
                Pos2::new(center.x - rx, center.y - ry),
                Pos2::new(center.x + rx, center.y - ry),
                Pos2::new(center.x + rx, center.y + ry),
                Pos2::new(center.x - rx, center.y + ry),
            ];
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        7 => {
            // Ti → inverted triangle (point down)
            let pts = vec![
                Pos2::new(center.x, center.y + ry),
                Pos2::new(center.x + rx, center.y - ry),
                Pos2::new(center.x - rx, center.y - ry),
            ];
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
        _ => {
            // Chromatic → small muted circle
            let n = 16;
            let sr = rx * 0.5;
            let sy = ry * 0.5;
            let mut pts = Vec::with_capacity(n);
            for i in 0..n {
                let angle = 2.0 * std::f32::consts::PI * (i as f32 / n as f32);
                pts.push(Pos2::new(
                    center.x + sr * angle.cos(),
                    center.y + sy * angle.sin(),
                ));
            }
            painter.add(egui::Shape::convex_polygon(pts, color, Stroke::NONE));
        }
    }
}

/// Draw a note on the grand staff: staff position + ledger lines + shaped notehead + sharp indicator.
pub fn draw_note_on_staff(
    painter: &egui::Painter,
    rect: Rect,
    midi: u8,
    key_root: u8,
    color: Color32,
    alpha: f32,
) {
    let center_y = rect.center().y;
    let center_x = rect.center().x;
    let staff_pos = music::midi_to_staff_pos(midi);
    let y = pos_to_y(staff_pos, center_y);
    let degree = music::scale_degree(midi, key_root);

    let c = Color32::from_rgba_premultiplied(
        (color.r() as f32 * alpha) as u8,
        (color.g() as f32 * alpha) as u8,
        (color.b() as f32 * alpha) as u8,
        (255.0 * alpha) as u8,
    );

    // Ledger lines
    draw_ledger_lines(painter, center_x, center_y, staff_pos);

    // Sharp indicator
    if music::IS_SHARP[midi as usize % 12] {
        painter.text(
            Pos2::new(center_x - NOTE_RX - 8.0, y),
            egui::Align2::RIGHT_CENTER,
            "\u{266f}",
            egui::FontId::proportional(14.0),
            c,
        );
    }

    // Shaped notehead
    draw_shaped_note(painter, Pos2::new(center_x, y), NOTE_RX, NOTE_RY, degree, c);
}

/// Draw a hand position on the grand staff: all notes laid out left-to-right,
/// completed notes dimmed green, current note bright, upcoming notes dimmed grey.
pub fn draw_hand_position(
    painter: &egui::Painter,
    rect: Rect,
    position: &[u8],
    current_finger: usize,
    key_root: u8,
) {
    if position.is_empty() {
        return;
    }

    let center_y = rect.center().y;
    let center_x = rect.center().x;
    let spacing = NOTE_RX * 3.5;
    let total_width = (position.len() as f32 - 1.0) * spacing;
    let start_x = center_x - total_width / 2.0;

    for (i, &midi) in position.iter().enumerate() {
        let x = start_x + i as f32 * spacing;
        let staff_pos = music::midi_to_staff_pos(midi);
        let y = pos_to_y(staff_pos, center_y);
        let degree = music::scale_degree(midi, key_root);

        let (color, alpha) = if i == current_finger {
            // Current target — bright blue
            (Color32::from_rgb(37, 99, 235), 1.0)
        } else if i < current_finger {
            // Already played — green, dimmed
            (Color32::from_rgb(40, 167, 69), 0.35)
        } else {
            // Upcoming — dark grey, dimmed
            (Color32::from_rgb(90, 90, 90), 0.4)
        };

        let c = Color32::from_rgba_premultiplied(
            (color.r() as f32 * alpha) as u8,
            (color.g() as f32 * alpha) as u8,
            (color.b() as f32 * alpha) as u8,
            (255.0 * alpha) as u8,
        );

        // Ledger lines
        draw_ledger_lines(painter, x, center_y, staff_pos);

        // Sharp indicator
        if music::IS_SHARP[midi as usize % 12] {
            painter.text(
                Pos2::new(x - NOTE_RX - 8.0, y),
                egui::Align2::RIGHT_CENTER,
                "\u{266f}",
                egui::FontId::proportional(14.0),
                c,
            );
        }

        // Shaped notehead
        draw_shaped_note(painter, Pos2::new(x, y), NOTE_RX, NOTE_RY, degree, c);

        // Finger number below/above the note
        let finger_num = i + 1;
        let label_y = if staff_pos < MIDDLE_C_POS {
            y + NOTE_RY + 10.0
        } else {
            y - NOTE_RY - 10.0
        };
        let label_alpha = if i == current_finger { 0.9 } else { 0.3 };
        let label_color = Color32::from_rgba_premultiplied(
            (90.0 * label_alpha) as u8,
            (90.0 * label_alpha) as u8,
            (90.0 * label_alpha) as u8,
            (255.0 * label_alpha) as u8,
        );
        painter.text(
            Pos2::new(x, label_y),
            egui::Align2::CENTER_CENTER,
            format!("{finger_num}"),
            egui::FontId::monospace(10.0),
            label_color,
        );
    }
}

/// Draw the shape-to-solfege legend bar at the top of the staff area.
pub fn draw_legend(painter: &egui::Painter, rect: Rect, key_root: u8) {
    let y = rect.top() + 14.0;
    let start_x = rect.left() + 20.0;
    let spacing = 52.0;
    let note_color = Color32::from_rgb(58, 58, 58);
    let text_color = LABEL_COLOR;

    for d in 1u8..=7 {
        let ox = start_x + (d as f32 - 1.0) * spacing;
        draw_shaped_note(
            painter,
            Pos2::new(ox + 7.0, y),
            6.0,
            5.0,
            d,
            note_color,
        );
        painter.text(
            Pos2::new(ox + 16.0, y),
            egui::Align2::LEFT_CENTER,
            music::SOLFEGE[(d - 1) as usize],
            egui::FontId::monospace(10.0),
            text_color,
        );
    }

    let _ = key_root; // key_root available for future use
}
