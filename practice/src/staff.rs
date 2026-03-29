/// Scrolling staff renderer — no clefs, clean notation.

use eframe::egui;
use crate::music;
use crate::drill::TimedChallenge;

const LINE_SP: f32 = 11.0;
const STAFF_LINES: usize = 5;
const STAFF_GAP: f32 = 48.0;
const NOTE_RX: f32 = 7.5;
const NOTE_RY: f32 = 5.0;
const STEM_H: f32 = 32.0;
const LEDGER_EXT: f32 = 5.0;
const CHALLENGE_GAP: f32 = 110.0;
const PLAYHEAD_FRAC: f32 = 0.20;

const STAFF_C: egui::Color32 = egui::Color32::from_rgb(190, 185, 175);
const NOTE_C: egui::Color32 = egui::Color32::from_rgb(50, 50, 50);
const ACTIVE_C: egui::Color32 = egui::Color32::from_rgb(37, 99, 235);
const DONE_C: egui::Color32 = egui::Color32::from_rgb(40, 167, 69);
const UPCOMING_C: egui::Color32 = egui::Color32::from_rgb(160, 160, 160);
const PLAYHEAD_GLOW: egui::Color32 = egui::Color32::from_rgba_premultiplied(37, 99, 235, 50);
const LABEL_C: egui::Color32 = egui::Color32::from_rgb(100, 100, 100);
const TARGET_BG: egui::Color32 = egui::Color32::from_rgba_premultiplied(37, 99, 235, 15);

const DIAT: [i32; 12] = [0, 0, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6];

fn midi_staff_pos(midi: i32) -> i32 {
    (midi / 12 - 5) * 7 + DIAT[(midi % 12) as usize]
}

pub struct Layout {
    pub treble_top: f32,
    pub bass_top: f32,
}

impl Layout {
    pub fn new(rect: &egui::Rect) -> Self {
        let treble_h = (STAFF_LINES - 1) as f32 * LINE_SP;
        let bass_h = treble_h;
        let total = treble_h + STAFF_GAP + bass_h;
        let top = rect.center().y - total / 2.0;
        Self { treble_top: top, bass_top: top + treble_h + STAFF_GAP }
    }
    fn treble_y(&self, pos: i32) -> f32 { self.treble_top + (10 - pos) as f32 * (LINE_SP / 2.0) }
    fn bass_y(&self, pos: i32) -> f32 { self.bass_top + (-2 - pos) as f32 * (LINE_SP / 2.0) }
}

fn notehead(p: &egui::Painter, x: f32, y: f32, color: egui::Color32) {
    let tilt = -0.3f32;
    let pts: Vec<egui::Pos2> = (0..16).map(|i| {
        let a = std::f32::consts::TAU * i as f32 / 16.0;
        let ex = NOTE_RX * a.cos();
        let ey = NOTE_RY * a.sin();
        egui::Pos2::new(x + ex * tilt.cos() - ey * tilt.sin(), y + ex * tilt.sin() + ey * tilt.cos())
    }).collect();
    p.add(egui::Shape::convex_polygon(pts, color, egui::Stroke::NONE));
}

fn ledgers(p: &egui::Painter, x: f32, pos: i32, treble: bool, ly: &Layout) {
    let w = NOTE_RX + LEDGER_EXT;
    let s = egui::Stroke::new(1.0, STAFF_C);
    if treble {
        let mut pp = 0; while pp >= pos { p.line_segment([egui::Pos2::new(x-w, ly.treble_y(pp)), egui::Pos2::new(x+w, ly.treble_y(pp))], s); pp -= 2; }
        let mut pp = 12; while pp <= pos+1 { p.line_segment([egui::Pos2::new(x-w, ly.treble_y(pp)), egui::Pos2::new(x+w, ly.treble_y(pp))], s); pp += 2; }
    } else {
        let mut pp = 0; while pp <= pos+1 { p.line_segment([egui::Pos2::new(x-w, ly.bass_y(pp)), egui::Pos2::new(x+w, ly.bass_y(pp))], s); pp += 2; }
        let mut pp = -12; while pp >= pos { p.line_segment([egui::Pos2::new(x-w, ly.bass_y(pp)), egui::Pos2::new(x+w, ly.bass_y(pp))], s); pp -= 2; }
    }
}

fn draw_note(p: &egui::Painter, ly: &Layout, x: f32, midi: i32, color: egui::Color32) {
    let pos = midi_staff_pos(midi);
    let on_treble = midi >= 60;
    let y = if on_treble { ly.treble_y(pos) } else { ly.bass_y(pos) };
    ledgers(p, x, pos, on_treble, ly);
    notehead(p, x, y, color);
    if on_treble {
        p.line_segment([egui::Pos2::new(x+NOTE_RX-0.5, y), egui::Pos2::new(x+NOTE_RX-0.5, y-STEM_H)], egui::Stroke::new(1.2, color));
    } else {
        p.line_segment([egui::Pos2::new(x-NOTE_RX+0.5, y), egui::Pos2::new(x-NOTE_RX+0.5, y+STEM_H)], egui::Stroke::new(1.2, color));
    }
}

pub fn draw_staff_lines(p: &egui::Painter, ly: &Layout, w: f32) {
    let s = egui::Stroke::new(1.0, STAFF_C);
    for i in 0..STAFF_LINES {
        let y = ly.treble_top + i as f32 * LINE_SP;
        p.line_segment([egui::Pos2::new(0.0, y), egui::Pos2::new(w, y)], s);
    }
    for i in 0..STAFF_LINES {
        let y = ly.bass_top + i as f32 * LINE_SP;
        p.line_segment([egui::Pos2::new(0.0, y), egui::Pos2::new(w, y)], s);
    }
    // Brace connecting staves
    let top = ly.treble_top;
    let bot = ly.bass_top + (STAFF_LINES-1) as f32 * LINE_SP;
    p.line_segment([egui::Pos2::new(0.0, top), egui::Pos2::new(0.0, bot)], egui::Stroke::new(2.0, NOTE_C));
}

fn draw_challenge(p: &egui::Painter, ly: &Layout, x: f32, tc: &TimedChallenge, current: bool, past: bool) {
    let ch = &tc.challenge;
    if current {
        let top = ly.treble_top - 20.0;
        let bot = ly.bass_top + (STAFF_LINES-1) as f32 * LINE_SP + 20.0;
        p.rect_filled(egui::Rect::from_min_max(egui::Pos2::new(x-30.0, top), egui::Pos2::new(x+30.0, bot)), 6.0, TARGET_BG);
    }
    let label_y = ly.treble_top - 8.0;
    let (lc, ls) = if current { (ACTIVE_C, 13.0) } else if past { (egui::Color32::from_rgba_premultiplied(100,100,100,60), 11.0) } else { (LABEL_C, 11.0) };
    p.text(egui::Pos2::new(x, label_y), egui::Align2::CENTER_BOTTOM, &ch.label, egui::FontId::proportional(ls), lc);

    for (i, &midi) in ch.notes.iter().enumerate() {
        let color = if past { egui::Color32::from_rgba_premultiplied(DONE_C.r(), DONE_C.g(), DONE_C.b(), 40) }
            else if current && i == ch.next_finger { ACTIVE_C }
            else if current && i < ch.next_finger { egui::Color32::from_rgba_premultiplied(DONE_C.r(), DONE_C.g(), DONE_C.b(), 120) }
            else if current { NOTE_C }
            else { UPCOMING_C };
        draw_note(p, ly, x, midi as i32, color);
        if current && i == ch.next_finger {
            let y = if midi >= 60 { ly.treble_y(midi_staff_pos(midi as i32)) } else { ly.bass_y(midi_staff_pos(midi as i32)) };
            let nm = music::NOTE_NAMES[midi as usize % 12];
            let oc = (midi as i32 / 12) - 1;
            p.text(egui::Pos2::new(x + NOTE_RX + 10.0, y), egui::Align2::LEFT_CENTER, format!("{nm}{oc}"), egui::FontId::monospace(11.0), ACTIVE_C);
        }
    }
}

pub fn render_scrolling(p: &egui::Painter, rect: egui::Rect, timeline: &[TimedChallenge], current_idx: usize, elapsed: f32, _key_root: u8) {
    let ly = Layout::new(&rect);
    draw_staff_lines(p, &ly, rect.width());
    let ph_x = rect.left() + rect.width() * PLAYHEAD_FRAC;
    let frac = timeline.get(current_idx).map_or(0.0, |tc| ((elapsed - tc.start_time) / tc.duration).clamp(0.0, 1.0));
    for (i, tc) in timeline.iter().enumerate() {
        let x = ph_x + (i as f32 - current_idx as f32 - frac) * CHALLENGE_GAP;
        if x < rect.left() - 40.0 || x > rect.right() + 40.0 { continue; }
        draw_challenge(p, &ly, x, tc, i == current_idx, i < current_idx);
    }
    let top = ly.treble_top - 20.0;
    let bot = ly.bass_top + (STAFF_LINES-1) as f32 * LINE_SP + 20.0;
    p.rect_filled(egui::Rect::from_min_max(egui::Pos2::new(ph_x-3.0, top), egui::Pos2::new(ph_x+3.0, bot)), 0.0, PLAYHEAD_GLOW);
    p.line_segment([egui::Pos2::new(ph_x, top), egui::Pos2::new(ph_x, bot)], egui::Stroke::new(2.0, ACTIVE_C));
}

pub fn draw_detected_note(p: &egui::Painter, rect: egui::Rect, midi: u8, _key_root: u8) {
    let ly = Layout::new(&rect);
    let x = rect.left() + rect.width() * PLAYHEAD_FRAC + 30.0;
    draw_note(p, &ly, x, midi as i32, egui::Color32::from_rgba_premultiplied(156, 163, 175, 90));
}
