from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, HRFlowable
from reportlab.lib.enums import TA_CENTER, TA_LEFT
import os
OUTPUT = "/home/james.clements/projects/harp/harp_practice_drill.pdf"

doc = SimpleDocTemplate(
    OUTPUT,
    pagesize=letter,
    leftMargin=0.65*inch,
    rightMargin=0.65*inch,
    topMargin=0.55*inch,
    bottomMargin=0.55*inch
)

title_style = ParagraphStyle('T', fontName='Helvetica-Bold',
    fontSize=16, textColor=colors.HexColor('#1a1a1a'),
    spaceAfter=10, alignment=TA_CENTER)
subtitle_style = ParagraphStyle('ST', fontName='Helvetica',
    fontSize=9, textColor=colors.HexColor('#666666'),
    spaceAfter=6, alignment=TA_CENTER)
section_style = ParagraphStyle('S', fontName='Helvetica-Bold',
    fontSize=9, textColor=colors.HexColor('#333333'),
    spaceAfter=3, spaceBefore=6)
obs_style = ParagraphStyle('O', fontName='Helvetica',
    fontSize=8, textColor=colors.HexColor('#444444'),
    spaceAfter=2, leading=11)

mode_colors = {
    'Ionian':      colors.HexColor('#d4edda'),
    'Dorian':      colors.HexColor('#d1ecf1'),
    'Phrygian':    colors.HexColor('#fff3cd'),
    'Lydian':      colors.HexColor('#e2d9f3'),
    'Mixolydian':  colors.HexColor('#fde8d8'),
    'Aeolian':     colors.HexColor('#f8d7da'),
    'Locrian':     colors.HexColor('#e2e3e5'),
}
mode_text_colors = {
    'Ionian':      colors.HexColor('#155724'),
    'Dorian':      colors.HexColor('#0c5460'),
    'Phrygian':    colors.HexColor('#856404'),
    'Lydian':      colors.HexColor('#4a235a'),
    'Mixolydian':  colors.HexColor('#7c3000'),
    'Aeolian':     colors.HexColor('#721c24'),
    'Locrian':     colors.HexColor('#383d41'),
}

def concept_table(rows):
    t = Table(rows, colWidths=[1.6*inch, 5.1*inch])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (0,-1), colors.HexColor('#f0f0f0')),
        ('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'),
        ('FONTNAME', (1,0), (1,-1), 'Helvetica'),
        ('FONTSIZE', (0,0), (-1,-1), 8),
        ('GRID', (0,0), (-1,-1), 0.4, colors.HexColor('#cccccc')),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('TOPPADDING', (0,0), (-1,-1), 3),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    return t

def drill_table(data, col_widths, storage_rows=None):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    s = [
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#2a2a2a')),
        ('TEXTCOLOR', (0,0), (-1,0), colors.white),
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,0), 7),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('FONTNAME', (0,1), (-1,-1), 'Helvetica'),
        ('FONTSIZE', (0,1), (-1,-1), 7),
        ('GRID', (0,0), (-1,-1), 0.4, colors.HexColor('#cccccc')),
        ('TOPPADDING', (0,0), (-1,-1), 1),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]
    for i, row in enumerate(data[1:], 1):
        mode = row[4]
        bg = mode_colors.get(mode, colors.white)
        tc = mode_text_colors.get(mode, colors.black)
        s.append(('BACKGROUND', (0,i), (-1,i), bg))
        s.append(('TEXTCOLOR', (4,i), (4,i), tc))
        s.append(('FONTNAME', (4,i), (4,i), 'Helvetica-Bold'))
    if storage_rows:
        for r in storage_rows:
            s.append(('FONTNAME', (0,r), (-1,r), 'Helvetica-Bold'))
            s.append(('LINEABOVE', (0,r), (-1,r), 1.2, colors.HexColor('#888888')))
            s.append(('LINEBELOW', (0,r), (-1,r), 1.2, colors.HexColor('#888888')))
    t.setStyle(TableStyle(s))
    return t

def legend_table():
    modes_list = ["Ionian", "Dorian", "Phrygian", "Lydian", "Mixolydian", "Aeolian", "Locrian"]
    legend_data = [modes_list]
    leg_w = [0.97*inch] * 7
    t = Table(legend_data, colWidths=leg_w)
    s = [
        ('FONTNAME', (0,0), (-1,-1), 'Helvetica-Bold'),
        ('FONTSIZE', (0,0), (-1,-1), 7.5),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('GRID', (0,0), (-1,-1), 0.4, colors.HexColor('#cccccc')),
    ]
    for i, m in enumerate(modes_list):
        s.append(('BACKGROUND', (i,0), (i,0), mode_colors[m]))
        s.append(('TEXTCOLOR', (i,0), (i,0), mode_text_colors[m]))
    t.setStyle(TableStyle(s))
    return t

story = []

# ════════════════════════════════════════════
# PAGE 1 - LEVER HARP
# ════════════════════════════════════════════

story.append(Paragraph("Lever Harp Practice Drill", title_style))
story.append(Paragraph("Descending Scale Drill  -  Eb Storage Cycle", subtitle_style))
story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#cccccc'), spaceAfter=6))

story.append(concept_table([
    ["Example strings",        "C4  B3  A3  G3  F3  E3  D3  C3  (any descending scale passage)"],
    ["Storage",                "All levers UP  =  Eb major  (all flats)"],
    ["Lever action",           "Engaging lever raises that string one half step"],
    ["Direction",              "One way only  -  levers engage CW order: B  E  A  D  G  C  F"],
    ["Return",                 "Physically flip all levers back up to reset to storage"],
]))
story.append(Spacer(1, 8))

lever_data = [
    ["Step", "Lever", "Key", "Degree", "Mode",       "Feel",         "Tetrachords",       "Notes"],
    ["1",    "storage","Eb", "6",      "Aeolian",    "Dark minor",   "W H W | W H W W",   "Natural minor scale"],
    ["2",    "B dn",  "Bb",  "2",      "Dorian",     "Funky minor",  "W H W | W W H W",   "Minor with raised 6th"],
    ["3",    "E dn",  "F",   "5",      "Mixolydian", "Bluesy",       "W W H | W W H W",   "Major with flat 7th"],
    ["4",    "A dn",  "C",   "1",      "Ionian",     "Bright major", "W W H | W W W H",   "Standard major scale"],
    ["5",    "D dn",  "G",   "4",      "Lydian",     "Ethereal",     "W W W | H W W H",   "Major with raised 4th"],
    ["6",    "G dn",  "D",   "7",      "Locrian",    "Unstable",     "H W W | H W W W",   "Diminished with flat 2nd"],
    ["7",    "C dn",  "A",   "3",      "Phrygian",   "Spanish",      "H W W | W H W W",   "Minor with flat 2nd"],
    ["8",    "F dn",  "E",   "6",      "Aeolian",    "Dark minor",   "W H W | W H W W",   "Natural minor scale"],
    ["9",    "F up",  "A",   "3",      "Phrygian",   "Spanish",      "H W W | W H W W",   "Minor with flat 2nd"],
    ["10",   "C up",  "D",   "7",      "Locrian",    "Unstable",     "H W W | H W W W",   "Diminished with flat 2nd"],
    ["11",   "G up",  "G",   "4",      "Lydian",     "Ethereal",     "W W W | H W W H",   "Major with raised 4th"],
    ["12",   "D up",  "C",   "1",      "Ionian",     "Bright major", "W W H | W W W H",   "Standard major scale"],
    ["13",   "A up",  "F",   "5",      "Mixolydian", "Bluesy",       "W W H | W W H W",   "Major with flat 7th"],
    ["14",   "E up",  "Bb",  "2",      "Dorian",     "Funky minor",  "W H W | W W H W",   "Minor with raised 6th"],
    ["15",   "B up",  "Eb",  "6",      "Aeolian",    "Dark minor",   "W H W | W H W W",   "Natural minor scale"],
    ["16",   "storage","Eb", "6",      "Aeolian",    "Dark minor",   "W H W | W H W W",   "Natural minor scale"],
]

story.append(drill_table(
    lever_data,
    [0.33*inch, 0.55*inch, 0.4*inch, 0.52*inch, 0.82*inch, 0.75*inch, 1.35*inch, 1.6*inch],
    storage_rows=[1, 16]
))
story.append(Spacer(1, 8))

story.append(Paragraph("Key Observations", section_style))
for o in [
    "All 7 modes visited in 8 stops with 7 lever moves",
    "Starts and ends on Aeolian  -  Ionian (C major) at step 4 (midpoint)",
    "Mode sequence: Aeolian > Dorian > Mixolydian > Ionian > Lydian > Locrian > Phrygian > Aeolian",
    "Lever order follows circle of fifths: B  E  A  D  G  C  F",
]:
    story.append(Paragraph("*  " + o, obs_style))

story.append(Spacer(1, 8))
story.append(Paragraph("Mode Color Key", section_style))
story.append(legend_table())

story.append(PageBreak())

# ════════════════════════════════════════════
# PAGE 2 - PEDAL HARP
# ════════════════════════════════════════════

story.append(Paragraph("Pedal Harp Practice Drill", title_style))
story.append(Paragraph("Descending Scale Drill  -  Storage Cycle", subtitle_style))
story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#cccccc'), spaceAfter=6))

story.append(concept_table([
    ["Example strings",        "C4  B3  A3  G3  F3  E3  D3  C3  (any descending scale passage)"],
    ["Storage",                "All pedals UP  =  Cb major  (all flats)"],
    ["Pedal action",           "3 positions: UP (flat)  -  MIDDLE (natural)  -  DOWN (sharp)"],
    ["CW direction (dn)",      "Sharp pedals: F  C  G  D  A  E  B"],
    ["CCW direction (up)",     "Flat pedals:  B  E  A  D  G  C  F"],
]))
story.append(Spacer(1, 8))

pedal_data = [
    ["Step", "Pedal", "Key", "Deg", "Mode",        "Feel",         "Tetrachords",      "Notes"],
    ["1",  "storage","Cb",  "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
    ["2",  "F dn",  "Gb",   "4", "Lydian",      "Ethereal",     "W W W | H W W H",  "Major with raised 4th"],
    ["3",  "C dn",  "Db",   "7", "Locrian",     "Unstable",     "H W W | H W W W",  "Diminished with flat 2nd"],
    ["4",  "G dn",  "Ab",   "3", "Phrygian",    "Spanish",      "H W W | W H W W",  "Minor with flat 2nd"],
    ["5",  "D dn",  "Eb",   "6", "Aeolian",     "Dark minor",   "W H W | W H W W",  "Natural minor scale"],
    ["6",  "A dn",  "Bb",   "2", "Dorian",      "Funky minor",  "W H W | W W H W",  "Minor with raised 6th"],
    ["7",  "E dn",  "F",    "5", "Mixolydian",  "Bluesy",       "W W H | W W H W",  "Major with flat 7th"],
    ["8",  "B dn",  "C",    "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
    ["9",  "F dn",  "G",    "4", "Lydian",      "Ethereal",     "W W W | H W W H",  "Major with raised 4th"],
    ["10", "C dn",  "D",    "7", "Locrian",     "Unstable",     "H W W | H W W W",  "Diminished with flat 2nd"],
    ["11", "G dn",  "A",    "3", "Phrygian",    "Spanish",      "H W W | W H W W",  "Minor with flat 2nd"],
    ["12", "D dn",  "E",    "6", "Aeolian",     "Dark minor",   "W H W | W H W W",  "Natural minor scale"],
    ["13", "A dn",  "B",    "2", "Dorian",      "Funky minor",  "W H W | W W H W",  "Minor with raised 6th"],
    ["14", "E dn",  "F#",   "5", "Mixolydian",  "Bluesy",       "W W H | W W H W",  "Major with flat 7th"],
    ["15", "B dn",  "C#",   "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
    ["16", "B up",  "F#",   "4", "Lydian",      "Ethereal",     "W W W | H W W H",  "Major with raised 4th"],
    ["17", "E up",  "B",    "7", "Locrian",     "Unstable",     "H W W | H W W W",  "Diminished with flat 2nd"],
    ["18", "A up",  "E",    "3", "Phrygian",    "Spanish",      "H W W | W H W W",  "Minor with flat 2nd"],
    ["19", "D up",  "A",    "6", "Aeolian",     "Dark minor",   "W H W | W H W W",  "Natural minor scale"],
    ["20", "G up",  "D",    "2", "Dorian",      "Funky minor",  "W H W | W W H W",  "Minor with raised 6th"],
    ["21", "C up",  "G",    "5", "Mixolydian",  "Bluesy",       "W W H | W W H W",  "Major with flat 7th"],
    ["22", "F up",  "C",    "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
    ["23", "B up",  "F",    "4", "Lydian",      "Ethereal",     "W W W | H W W H",  "Major with raised 4th"],
    ["24", "E up",  "Bb",   "7", "Locrian",     "Unstable",     "H W W | H W W W",  "Diminished with flat 2nd"],
    ["25", "A up",  "Eb",   "3", "Phrygian",    "Spanish",      "H W W | W H W W",  "Minor with flat 2nd"],
    ["26", "D up",  "Ab",   "6", "Aeolian",     "Dark minor",   "W H W | W H W W",  "Natural minor scale"],
    ["27", "G up",  "Db",   "2", "Dorian",      "Funky minor",  "W H W | W W H W",  "Minor with raised 6th"],
    ["28", "C up",  "Gb",   "5", "Mixolydian",  "Bluesy",       "W W H | W W H W",  "Major with flat 7th"],
    ["29", "F up",  "Cb",   "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
    ["30", "storage","Cb",  "1", "Ionian",      "Bright major", "W W H | W W W H",  "Standard major scale"],
]

cw = [0.28*inch, 0.44*inch, 0.33*inch, 0.28*inch, 0.65*inch, 0.62*inch, 1.0*inch, 1.65*inch]
story.append(drill_table(pedal_data, cw, storage_rows=[1, 16, 30]))
story.append(Spacer(1, 8))

story.append(Paragraph("Key Observations", section_style))
for o in [
    "Storage (Cb Ionian) at steps 1 and 30  -  C major (all naturals) at steps 8 and 22",
    "Ionian every 7 steps: 1, 8, 15, 22, 29  -  all 15 keys visited twice (dn then up)",
    "Total: 28 pedal moves  -  30 stops",
]:
    story.append(Paragraph("*  " + o, obs_style))

doc.build(story)
print("Done")
