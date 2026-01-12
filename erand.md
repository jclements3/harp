# Erard Concert Harp Design

**1901 Paris - 47 String Concert Pedal Harp**

Source: https://harpcanada.com/harpmaking/erard.htm

## Coordinate System

```
     Z (height)
     |
     |    / soundboard plane (65deg from X)
     |   /
     |  /
     | /
     +------------ X (string positions)
    /
   /
  Y (depth into soundbox, negative)
```

- **X-axis**: String positions (C1=0mm, G7=699.8mm)
- **Y-axis**: Depth into soundbox (negative direction)
- **Z-axis**: Height from floor (soundboard base to tuning pins)

## String Point Terminology

Each string has 5 key points along its path from top to bottom:

| Suffix | Point | Description |
|--------|-------|-------------|
| **p** | peg | tuning peg circle tangent to string (top) |
| **f** | flat | flat pin circle tangent to string |
| **n** | natural | natural disc prong (3 o'clock) tangent |
| **s** | sharp | sharp disc prong (3 o'clock) tangent |
| **b** | base | soundboard contact point (bottom) |

String path: `p → f → n → s → b`

Examples:
- `c1p` - C1 peg point (longest string, top)
- `c1b` - C1 base point (longest string, soundboard)
- `g7b` - G7 base point (shortest string, soundboard)
- `f3n` - F3 natural disc contact

## Dimensions

| Parameter | Value | Notes |
|-----------|-------|-------|
| Number of strings | 47 | C1 to G7 |
| Soundboard angle | 65° | from X-axis |
| Soundboard length | 1294.7 mm | |
| String spacing (total) | 717.0 mm | C1 to G7 |
| Finger gap | 14.0 mm | edge-to-edge |
| Longest string (C1) | 1514.9 mm | |
| Shortest string (G7) | 60.6 mm | |
| Pillar offset | 70 mm | beyond G7 |
| Harp width | 769.8 mm | |
| Soundbox depth | 250 mm | at origin |

## String Layout

| Range | Strings | Material | Winding |
|-------|---------|----------|---------|
| Bass | 1-9 (C1-D2) | Steel core | Bronze wound |
| Low-mid | 10-20 (E2-A3) | Nylon core | Nylon wound |
| Treble | 21-47 (B3-G7) | Nylon | None |

### Key Strings

| String | Note | Frequency | Length | Tension |
|--------|------|-----------|--------|---------|
| 1 | C1 | 32.7 Hz | 1514.9 mm | 234.4 N |
| 8 | C2 | 65.4 Hz | 1318.0 mm | 222.1 N |
| 15 | C3 | 130.8 Hz | 989.8 mm | 210.4 N |
| 22 | C4 (middle C) | 261.6 Hz | 570.6 mm | 120.5 N |
| 29 | C5 | 523.3 Hz | 328.2 mm | 102.1 N |
| 36 | C6 | 1046.5 Hz | 191.9 mm | 89.3 N |
| 43 | C7 | 2093.0 Hz | 101.0 mm | 60.4 N |
| 47 | G7 | 3136.0 Hz | 60.6 mm | 48.8 N |

## Soundbox Geometry

The soundbox is modeled as a **half-cone**:

- Apex located on the soundboard line at the G7 position
- Edge runs along the soundboard plane (65° from X-axis)
- Curved surface extends into -Y direction only
- Maximum depth of 250mm at the origin (bass end)
- Cross-section is semi-elliptical

### Half-Cone Parameters

```
Apex position: x = 769.8mm, y = 0, z = 1651mm (on soundboard line)
Cone axis: along soundboard line toward origin
Base radius: 250mm at Z=0
Angular extent: 180° (half-cone, -Y direction only)
```

## Design Calculations

The design follows classical harp proportions:

1. **String lengths**: Calculated from `L = sqrt(T/μ) / (2f)` where:
   - T = tension (N)
   - μ = linear mass density (kg/m)
   - f = frequency (Hz)

2. **Soundboard angle**: 65° provides optimal plucking angle and string break

3. **String spacing**: 14mm edge-to-edge gap allows comfortable finger placement

4. **Tension curve**: Relatively flat across the range (48-236 N) for consistent feel

## Files

| File | Description |
|------|-------------|
| `erand.json` | Machine-readable specification with all string data |
| `erand.py` | Multi-view SVG generator |
| `erand.svg` | Generated orthographic views |
| `erand.dxf` | Original DXF drawing reference |
| `erand.md` | This documentation |

## Generating the SVG

```bash
cd design
python3 erand.py
```

Output includes:
- **Side view (XZ)**: Shows all 47 strings, soundboard line, and neck curve
- **Front view (YZ)**: Shows soundbox depth and representative strings
- **Top view (XY)**: Shows string positions and soundbox outline
- **String table**: Key string data (one per octave)

## Color Coding

In the SVG views:
- **Red**: X-axis
- **Green**: Y-axis
- **Blue**: Z-axis
- **Brown**: Soundboard line
- **Gray**: Soundbox outline
- **Dark brown**: Bass strings (wound)
- **Dark gray**: Treble strings (plain nylon)
