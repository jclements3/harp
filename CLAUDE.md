# CLAUDE.md

This file provides guidance to Claude Code when working with code in this repository.

## Project Overview

CLEMENTS47 is a parametric 47-string pedal harp design using CadQuery for 3D CAD modeling. The design computes string positions, hardware dimensions, and generates manufacturable 3D models.

## Common Commands

```bash
# Generate harp SVG and BOM
python3 harp.py

# With physics analysis
python3 harp.py --physics

# BOM only (no SVG)
python3 harp.py --bom-only
```

## Architecture

### Core Modules

- `harp.py` - Main entry point, CLI argument handling
- `harp_models.py` - Data classes (Disc, FlatPin, TuningPin, String, Neck, Harp)
- `harp_geometry.py` - Geometry calculations, loads config from `harp.json`
- `harp_physics.py` - String tension, forces, vibration analysis
- `harp_bom.py` - Bill of materials generation
- `harp_renderer.py` - SVG rendering
- `harp_validation.py` - Constraint validation

### Configuration

`harp.json` - Simplified config with geometry and materials. Derived values (frequency, tension, diameter) computed from physics.

## Hardware Specs

### Flat Pins (on stainless steel plates)
| Range | Thread | Shaft |
|-------|--------|-------|
| Bass (1-9) | M10 | 8mm |
| Mid (10-28) | M8 | 6mm |
| Treble (29-47) | M6 | 5mm |

### Tuning Pins (on neck)
| Range | Thread | Diameter |
|-------|--------|----------|
| Bass (1-9) | M6 | 6mm |
| Mid (10-28) | M5 | 5mm |
| Treble (29-47) | M4 | 4mm |

### Discs (2 per string)
| Range | Prong | Thickness |
|-------|-------|-----------|
| Bass (1-9) | 5mm | 5mm |
| Mid (10-28) | 4mm | 4mm |
| Treble (29-47) | 3mm | 3mm |

## Key Constants

- Finger gap: 14mm perpendicular
- String tilt: -3° from vertical
- Neck clearance: 30mm above highest flat pin
- Natural disc rotation: 45°
- Sharp disc rotation: 90°
- Semitone ratio: 2^(1/12) ≈ 1.0595
