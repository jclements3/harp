# Erard Harp Mechanism - Parts Manufacturing Guide

Based on analysis of Fusion 360 reference models from concertharps.com.

## Overview

The pedal harp mechanism converts pedal motion to disc rotation, shortening strings by one or two semitones. Each string has two disc assemblies (natural and sharp positions).

```
PEDAL → ACTION ROD → BELL CRANK → DISC ROTATION → STRING SHORTENING
```

## Part 1: Disc Assembly (Monolithic)

**Manufacturing:** 316 Stainless Steel, CNC lathe from rod stock

Each disc assembly is a single turned piece:
- Stadium-shaped disc body
- Two prongs on front face (engage string)
- Axle on back face (mounts in plates)

### Specifications by Register

| Register | Strings | Disc Size (LxW) | Thickness | Prong Ø | Prong Length | Axle Ø | Axle Length |
|----------|---------|-----------------|-----------|---------|--------------|--------|-------------|
| Bass     | 1-9     | 56×10mm         | 2.5mm     | 3.0mm   | 23mm         | 4.0mm  | 20mm        |
| Mid      | 10-28   | 48×8mm          | 2.0mm     | 2.5mm   | 20mm         | 3.0mm  | 20mm        |
| Treble   | 29-38   | 18×6mm          | 1.5mm     | 2.0mm   | 6mm          | 2.0mm  | 20mm        |
| High     | 39-47   | 10×4mm          | 1.2mm     | 1.5mm   | 3mm          | 1.5mm  | 20mm        |

**Key dimension:** `prong_length = sphere_radius` (vibration clearance)

**Quantity:** 94 total (47 strings × 2 positions each)

**Generate plates:** `python3 disc_assembly.py`

---

## Part 2: Bell Crank (Erard Style)

**Manufacturing:** Brass, CNC machined or investment cast

Yoke-style lever that converts action rod linear motion to disc rotation.

```
    +--------+
   /    U     \   ← Fork wraps disc axle
  +----+  +----+
       |  |       ← Fork gap
  +----+  +----+
   \          /
    +--[O]--+     ← Pivot hole (center)
       |
       |          ← Arm to action rod
       O          ← Action rod connection hole
```

### Specifications by Register

| Register | Pivot Ø | Arm Length | Fork Length | Fork Gap | Thickness |
|----------|---------|------------|-------------|----------|-----------|
| Bass     | 3.0mm   | 18mm       | 12mm        | 5.0mm    | 4.0mm     |
| Mid      | 3.0mm   | 15mm       | 10mm        | 4.0mm    | 3.5mm     |
| Treble   | 2.5mm   | 12mm       | 8mm         | 3.0mm    | 3.0mm     |
| High     | 2.0mm   | 10mm       | 6mm         | 2.5mm    | 2.5mm     |

**Quantity:** 94 total

**Generate plates:** `python3 bell_crank.py`

---

## Part 3: Front Plate

**Manufacturing:** 316 Stainless Steel, waterjet/laser cut + drilled holes

Mounts disc axles with clearance fit. Holes are positioned at *n and *s locations.

| Spec | Value |
|------|-------|
| Thickness | 6mm |
| Material | 316 Stainless Steel |
| Hole clearance | +0.1mm over axle diameter |
| Total holes | 94 (47 natural + 47 sharp) |

### Hole Sizes

| Register | Axle Ø | Hole Ø |
|----------|--------|--------|
| Bass     | 4.0mm  | 4.1mm  |
| Mid      | 3.0mm  | 3.1mm  |
| Treble   | 2.0mm  | 2.1mm  |
| High     | 1.5mm  | 1.6mm  |

**Generate plate:** `python3 action_plates.py --svg`

---

## Part 4: Back Plate

**Manufacturing:** 316 Stainless Steel, waterjet/laser cut + drilled holes

Provides bearing surface for disc axle ends.

| Spec | Value |
|------|-------|
| Thickness | 4mm |
| Gap from front plate | 9mm |
| Material | 316 Stainless Steel |
| Hole fit | Tight (+0.05mm) |

---

## Part 5: Action Rods

**Manufacturing:** 316 Stainless Steel rod, cut to length, threaded ends

Long rods running from pedal mechanism through pillar to bell cranks.

| Spec | Value |
|------|-------|
| Diameter | 3-4mm |
| Material | 316 Stainless Steel |
| Groups | 7 (one per note: C,D,E,F,G,A,B) |
| End treatment | M3 or M4 thread for adjustment |

---

## Part 6: Pivot Pins

**Manufacturing:** 316 Stainless Steel, turned on lathe

Small pins that bell cranks rotate around.

| Register | Diameter | Length |
|----------|----------|--------|
| Bass     | 3.0mm    | 10mm   |
| Mid      | 3.0mm    | 9mm    |
| Treble   | 2.5mm    | 8mm    |
| High     | 2.0mm    | 7mm    |

---

## Part 7: Return Springs

**Manufacturing:** Purchase standard torsion springs

Return disc to natural position when pedal released.

| Spec | Value |
|------|-------|
| Type | Torsion spring |
| Torque | 5 N·mm |
| Quantity | 94 |

---

## Assembly Stack-up

```
STRING SIDE (z=0)
    ↓
    Prong tips (tangent to vibration sphere)
    ↓
    Disc body (z = -sphere_radius - thickness/2)
    ↓
    Front plate (6mm thick)
    ↓
    Clearance (1mm)
    ↓
    Bell crank zone (8mm)
    ↓
    Clearance (1mm)
    ↓
    Back plate engagement (4mm)
BACK SIDE
```

Total axle length: 6 + 1 + 8 + 1 + 4 = **20mm** (all registers)

---

## File Summary

| File | Purpose |
|------|---------|
| `disc_assembly.py` | Monolithic disc manufacturing plates |
| `bell_crank.py` | Bell crank manufacturing plates |
| `action_plates.py` | Front/back plate specs and visualization |
| `erand.json` | Master configuration with all dimensions |
| `erand.py` | Parametric position calculations |

---

## Reference Videos

Downloaded from concertharps.com (Germán Ocaña):
- `harp_mechanism_fusion360.mp4` - Complete mechanism exploded view
- `harp_mechanism_how_works.mp4` - Lyon & Healy mechanism operation
- `harp_mechanism_erard.mp4` - Erard mechanism (Spanish)
