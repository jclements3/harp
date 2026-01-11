#!/usr/bin/env python3
"""
Half cone geometry test for Erard harp soundbox

The half-cone has:
- Apex at (apex_x, 0, g7f) on the soundboard line
- Edge along soundboard line at 65° from X-axis (in Y=0 plane)
- Curves into -Y, with max depth 250mm at Z=0
"""

import cadquery as cq
import json
import math
from pathlib import Path

# Load string data for g7f reference
with open(Path(__file__).parent / "erard.json") as f:
    data = json.load(f)
strings = data['strings']

# Harp width (G7 x position + pillar offset)
pillar_width = 70
g7 = strings[-1]
harp_width = g7['x_pos'] + pillar_width  # ~770mm

# G7 flat (top of G7 string)
sb_angle = 65
tan_sb = math.tan(math.radians(sb_angle))
g7b = harp_width * tan_sb
g7f = g7b + g7['length']

# Apex on soundboard line at z = g7f
apex_z = g7f
apex_x = g7f / tan_sb  # x = z / tan(65°) to stay on soundboard line
apex_y = 0

# At Z=0, the half cone depth is 250mm
depth_at_z0 = 250  # mm

# Distance from apex to origin along soundboard line
dist_apex_to_origin = math.sqrt(apex_x**2 + apex_z**2)

print(f"Half Cone Parameters:")
print(f"  Apex: ({apex_x:.1f}, {apex_y:.1f}, {apex_z:.1f})")
print(f"  Soundboard angle: {sb_angle}° from X-axis")
print(f"  Distance apex to origin: {dist_apex_to_origin:.1f}mm")
print(f"  Depth at Z=0: {depth_at_z0}mm")

# Build the half-cone by creating a triangular profile in XZ plane at Y=0,
# then for each point on the edge, add a semicircle into -Y

# The edge line goes from apex (apex_x, 0, apex_z) to origin (0, 0, 0)
# This line is at 65° from the X-axis

# Create the shape by lofting semicircular sections along the edge
def make_section(t):
    """
    Create a semicircular section at parameter t along the edge.
    t=0 is apex, t=1 is origin.
    The semicircle is in a plane perpendicular to XZ, curving into -Y.
    """
    # Point on the edge (soundboard line)
    x = apex_x * (1 - t)
    z = apex_z * (1 - t)

    # Depth at this point (proportional to distance from apex)
    depth = depth_at_z0 * t

    if depth < 1:
        return None

    # Semicircle centered at (x, 0, z), curving into -Y
    # The semicircle lies in the plane perpendicular to the soundboard line
    # containing the Y-axis
    n_pts = 24
    points = []

    for i in range(n_pts + 1):
        theta = math.pi * i / n_pts  # 0 to pi
        # Semicircle in a vertical plane perpendicular to XZ at point (x, z)
        # The plane contains the Y-axis direction
        py = -depth * math.sin(theta)  # into -Y
        # The "width" of the semicircle spreads along a direction in XZ plane
        # perpendicular to the soundboard line
        perp_x = math.sin(math.radians(sb_angle))   # sin(65°)
        perp_z = -math.cos(math.radians(sb_angle))  # -cos(65°)
        spread = depth * math.cos(theta)

        px = x + spread * perp_x
        pz = z + spread * perp_z

        points.append((px, py, pz))

    return points

# Create sections
n_sections = 30
loft_wires = []

for i in range(1, n_sections + 1):
    t = i / n_sections
    pts = make_section(t)
    if pts:
        wire = cq.Workplane("XY").polyline(pts).close().val()
        loft_wires.append(wire)

print(f"  Created {len(loft_wires)} sections")

try:
    half_cone = cq.Workplane("XY").add(cq.Solid.makeLoft(loft_wires))
    print("  Loft successful")
except Exception as e:
    print(f"  Loft failed: {e}")

# Add X and Z axes for reference - starting from origin
axis_length = 1000  # mm

# Origin marker - small sphere at (0,0,0)
origin_marker = cq.Workplane("XY").sphere(20)

# X axis - line with arrow/marker at end
x_axis = cq.Workplane("XY").moveTo(0, 0).lineTo(axis_length, 0).val()
x_marker = cq.Workplane("XY").moveTo(axis_length, 0).sphere(30)  # bigger sphere at +X end

# Z axis - line with arrow/marker at end
z_axis = cq.Workplane("XZ").moveTo(0, 0).lineTo(0, axis_length).val()
z_marker = cq.Workplane("XZ").moveTo(0, axis_length).sphere(15)  # smaller sphere at +Z end

# Combine cone with axes and origin marker
assembly = cq.Workplane("XY").add(half_cone.val())
assembly = assembly.add(origin_marker.val())
assembly = assembly.add(x_axis)
assembly = assembly.add(x_marker.val())
assembly = assembly.add(z_axis)
assembly = assembly.add(z_marker.val())

# Print edge endpoints
print(f"  Edge: from apex ({apex_x:.1f}, 0, {apex_z:.1f}) to origin (0, 0, 0)")
print(f"  Edge angle: {sb_angle}° from X-axis")

# Export
cq.exporters.export(half_cone, str(Path(__file__).parent / "half_cone.step"))
print("Exported half_cone.step")

# Export SVG views
from cadquery import exporters

output_dir = Path(__file__).parent

# Side view - X horizontal right, Z vertical up
assembly_rotated = (assembly
    .rotate((0, 0, 0), (0, 1, 0), -90)
    .rotate((0, 0, 0), (1, 0, 0), 180)
    .rotate((0, 0, 0), (0, 1, 0), 180))
exporters.export(assembly_rotated, str(output_dir / "half_cone_side.svg"),
                 opt={"projectionDir": (0, 1, 0), "width": 800, "height": 600})
print("Exported half_cone_side.svg")

# Top view (looking down -Z)
exporters.export(half_cone, str(output_dir / "half_cone_top.svg"),
                 opt={"projectionDir": (0, 0, -1), "width": 800, "height": 400})
print("Exported half_cone_top.svg")

# Front view (looking along -X)
exporters.export(half_cone, str(output_dir / "half_cone_front.svg"),
                 opt={"projectionDir": (-1, 0, 0), "width": 400, "height": 600})
print("Exported half_cone_front.svg")

print("\nDone!")
