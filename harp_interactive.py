#!/usr/bin/env python3
"""
HARP_INTERACTIVE.PY - Interactive node/control point editor using matplotlib

Drag nodes and control points directly on the visualization.
Curves update in real-time as you drag.
"""

import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, PathPatch
from matplotlib.path import Path
from matplotlib.widgets import Button
from matplotlib.lines import Line2D

# Config file path
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "node_config.json")


def load_config():
    """Load config from JSON file."""
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, 'r') as f:
            return json.load(f)
    return get_defaults()


def get_defaults():
    """Return default config values."""
    return {
        "plate": {
            "p5_raise": 70,
            "left_bulge": -25,
            "right_bulge": 25,
            "p1_1_extend": 50,
            "p2_2_extend": 50,
            "p2_1_extend": 150,
            "p2_1_drop": 80,
            "p4_x_offset": -15
        },
        "neck": {
            "neck_margin": 20,
            "n4_x_offset": -20,
            "neck_left_bulge": -20,
            "neck_right_bulge": 25,
            "n1_1_extend": 50,
            "n2_2_extend": 50,
            "n2_1_extend": 150,
            "n2_1_drop": 80,
            "c1_extend": 100,
            "g7_extend": 50
        }
    }


def save_config(config):
    """Save config to JSON file."""
    with open(CONFIG_PATH, 'w') as f:
        json.dump(config, f, indent=2)


class DraggablePoint:
    """A draggable point on the plot."""

    def __init__(self, ax, x, y, name, color='red', size=10, is_control=False, parent=None):
        self.ax = ax
        self.name = name
        self.is_control = is_control
        self.parent = parent  # Parent node for control points
        self.x = x
        self.y = y

        marker = 'o' if not is_control else 's'
        self.point, = ax.plot(x, y, marker, color=color, markersize=size,
                              picker=5, zorder=10)
        self.label = ax.annotate(name, (x, y), textcoords="offset points",
                                  xytext=(5, 5), fontsize=8, zorder=11)

        # Line to parent (for control points)
        self.control_line = None
        if parent:
            self.control_line, = ax.plot([parent.x, x], [parent.y, y],
                                         '--', color='gray', linewidth=0.5, zorder=5)

    def update_position(self, x, y):
        """Update the point position."""
        self.x = x
        self.y = y
        self.point.set_data([x], [y])
        self.label.set_position((x, y))

        if self.control_line and self.parent:
            self.control_line.set_data([self.parent.x, x], [self.parent.y, y])

    def update_control_line(self):
        """Update control line to parent."""
        if self.control_line and self.parent:
            self.control_line.set_data([self.parent.x, self.x], [self.parent.y, self.y])


class InteractiveHarpEditor:
    """Interactive editor for harp plate and neck curves."""

    def __init__(self):
        self.config = load_config()
        self.fig, self.ax = plt.subplots(1, 1, figsize=(14, 10))
        self.ax.set_aspect('equal')
        self.ax.set_title('Harp Node Editor - Drag nodes and control points')
        self.ax.grid(True, alpha=0.3)

        # Store draggable points
        self.points = {}
        self.selected_point = None
        self.press_data = None

        # Curve patches
        self.plate_patch = None
        self.neck_patch = None
        self.hardware_drawn = False  # Flag to draw hardware only once

        # Load harp geometry to get base positions
        self._load_harp_geometry()

        # Create points
        self._create_points()

        # Draw initial curves
        self._update_curves()

        # Set axis limits
        self._set_axis_limits()

        # Connect events
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)

        # Add buttons
        self._add_buttons()

    def _load_harp_geometry(self):
        """Load harp geometry from harp.json to get base positions."""
        # Import harp modules
        import sys
        sys.path.insert(0, os.path.dirname(__file__))
        from harp_geometry import load_harp_from_json

        harp_json = os.path.join(os.path.dirname(__file__), "harp.json")
        self.harp = load_harp_from_json(harp_json)

        # Get tuning pin positions for neckpeg
        c1_tp = self.harp.strings[0].tuning_pin
        g7_tp = self.harp.strings[-1].tuning_pin
        peg_r = c1_tp.diameter_mm / 2

        self.neckpeg_c1 = (c1_tp.x_mm, c1_tp.z_mm - peg_r)
        self.neckpeg_g7 = (g7_tp.x_mm, g7_tp.z_mm - peg_r)

        # Neckpeg direction
        dx = self.neckpeg_g7[0] - self.neckpeg_c1[0]
        dz = self.neckpeg_g7[1] - self.neckpeg_c1[1]
        length = math.sqrt(dx**2 + dz**2)
        self.neckpeg_ux = dx / length
        self.neckpeg_uz = dz / length

        # Get flat pin positions for plate
        c1_fp = self.harp.strings[0].flat_pin
        g7_fp = self.harp.strings[-1].flat_pin

        # Calculate disc positions for contour
        self.disc_positions = []
        for s in self.harp.strings:
            self.disc_positions.append((s.natural_disc.x_mm, s.natural_disc.z_mm))
            self.disc_positions.append((s.sharp_disc.x_mm, s.sharp_disc.z_mm))

        # Store tuning pin positions
        self.tuning_pins = []
        for s in self.harp.strings:
            tp = s.tuning_pin
            self.tuning_pins.append((tp.x_mm, tp.z_mm, tp.diameter_mm))

        # Store flat pin positions
        self.flat_pins = []
        for s in self.harp.strings:
            fp = s.flat_pin
            self.flat_pins.append((fp.x_mm, fp.z_mm, fp.shaft_diameter_mm))

        # Base positions for plate corners
        plate_margin = 15
        self.p1_base = (c1_fp.x_mm, max(d[1] for d in self.disc_positions) + plate_margin)
        self.p2_base = (g7_fp.x_mm + plate_margin,
                        max(d[1] for d in self.disc_positions if d[0] > g7_fp.x_mm - 100) + plate_margin)

        # Bottom trough
        self.bot_trough_z = min(d[1] for d in self.disc_positions) - plate_margin

    def _create_points(self):
        """Create all draggable points."""
        pc = self.config["plate"]
        nc = self.config["neck"]

        # Calculate actual positions from config
        # Plate nodes
        p1 = self.p1_base
        p2 = self.p2_base
        p4 = (p1[0] + pc["p4_x_offset"], self.bot_trough_z)
        p5_x = (p2[0] + p4[0]) / 2
        p5 = (p5_x, self.bot_trough_z + pc["p5_raise"])

        # Create plate nodes
        self.points['p1'] = DraggablePoint(self.ax, p1[0], p1[1], 'p1', 'red', 12)
        self.points['p2'] = DraggablePoint(self.ax, p2[0], p2[1], 'p2', 'red', 12)
        self.points['p4'] = DraggablePoint(self.ax, p4[0], p4[1], 'p4', 'red', 12)
        self.points['p5'] = DraggablePoint(self.ax, p5[0], p5[1], 'p5', 'red', 12)

        # Plate control points
        # p1-1: leaving p1 along neckpeg
        p1_1 = (p1[0] + self.neckpeg_ux * pc["p1_1_extend"],
                p1[1] + self.neckpeg_uz * pc["p1_1_extend"])
        self.points['p1-1'] = DraggablePoint(self.ax, p1_1[0], p1_1[1], 'p1-1',
                                              'pink', 8, is_control=True, parent=self.points['p1'])

        # p2 control points (need direction calculation)
        p2_1_dir_x = self.neckpeg_ux * pc["p2_1_extend"]
        p2_1_dir_z = self.neckpeg_uz * pc["p2_1_extend"] - pc["p2_1_drop"]
        p2_1_dir_len = math.sqrt(p2_1_dir_x**2 + p2_1_dir_z**2)
        p2_1_dir_ux = p2_1_dir_x / p2_1_dir_len
        p2_1_dir_uz = p2_1_dir_z / p2_1_dir_len

        p2_2 = (p2[0] - p2_1_dir_ux * pc["p2_2_extend"],
                p2[1] - p2_1_dir_uz * pc["p2_2_extend"])
        p2_1 = (p2[0] + p2_1_dir_ux * pc["p2_1_extend"],
                p2[1] + p2_1_dir_uz * pc["p2_1_extend"])

        self.points['p2-2'] = DraggablePoint(self.ax, p2_2[0], p2_2[1], 'p2-2',
                                              'pink', 8, is_control=True, parent=self.points['p2'])
        self.points['p2-1'] = DraggablePoint(self.ax, p2_1[0], p2_1[1], 'p2-1',
                                              'pink', 8, is_control=True, parent=self.points['p2'])

        # p5 control points (horizontal)
        p5_2 = ((p2[0] + p5[0]) / 2, p5[1])
        p5_1 = ((p4[0] + p5[0]) / 2, p5[1])
        self.points['p5-2'] = DraggablePoint(self.ax, p5_2[0], p5_2[1], 'p5-2',
                                              'pink', 8, is_control=True, parent=self.points['p5'])
        self.points['p5-1'] = DraggablePoint(self.ax, p5_1[0], p5_1[1], 'p5-1',
                                              'pink', 8, is_control=True, parent=self.points['p5'])

        # p4 control points
        p4_2 = ((p5[0] + p4[0]) / 2, p4[1])
        p4_1 = (p4[0] + pc["left_bulge"], p4[1])
        self.points['p4-2'] = DraggablePoint(self.ax, p4_2[0], p4_2[1], 'p4-2',
                                              'pink', 8, is_control=True, parent=self.points['p4'])
        self.points['p4-1'] = DraggablePoint(self.ax, p4_1[0], p4_1[1], 'p4-1',
                                              'pink', 8, is_control=True, parent=self.points['p4'])

        # p1-2 (arriving)
        p1_2 = (p1[0] + pc["left_bulge"], p1[1])
        self.points['p1-2'] = DraggablePoint(self.ax, p1_2[0], p1_2[1], 'p1-2',
                                              'pink', 8, is_control=True, parent=self.points['p1'])

        # Neck nodes
        n1 = self.neckpeg_c1
        n2 = self.neckpeg_g7
        n4 = (p4[0] - nc["neck_margin"] + nc["n4_x_offset"], p4[1] - nc["neck_margin"])
        n5 = (p5[0], p5[1] - nc["neck_margin"])

        self.points['n1'] = DraggablePoint(self.ax, n1[0], n1[1], 'n1', 'blue', 12)
        self.points['n2'] = DraggablePoint(self.ax, n2[0], n2[1], 'n2', 'blue', 12)
        self.points['n4'] = DraggablePoint(self.ax, n4[0], n4[1], 'n4', 'blue', 12)
        self.points['n5'] = DraggablePoint(self.ax, n5[0], n5[1], 'n5', 'blue', 12)

        # Neck control points
        n1_1 = (n1[0] + self.neckpeg_ux * nc["n1_1_extend"],
                n1[1] + self.neckpeg_uz * nc["n1_1_extend"])
        self.points['n1-1'] = DraggablePoint(self.ax, n1_1[0], n1_1[1], 'n1-1',
                                              'cyan', 8, is_control=True, parent=self.points['n1'])

        # n2 control points
        n2_1_dir_x = self.neckpeg_ux * nc["n2_1_extend"]
        n2_1_dir_z = self.neckpeg_uz * nc["n2_1_extend"] - nc["n2_1_drop"]
        n2_1_dir_len = math.sqrt(n2_1_dir_x**2 + n2_1_dir_z**2)
        n2_1_dir_ux = n2_1_dir_x / n2_1_dir_len
        n2_1_dir_uz = n2_1_dir_z / n2_1_dir_len

        n2_2 = (n2[0] - n2_1_dir_ux * nc["n2_2_extend"],
                n2[1] - n2_1_dir_uz * nc["n2_2_extend"])
        n2_1 = (n2[0] + n2_1_dir_ux * nc["n2_1_extend"],
                n2[1] + n2_1_dir_uz * nc["n2_1_extend"])

        self.points['n2-2'] = DraggablePoint(self.ax, n2_2[0], n2_2[1], 'n2-2',
                                              'cyan', 8, is_control=True, parent=self.points['n2'])
        self.points['n2-1'] = DraggablePoint(self.ax, n2_1[0], n2_1[1], 'n2-1',
                                              'cyan', 8, is_control=True, parent=self.points['n2'])

        # n5 control points
        n5_2 = ((n2[0] + n5[0]) / 2, n5[1])
        n5_1 = ((n4[0] + n5[0]) / 2, n5[1])
        self.points['n5-2'] = DraggablePoint(self.ax, n5_2[0], n5_2[1], 'n5-2',
                                              'cyan', 8, is_control=True, parent=self.points['n5'])
        self.points['n5-1'] = DraggablePoint(self.ax, n5_1[0], n5_1[1], 'n5-1',
                                              'cyan', 8, is_control=True, parent=self.points['n5'])

        # n4 control points
        n4_2 = ((n5[0] + n4[0]) / 2, n4[1])
        n4_1 = (n4[0] + nc["neck_left_bulge"], (n1[1] + n4[1]) / 2)
        self.points['n4-2'] = DraggablePoint(self.ax, n4_2[0], n4_2[1], 'n4-2',
                                              'cyan', 8, is_control=True, parent=self.points['n4'])
        self.points['n4-1'] = DraggablePoint(self.ax, n4_1[0], n4_1[1], 'n4-1',
                                              'cyan', 8, is_control=True, parent=self.points['n4'])

        # n1-2 (arriving)
        n1_2 = (n1[0] - self.neckpeg_ux * nc["c1_extend"],
                n1[1] - self.neckpeg_uz * nc["c1_extend"])
        self.points['n1-2'] = DraggablePoint(self.ax, n1_2[0], n1_2[1], 'n1-2',
                                              'cyan', 8, is_control=True, parent=self.points['n1'])

    def _bezier_curve(self, p0, p1, p2, p3, num_points=50):
        """Generate points along a cubic bezier curve."""
        t = np.linspace(0, 1, num_points)
        x = (1-t)**3 * p0[0] + 3*(1-t)**2*t * p1[0] + 3*(1-t)*t**2 * p2[0] + t**3 * p3[0]
        y = (1-t)**3 * p0[1] + 3*(1-t)**2*t * p1[1] + 3*(1-t)*t**2 * p2[1] + t**3 * p3[1]
        return x, y

    def _update_curves(self):
        """Update the plate and neck curves based on current point positions."""
        # Remove old curves
        if self.plate_patch:
            self.plate_patch.remove()
        if self.neck_patch:
            self.neck_patch.remove()

        # Get point positions
        p1 = (self.points['p1'].x, self.points['p1'].y)
        p2 = (self.points['p2'].x, self.points['p2'].y)
        p4 = (self.points['p4'].x, self.points['p4'].y)
        p5 = (self.points['p5'].x, self.points['p5'].y)

        p1_1 = (self.points['p1-1'].x, self.points['p1-1'].y)
        p1_2 = (self.points['p1-2'].x, self.points['p1-2'].y)
        p2_1 = (self.points['p2-1'].x, self.points['p2-1'].y)
        p2_2 = (self.points['p2-2'].x, self.points['p2-2'].y)
        p4_1 = (self.points['p4-1'].x, self.points['p4-1'].y)
        p4_2 = (self.points['p4-2'].x, self.points['p4-2'].y)
        p5_1 = (self.points['p5-1'].x, self.points['p5-1'].y)
        p5_2 = (self.points['p5-2'].x, self.points['p5-2'].y)

        # Build plate path
        plate_x, plate_y = [], []

        # p1 -> p2
        x, y = self._bezier_curve(p1, p1_1, p2_2, p2)
        plate_x.extend(x)
        plate_y.extend(y)

        # p2 -> p5
        x, y = self._bezier_curve(p2, p2_1, p5_2, p5)
        plate_x.extend(x)
        plate_y.extend(y)

        # p5 -> p4
        x, y = self._bezier_curve(p5, p5_1, p4_2, p4)
        plate_x.extend(x)
        plate_y.extend(y)

        # p4 -> p1
        x, y = self._bezier_curve(p4, p4_1, p1_2, p1)
        plate_x.extend(x)
        plate_y.extend(y)

        self.plate_patch, = self.ax.plot(plate_x, plate_y, '-', color='gold',
                                          linewidth=2, zorder=1, label='Plate')

        # Neck path
        n1 = (self.points['n1'].x, self.points['n1'].y)
        n2 = (self.points['n2'].x, self.points['n2'].y)
        n4 = (self.points['n4'].x, self.points['n4'].y)
        n5 = (self.points['n5'].x, self.points['n5'].y)

        n1_1 = (self.points['n1-1'].x, self.points['n1-1'].y)
        n1_2 = (self.points['n1-2'].x, self.points['n1-2'].y)
        n2_1 = (self.points['n2-1'].x, self.points['n2-1'].y)
        n2_2 = (self.points['n2-2'].x, self.points['n2-2'].y)
        n4_1 = (self.points['n4-1'].x, self.points['n4-1'].y)
        n4_2 = (self.points['n4-2'].x, self.points['n4-2'].y)
        n5_1 = (self.points['n5-1'].x, self.points['n5-1'].y)
        n5_2 = (self.points['n5-2'].x, self.points['n5-2'].y)

        neck_x, neck_y = [], []

        # n1 -> n2
        x, y = self._bezier_curve(n1, n1_1, n2_2, n2)
        neck_x.extend(x)
        neck_y.extend(y)

        # n2 -> n5
        x, y = self._bezier_curve(n2, n2_1, n5_2, n5)
        neck_x.extend(x)
        neck_y.extend(y)

        # n5 -> n4
        x, y = self._bezier_curve(n5, n5_1, n4_2, n4)
        neck_x.extend(x)
        neck_y.extend(y)

        # n4 -> n1
        x, y = self._bezier_curve(n4, n4_1, n1_2, n1)
        neck_x.extend(x)
        neck_y.extend(y)

        self.neck_patch, = self.ax.plot(neck_x, neck_y, '-', color='black',
                                         linewidth=2, zorder=2, label='Neck')

        # Draw hardware only once (discs, tuning pins, flat pins)
        if not self.hardware_drawn:
            # Draw discs as reference (small circles)
            for dx, dz in self.disc_positions:
                circle = Circle((dx, dz), 8, fill=True, facecolor='#E8E8E8',
                               edgecolor='gray', alpha=0.5, zorder=0)
                self.ax.add_patch(circle)

            # Draw tuning pins (circles at top)
            for tx, tz, diam in self.tuning_pins:
                circle = Circle((tx, tz), diam/2, fill=True, facecolor='#4488CC',
                               edgecolor='#224466', alpha=0.7, zorder=0)
                self.ax.add_patch(circle)

            # Draw flat pins (smaller circles)
            for fx, fz, diam in self.flat_pins:
                circle = Circle((fx, fz), diam/2, fill=True, facecolor='#88CC44',
                               edgecolor='#446622', alpha=0.7, zorder=0)
                self.ax.add_patch(circle)

            self.hardware_drawn = True

        self.fig.canvas.draw_idle()

    def _set_axis_limits(self):
        """Set axis limits based on point positions."""
        all_x = [p.x for p in self.points.values()]
        all_y = [p.y for p in self.points.values()]

        margin = 50
        self.ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
        self.ax.set_ylim(min(all_y) - margin, max(all_y) + margin)

    def _on_press(self, event):
        """Handle mouse press."""
        if event.inaxes != self.ax:
            return

        # Find closest point
        min_dist = float('inf')
        closest = None

        for name, point in self.points.items():
            dist = math.sqrt((point.x - event.xdata)**2 + (point.y - event.ydata)**2)
            if dist < min_dist and dist < 20:  # 20 pixel threshold
                min_dist = dist
                closest = name

        if closest:
            self.selected_point = closest
            self.press_data = (event.xdata, event.ydata,
                              self.points[closest].x, self.points[closest].y)

    def _on_release(self, event):
        """Handle mouse release."""
        self.selected_point = None
        self.press_data = None

    def _on_motion(self, event):
        """Handle mouse motion."""
        if self.selected_point is None or event.inaxes != self.ax:
            return

        # Calculate new position
        dx = event.xdata - self.press_data[0]
        dy = event.ydata - self.press_data[1]
        new_x = self.press_data[2] + dx
        new_y = self.press_data[3] + dy

        # Update point position
        self.points[self.selected_point].update_position(new_x, new_y)

        # Update control lines for control points attached to this node
        for name, point in self.points.items():
            if point.parent and point.parent.name == self.selected_point:
                point.update_control_line()

        # Update curves
        self._update_curves()

    def _add_buttons(self):
        """Add Save and Reset buttons."""
        # Save button - store as instance attribute to prevent garbage collection
        self.ax_save = plt.axes([0.7, 0.02, 0.1, 0.04])
        self.btn_save = Button(self.ax_save, 'Save')
        self.btn_save.on_clicked(self._save_config)

        # Reset button
        self.ax_reset = plt.axes([0.81, 0.02, 0.1, 0.04])
        self.btn_reset = Button(self.ax_reset, 'Reset')
        self.btn_reset.on_clicked(self._reset_config)

        # Print button
        self.ax_print = plt.axes([0.59, 0.02, 0.1, 0.04])
        self.btn_print = Button(self.ax_print, 'Print Config')
        self.btn_print.on_clicked(self._print_config)

    def _calculate_config_from_points(self):
        """Calculate config values from current point positions."""
        p1 = (self.points['p1'].x, self.points['p1'].y)
        p2 = (self.points['p2'].x, self.points['p2'].y)
        p4 = (self.points['p4'].x, self.points['p4'].y)
        p5 = (self.points['p5'].x, self.points['p5'].y)

        p1_1 = (self.points['p1-1'].x, self.points['p1-1'].y)
        p2_1 = (self.points['p2-1'].x, self.points['p2-1'].y)
        p2_2 = (self.points['p2-2'].x, self.points['p2-2'].y)
        p4_1 = (self.points['p4-1'].x, self.points['p4-1'].y)

        n1 = (self.points['n1'].x, self.points['n1'].y)
        n2 = (self.points['n2'].x, self.points['n2'].y)
        n4 = (self.points['n4'].x, self.points['n4'].y)
        n5 = (self.points['n5'].x, self.points['n5'].y)

        n1_1 = (self.points['n1-1'].x, self.points['n1-1'].y)
        n1_2 = (self.points['n1-2'].x, self.points['n1-2'].y)
        n2_1 = (self.points['n2-1'].x, self.points['n2-1'].y)
        n2_2 = (self.points['n2-2'].x, self.points['n2-2'].y)
        n4_1 = (self.points['n4-1'].x, self.points['n4-1'].y)

        # Calculate plate config
        p5_raise = p5[1] - self.bot_trough_z
        p4_x_offset = p4[0] - self.p1_base[0]
        left_bulge = p4_1[0] - p4[0]

        # p1_1_extend: distance along neckpeg from p1
        p1_1_dx = p1_1[0] - p1[0]
        p1_1_dz = p1_1[1] - p1[1]
        p1_1_extend = math.sqrt(p1_1_dx**2 + p1_1_dz**2)

        # p2 control points
        p2_1_dx = p2_1[0] - p2[0]
        p2_1_dz = p2_1[1] - p2[1]
        p2_1_extend = math.sqrt(p2_1_dx**2 + p2_1_dz**2)

        p2_2_dx = p2[0] - p2_2[0]
        p2_2_dz = p2[1] - p2_2[1]
        p2_2_extend = math.sqrt(p2_2_dx**2 + p2_2_dz**2)

        # Calculate drop from neckpeg line
        neckpeg_z_at_p2_1 = p2[1] + self.neckpeg_uz * p2_1_extend
        p2_1_drop = neckpeg_z_at_p2_1 - p2_1[1]

        # Neck config
        nc = self.config["neck"]
        neck_margin = p5[1] - n5[1]
        n4_x_offset = n4[0] - (p4[0] - neck_margin)
        neck_left_bulge = n4_1[0] - n4[0]

        n1_1_dx = n1_1[0] - n1[0]
        n1_1_dz = n1_1[1] - n1[1]
        n1_1_extend = math.sqrt(n1_1_dx**2 + n1_1_dz**2)

        n2_1_dx = n2_1[0] - n2[0]
        n2_1_dz = n2_1[1] - n2[1]
        n2_1_extend = math.sqrt(n2_1_dx**2 + n2_1_dz**2)

        n2_2_dx = n2[0] - n2_2[0]
        n2_2_dz = n2[1] - n2_2[1]
        n2_2_extend = math.sqrt(n2_2_dx**2 + n2_2_dz**2)

        neckpeg_z_at_n2_1 = n2[1] + self.neckpeg_uz * n2_1_extend
        n2_1_drop = neckpeg_z_at_n2_1 - n2_1[1]

        c1_dx = n1[0] - n1_2[0]
        c1_dz = n1[1] - n1_2[1]
        c1_extend = math.sqrt(c1_dx**2 + c1_dz**2)

        return {
            "plate": {
                "p5_raise": round(p5_raise, 1),
                "left_bulge": round(left_bulge, 1),
                "right_bulge": self.config["plate"]["right_bulge"],
                "p1_1_extend": round(p1_1_extend, 1),
                "p2_2_extend": round(p2_2_extend, 1),
                "p2_1_extend": round(p2_1_extend, 1),
                "p2_1_drop": round(p2_1_drop, 1),
                "p4_x_offset": round(p4_x_offset, 1)
            },
            "neck": {
                "neck_margin": round(neck_margin, 1),
                "n4_x_offset": round(n4_x_offset, 1),
                "neck_left_bulge": round(neck_left_bulge, 1),
                "neck_right_bulge": nc["neck_right_bulge"],
                "n1_1_extend": round(n1_1_extend, 1),
                "n2_2_extend": round(n2_2_extend, 1),
                "n2_1_extend": round(n2_1_extend, 1),
                "n2_1_drop": round(n2_1_drop, 1),
                "c1_extend": round(c1_extend, 1),
                "g7_extend": nc["g7_extend"]
            }
        }

    def _save_config(self, event):
        """Save current positions to config."""
        config = self._calculate_config_from_points()
        save_config(config)
        print("Config saved to", CONFIG_PATH)
        print(json.dumps(config, indent=2))

    def _reset_config(self, event):
        """Reset to default config."""
        self.config = get_defaults()
        save_config(self.config)
        # Recreate points
        for point in list(self.points.values()):
            point.point.remove()
            point.label.remove()
            if point.control_line:
                point.control_line.remove()
        self.points.clear()
        self._create_points()
        self._update_curves()
        print("Reset to defaults")

    def _print_config(self, event):
        """Print current config values."""
        config = self._calculate_config_from_points()
        print("\nCurrent config:")
        print(json.dumps(config, indent=2))

    def show(self):
        """Show the interactive editor."""
        plt.legend()
        plt.show()


def main():
    editor = InteractiveHarpEditor()
    editor.show()


if __name__ == "__main__":
    main()
