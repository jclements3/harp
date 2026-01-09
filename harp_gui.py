#!/usr/bin/env python3
"""
HARP_GUI.PY - Tkinter GUI for adjusting harp node parameters

Provides sliders for adjusting plate and neck control points,
with Apply button to regenerate SVG and Save button to persist config.
"""

import json
import os
import subprocess
import tkinter as tk
from tkinter import ttk, messagebox

# Config file path
CONFIG_PATH = os.path.join(os.path.dirname(__file__), "node_config.json")
HARP_SCRIPT = os.path.join(os.path.dirname(__file__), "harp.py")


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


class HarpNodeEditor:
    """GUI for editing harp node parameters."""

    def __init__(self, root):
        self.root = root
        self.root.title("Harp Node Editor")
        self.config = load_config()
        self.sliders = {}

        # Create main frame with scrollbar
        main_frame = ttk.Frame(root, padding="10")
        main_frame.grid(row=0, column=0, sticky="nsew")

        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # Create notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=0, column=0, sticky="nsew", pady=(0, 10))

        # Plate tab
        plate_frame = ttk.Frame(notebook, padding="10")
        notebook.add(plate_frame, text="Plate")
        self._create_plate_controls(plate_frame)

        # Neck tab
        neck_frame = ttk.Frame(notebook, padding="10")
        notebook.add(neck_frame, text="Neck")
        self._create_neck_controls(neck_frame)

        # Button frame
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=1, column=0, sticky="ew")

        ttk.Button(button_frame, text="Apply", command=self.apply_changes).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save", command=self.save_changes).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Reset to Defaults", command=self.reset_defaults).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Open SVG", command=self.open_svg).pack(side=tk.RIGHT, padx=5)

        # Status label
        self.status_var = tk.StringVar(value="Ready")
        status_label = ttk.Label(main_frame, textvariable=self.status_var)
        status_label.grid(row=2, column=0, sticky="w", pady=(10, 0))

    def _create_slider(self, parent, row, label, section, key, min_val, max_val):
        """Create a labeled slider with value display."""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=2)

        # Get current value
        current = self.config[section][key]

        # Variable for slider
        var = tk.DoubleVar(value=current)

        # Value label
        value_label = ttk.Label(parent, text=f"{current:.0f}")
        value_label.grid(row=row, column=2, sticky="w", padx=(5, 0))

        # Slider
        def on_change(val):
            value_label.config(text=f"{float(val):.0f}")

        slider = ttk.Scale(parent, from_=min_val, to=max_val, variable=var,
                          orient=tk.HORIZONTAL, length=200, command=on_change)
        slider.grid(row=row, column=1, sticky="ew", padx=5)

        # Store reference
        self.sliders[(section, key)] = var

        return var

    def _create_plate_controls(self, parent):
        """Create plate parameter controls."""
        row = 0

        ttk.Label(parent, text="P5 (Bottom Middle)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "p5_raise (height above trough)", "plate", "p5_raise", 0, 150)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="P1 (Top Left)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "p1_1_extend", "plate", "p1_1_extend", 0, 150)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="P2 (Top Right)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "p2_2_extend", "plate", "p2_2_extend", 0, 150)
        row += 1
        self._create_slider(parent, row, "p2_1_extend", "plate", "p2_1_extend", 0, 300)
        row += 1
        self._create_slider(parent, row, "p2_1_drop", "plate", "p2_1_drop", 0, 200)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="P4 (Bottom Left)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "p4_x_offset", "plate", "p4_x_offset", -50, 50)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="Bulge (Curve Shape)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "left_bulge", "plate", "left_bulge", -100, 50)
        row += 1
        self._create_slider(parent, row, "right_bulge", "plate", "right_bulge", -50, 100)

    def _create_neck_controls(self, parent):
        """Create neck parameter controls."""
        row = 0

        ttk.Label(parent, text="Neck Margin", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "neck_margin", "neck", "neck_margin", 0, 50)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="N1 (Top Left / C1)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "n1_1_extend", "neck", "n1_1_extend", 0, 150)
        row += 1
        self._create_slider(parent, row, "c1_extend", "neck", "c1_extend", 0, 200)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="N2 (Top Right / G7)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "n2_2_extend", "neck", "n2_2_extend", 0, 150)
        row += 1
        self._create_slider(parent, row, "n2_1_extend", "neck", "n2_1_extend", 0, 300)
        row += 1
        self._create_slider(parent, row, "n2_1_drop", "neck", "n2_1_drop", 0, 200)
        row += 1
        self._create_slider(parent, row, "g7_extend", "neck", "g7_extend", 0, 150)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="N4 (Bottom Left)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "n4_x_offset", "neck", "n4_x_offset", -50, 50)
        row += 1

        ttk.Separator(parent, orient='horizontal').grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="Bulge (Curve Shape)", font=('TkDefaultFont', 10, 'bold')).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 5))
        row += 1
        self._create_slider(parent, row, "neck_left_bulge", "neck", "neck_left_bulge", -100, 50)
        row += 1
        self._create_slider(parent, row, "neck_right_bulge", "neck", "neck_right_bulge", -50, 100)

    def _get_current_values(self):
        """Get current values from all sliders."""
        config = {"plate": {}, "neck": {}}
        for (section, key), var in self.sliders.items():
            config[section][key] = var.get()
        return config

    def apply_changes(self):
        """Apply changes and regenerate SVG."""
        self.status_var.set("Applying changes...")
        self.root.update()

        # Get current values and save to config
        self.config = self._get_current_values()
        save_config(self.config)

        # Run harp.py to regenerate SVG
        try:
            result = subprocess.run(
                ["python3", HARP_SCRIPT, "--skip-validation"],
                capture_output=True, text=True, cwd=os.path.dirname(__file__)
            )
            if result.returncode == 0:
                self.status_var.set("SVG regenerated successfully!")
            else:
                self.status_var.set(f"Error: {result.stderr[:50]}...")
                messagebox.showerror("Error", result.stderr)
        except Exception as e:
            self.status_var.set(f"Error: {str(e)[:50]}...")
            messagebox.showerror("Error", str(e))

    def save_changes(self):
        """Save current values to config file."""
        self.config = self._get_current_values()
        save_config(self.config)
        self.status_var.set(f"Saved to {CONFIG_PATH}")

    def reset_defaults(self):
        """Reset all sliders to default values."""
        defaults = get_defaults()
        for (section, key), var in self.sliders.items():
            var.set(defaults[section][key])
        self.status_var.set("Reset to defaults (click Apply to regenerate)")

    def open_svg(self):
        """Open the generated SVG in default viewer."""
        svg_path = os.path.join(os.path.dirname(__file__), "harp.svg")
        if os.path.exists(svg_path):
            try:
                subprocess.run(["xdg-open", svg_path])
                self.status_var.set("Opened harp.svg")
            except Exception as e:
                self.status_var.set(f"Error opening SVG: {e}")
        else:
            self.status_var.set("SVG not found - click Apply first")


def main():
    root = tk.Tk()
    app = HarpNodeEditor(root)
    root.mainloop()


if __name__ == "__main__":
    main()
