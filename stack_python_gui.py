import tkinter as tk
from tkinter import filedialog, messagebox # https://docs.python.org/3/library/tkinter.html

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backend_bases import key_press_handler # https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_tk_sgskip.html
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import re
import mplcursors

dataset = None
start_idx = 6 # hardcoded for csv file
bar_chart_frame = None  # Global reference to bar chart frame
pie_chart_container = None  # Global reference to pie chart frame
lipid_checkboxes = {}  # Dict to store {lipid_name: BooleanVar}
plot_button = None  # Global reference to plot button

def load_csv():
    global dataset
    file_path = filedialog.askopenfilename(title="Select CSV files", filetypes=[("CSV files", "*.csv")])

    if not file_path:
        return

    dataset = pd.read_csv(file_path)

    # detect lipid classes
    all_cols = dataset.columns[start_idx:]

    numeric_cols = [
        c for c in all_cols
        if pd.api.types.is_numeric_dtype(dataset[c])
    ]

    lipids = sorted({
        get_lipid_class(c)
        for c in numeric_cols
        if get_lipid_class(c) is not None
    })

    # Clear existing checkboxes
    global lipid_checkboxes
    for widget in lipid_checkbox_frame.winfo_children():
        widget.destroy()
    lipid_checkboxes.clear()

    # Create checkbox for each lipid
    for lipid in lipids:
        var = tk.BooleanVar(value=True)  # Auto-select all by default
        lipid_checkboxes[lipid] = var

        cb = tk.Checkbutton(
            lipid_checkbox_frame,
            text=lipid,
            variable=var,
            font=("Arial", 13),
            bg="white",
            anchor="w",
            activebackground="white",
            relief=tk.FLAT,
            highlightthickness=0,
            bd=0,
            padx=10,
            pady=5
        )
        cb.pack(fill=tk.X, padx=2, pady=1)

    # Enable and update plot button styling to match Load CSV
    global plot_button
    plot_button.config(state=tk.NORMAL, bg="SystemButtonFace", fg="black")

def get_lipid_class(col_name):
    match = re.match(r"[A-Za-z]+", col_name)
    return match.group(0) if match else None

def format_columns(items, n_rows=10, col_width=22):
    """
    items: list of (label, value) tuples
    n_rows: max rows per column
    col_width: character width of each column
    """
    lines = []
    for i in range(0, len(items), n_rows):
        chunk = items[i:i + n_rows]
        col = [
            f"{k}: {v:.3f}".ljust(col_width)
            for k, v in chunk
        ]
        lines.append(col)

    # transpose columns → rows
    return "\n".join(
        "".join(row) for row in zip(*lines)
    )

def plot_averages():
    global dataset, bar_chart_frame

    if dataset is None:
        messagebox.showerror("Error", "Load a CSV first.")
        return

    df = dataset

    # numeric lipid columns
    all_numeric_cols = df.select_dtypes(include=np.number).columns
    numeric_cols = all_numeric_cols[start_idx:]

    # lipid class : list of columns
    lipid_groups = {}
    for col in numeric_cols:
        lipid = get_lipid_class(col)
        lipid_groups.setdefault(lipid, []).append(col)

    # Get selected lipids from checkboxes
    selected_lipids = [lipid for lipid, var in lipid_checkboxes.items() if var.get()]

    if not selected_lipids:
        messagebox.showerror("Error", "Select at least one lipid to display.")
        return

    # Sort lipids by total value (largest first, so they appear at bottom)
    lipid_totals = {}
    for lipid in selected_lipids:
        cols = lipid_groups.get(lipid, [])
        if cols:
            total_value = df[cols].mean().sum()
            lipid_totals[lipid] = total_value
        else:
            lipid_totals[lipid] = 0

    selected_lipids = sorted(selected_lipids, key=lambda x: lipid_totals[x], reverse=True)

    fig = Figure(figsize=(4, 6))
    ax = fig.add_subplot(111)

    bottom = 0
    bars = []

    cmap = plt.get_cmap("tab20")
    colors = cmap(np.arange(len(selected_lipids)))

    for lipid, color in zip(selected_lipids, colors):
        cols = lipid_groups.get(lipid, [])

        if not cols:
            continue
            
        species_means = df[cols].mean()   # CE(16:0), CE(16:1), ...
        total_value = species_means.sum()

        bar = ax.bar(
            [0],
            total_value,
            bottom=bottom,
            label=lipid,
            color=color,
            picker=50  # Very large pick tolerance for easy clicking
        )

        rect = bar[0]
        rect.lipid_class = lipid
        rect.avg_data = species_means 

        bars.append(rect)
        bottom += total_value


    ax.set_xticks([0])
    ax.set_xticklabels(["All Mice"])
    ax.set_ylabel("Average lipid abundance")
    ax.set_title("Lipid class composition (all mice)")

    ax.legend(
        title="Lipid class",
        bbox_to_anchor=(1.05, 1),
        loc="upper left"
    )

    # Adjust layout to prevent legend cutoff
    fig.tight_layout(rect=[0, 0, 0.85, 1])

    # clear old plot
    for widget in bar_chart_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=bar_chart_frame)

    # Use button_press_event for immediate response
    def on_canvas_click(event):
        if event.inaxes == ax and event.xdata is not None and event.ydata is not None:
            # Find which bar was clicked based on y-coordinate
            y_pos = event.ydata
            cumulative = 0
            for bar_rect in bars:
                bar_height = bar_rect.get_height()
                if cumulative <= y_pos <= cumulative + bar_height:
                    lipid_class = bar_rect.lipid_class
                    species_data = bar_rect.avg_data
                    show_heatmap(lipid_class, species_data)
                    break
                cumulative += bar_height

    canvas.mpl_connect("button_press_event", on_canvas_click)
    canvas.draw()
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.config(cursor="hand2")  # Show hand cursor to indicate clickable
    canvas_widget.pack(fill=tk.BOTH, expand=True)

def show_heatmap(lipid_class, species_data):
    """Display heatmap showing species breakdown for selected lipid class"""
    global pie_chart_container

    # Clear existing widgets
    for widget in pie_chart_container.winfo_children():
        widget.destroy()

    # Create close button frame at top
    button_frame = tk.Frame(pie_chart_container, bg="white")
    button_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

    close_btn = tk.Button(
        button_frame,
        text="✕ Close",
        command=close_heatmap,
        font=("Arial", 10),
        relief=tk.FLAT,
        bg="#f0f0f0"
    )
    close_btn.pack(side=tk.RIGHT)

    # Sort species by abundance (highest first)
    species_data = species_data.sort_values(ascending=False)

    labels = species_data.index.tolist()
    values = species_data.values

    # Reshape data for heatmap (n_species x 1 matrix)
    data_2d = values.reshape(-1, 1)

    # Calculate figure height - no cap, allow scrolling for many species
    n_species = len(labels)
    fig_height = max(6, n_species * 0.35)  # No upper limit

    # Create scrollable container for heatmap
    scroll_container = tk.Frame(pie_chart_container, bg="white")
    scroll_container.pack(fill=tk.BOTH, expand=True)

    # Create canvas for scrolling
    heatmap_canvas = tk.Canvas(scroll_container, bg="white", highlightthickness=0)
    heatmap_scrollbar = tk.Scrollbar(scroll_container, orient="vertical", command=heatmap_canvas.yview)

    # Create frame inside canvas to hold the matplotlib figure
    scrollable_heatmap_frame = tk.Frame(heatmap_canvas, bg="white")

    scrollable_heatmap_frame.bind(
        "<Configure>",
        lambda e: heatmap_canvas.configure(scrollregion=heatmap_canvas.bbox("all"))
    )

    heatmap_canvas.create_window((0, 0), window=scrollable_heatmap_frame, anchor="nw")
    heatmap_canvas.configure(yscrollcommand=heatmap_scrollbar.set)

    # Pack scrollbar and canvas
    heatmap_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    heatmap_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    # Create heatmap (50% width for 50/50 split)
    fig = Figure(figsize=(6, fig_height))
    ax = fig.add_subplot(111)

    # Create heatmap with YlOrRd colormap
    im = ax.imshow(data_2d, aspect='auto', cmap='YlOrRd', vmin=0)

    # Set ticks and labels
    ax.set_yticks(np.arange(n_species))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xticks([0])
    ax.set_xticklabels(['Abundance'], fontsize=10)

    # Add value annotations in each cell
    for i in range(n_species):
        value = data_2d[i, 0]
        # Use white text for dark cells, black for light cells
        text_color = 'white' if value > data_2d.max() * 0.5 else 'black'
        ax.text(0, i, f'{value:.3f}',
                ha="center", va="center",
                color=text_color, fontsize=9, fontweight='bold')

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, label='Abundance', pad=0.02)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(f"Species Breakdown: {lipid_class}", fontsize=12, fontweight='bold', pad=10)

    # Adjust layout to fit everything
    fig.tight_layout()

    # Create matplotlib canvas and embed in scrollable frame
    mpl_canvas = FigureCanvasTkAgg(fig, master=scrollable_heatmap_frame)
    mpl_canvas.draw()
    mpl_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Enable mouse wheel scrolling for heatmap
    def _on_heatmap_mousewheel(event):
        heatmap_canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def _bind_heatmap_mousewheel(event):
        heatmap_canvas.bind_all("<MouseWheel>", _on_heatmap_mousewheel)

    def _unbind_heatmap_mousewheel(event):
        heatmap_canvas.unbind_all("<MouseWheel>")

    heatmap_canvas.bind("<Enter>", _bind_heatmap_mousewheel)
    heatmap_canvas.bind("<Leave>", _unbind_heatmap_mousewheel)

    # Show the heatmap container (50% split with bar chart)
    pie_chart_container.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10)

def close_heatmap():
    """Hide the heatmap and return to single bar chart view"""
    global pie_chart_container

    # Hide the container
    pie_chart_container.pack_forget()

    # Clear widgets
    for widget in pie_chart_container.winfo_children():
        widget.destroy()

application = tk.Tk()
strand_var = tk.StringVar()

application.title("Data Visualization Tool")
application.geometry("1790x950")
application.resizable(False, False)
application.configure(bg="white")

# Main container with left controls and right plots
main_container = tk.Frame(application, bg="white")
main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Left control panel
control_panel = tk.Frame(main_container, bg="white", width=250)
control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
control_panel.pack_propagate(False)  # Maintain fixed width

# Load CSV button at top
tk.Button(
    control_panel,
    text="Load CSV",
    command=load_csv,
    width=20,
    height=2,
    font=("Arial", 10, "bold")
).pack(pady=(0, 10))

# Lipid selection label
tk.Label(
    control_panel,
    text="Select Lipids to Display:",
    font=("Arial", 10, "bold"),
    bg="white"
).pack(pady=(0, 5))

# Lipid selection with scrollable checkboxes
checkbox_outer_frame = tk.Frame(control_panel, bg="white")
checkbox_outer_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

# Create canvas with scrollbar
canvas = tk.Canvas(checkbox_outer_frame, bg="white", highlightthickness=0)
scrollbar = tk.Scrollbar(checkbox_outer_frame, orient="vertical", command=canvas.yview)
scrollable_frame = tk.Frame(canvas, bg="white")

scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
canvas.configure(yscrollcommand=scrollbar.set)

# Enable mouse wheel scrolling (bind to canvas only, not all widgets)
def _on_mousewheel(event):
    canvas.yview_scroll(int(-1*(event.delta/120)), "units")

def _bind_mousewheel(event):
    canvas.bind_all("<MouseWheel>", _on_mousewheel)

def _unbind_mousewheel(event):
    canvas.unbind_all("<MouseWheel>")

canvas.bind("<Enter>", _bind_mousewheel)
canvas.bind("<Leave>", _unbind_mousewheel)

canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

# Store reference to scrollable frame for adding checkboxes
lipid_checkbox_frame = scrollable_frame

# Plot button at bottom (disabled initially with white text)
plot_button = tk.Button(
    control_panel,
    text="Plot",
    command=plot_averages,
    width=20,
    height=2,
    font=("Arial", 10, "bold"),
    state=tk.DISABLED,
    fg="white",
    disabledforeground="white"
)
plot_button.pack()

# Right side: plots container
plots_container = tk.Frame(main_container, bg="white")
plots_container.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

# Left side: bar chart (always visible)
bar_chart_frame = tk.Frame(plots_container, bg="white")
bar_chart_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

# Right side: pie chart (initially hidden)
pie_chart_container = tk.Frame(plots_container, bg="white")
# Don't pack initially - will be shown when bar is clicked

application.mainloop()