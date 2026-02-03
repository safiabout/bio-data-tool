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

    available_listbox.delete(0, tk.END)
    selected_listbox.delete(0, tk.END)

    for lipid in lipids:
        available_listbox.insert(tk.END, lipid)

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
    global dataset

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

    selected_lipids = selected_listbox.get(0, tk.END)

    if not selected_lipids:
        messagebox.showerror("Error", "Select at least one lipid to display.")
        return

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
            picker=True
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

    # clear old plot
    for widget in plot_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.mpl_connect("pick_event", on_bar_click)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    cursor = mplcursors.cursor(bars, hover=True)

    @cursor.connect("add")
    def on_add(sel):
        rect = sel.artist
        lipid = rect.lipid_class

        items = list(rect.avg_data.items())

        body = format_columns(
            items,
            n_rows=16,       # max lines per column
            col_width=24    # spacing between columns
        )

        sel.annotation.set_text(f"{lipid}\n{body}")
        sel.annotation.get_bbox_patch().set_alpha(0.9)
        sel.annotation.get_bbox_patch().set_boxstyle("round,pad=0.3")

def move_items(src, dst):
    selections = list(src.curselection())
    for i in reversed(selections):
        item = src.get(i)
        dst.insert(tk.END, item)
        src.delete(i)

def on_bar_click(event):
    lipid_class = event.artist.get_label()
    value = event.artist.get_height()

    messagebox.showinfo(
        "Lipid Class Selected",
        f"{lipid_class}\nTotal: {value:.2f}"
    )

    # NEXT STEP:
    # open new window with heatmap of subtypes and  within this class?

application = tk.Tk()
strand_var = tk.StringVar()

application.title("Data Visualization Tool")
application.geometry("1790x950")
application.resizable(False, False)
application.configure(bg="white")

list_frame = tk.Frame(application)
list_frame.pack(pady=10)

tk.Label(list_frame, text="Available Lipids").grid(row=0, column=0)
tk.Label(list_frame, text="Displayed Lipids").grid(row=0, column=2)

available_listbox = tk.Listbox(
    list_frame,
    selectmode=tk.MULTIPLE,
    height=10,
    exportselection=False
)
available_listbox.grid(row=1, column=0, padx=10)

selected_listbox = tk.Listbox(
    list_frame,
    selectmode=tk.MULTIPLE,
    height=10,
    exportselection=False
)
selected_listbox.grid(row=1, column=2, padx=10)

top_frame = tk.Frame(application)
top_frame.pack(pady=10)

tk.Button(top_frame, text="Load CSV", command=load_csv).grid(row=0, column=0, padx=5)

tk.Button(top_frame, text="Plot", command=plot_averages).grid(row=0, column=8, padx=10)

btn_frame = tk.Frame(list_frame)
btn_frame.grid(row=1, column=1)

tk.Button(btn_frame, text="→", command=lambda: move_items(
    available_listbox, selected_listbox
)).pack(pady=5)

tk.Button(btn_frame, text="←", command=lambda: move_items(
    selected_listbox, available_listbox
)).pack(pady=5)

plot_frame = tk.Frame(application, bg="white")
plot_frame.pack(fill=tk.BOTH, expand=True)

application.mainloop()