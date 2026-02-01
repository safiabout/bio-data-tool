import tkinter as tk
from tkinter import filedialog, messagebox # https://docs.python.org/3/library/tkinter.html

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backend_bases import key_press_handler # https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_tk_sgskip.html
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

datasets = {}
start_idx = 0
window_size = 10

def load_csv():
    global datasets, dataset_var, strand_var
    file_paths = filedialog.askopenfilenames(title="Select CSV files", filetypes=[("CSV files", "*.csv")], multiple=True)

    for path in file_paths:
        name = path.split("/")[-1].replace(".csv","")  # filename as friendly name
        datasets[name] = pd.read_csv(path)
    
    # Update the dataset dropdown
    dataset_menu['menu'].delete(0, 'end')  # clear old items
    for name in datasets.keys():
        dataset_menu['menu'].add_command(label=name, command=tk._setit(dataset_var, name))

    if datasets:
        first_dataset = list(datasets.keys())[0]
        dataset_var.set(first_dataset)

        update_strand_dropdown()
        plot_averages()
        messagebox.showinfo("Success", f"{len(datasets)} datasets loaded!")

def update_strand_dropdown():
    global strain_menu, strain_var
    dataset_name = dataset_var.get()

    df = datasets[dataset_name]
    df['Strain Only'] = df['Strains'].str.split('_').str[0]

    strains = df['Strain Only'].unique()
    strain_var.set(strains[0])
    strain_menu['menu'].delete(0, 'end')

    for s in strains:
        strain_menu['menu'].add_command(label=s, command=tk._setit(strain_var, s))

def plot_averages():
    global dataset_var, start_idx, window_size, strain
    if dataset_var is None:
        messagebox.showerror("Error", "Load a CSV first.")
        return
    
    df = datasets[dataset_var.get()]
    
    df['Strain Only'] = df['Strains'].str.split('_').str[0]
    
    strain = strain_var.get()
    if strain not in df['Strain Only'].unique():
        strain_var.set(df['Strain Only'].unique()[0])
        strain = strain_var.get()
    
    try:
        window_size = int(window_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Window size must be an integer.")
        return

    subset = df[df['Strain Only'] == strain_var.get()]
    if subset.empty:
        messagebox.showerror("Error", "No rows found for that strain.")
        return

    averages = subset.mean(numeric_only=True)
    n_proteins = len(averages)

    if start_idx >= n_proteins:
        start_idx = max(0, n_proteins - window_size)

    end_idx = min(start_idx + window_size, n_proteins)
    plot_data = averages.iloc[start_idx:end_idx]

    fig = plt.Figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    
    bars = ax.bar(plot_data.index, plot_data.values)

    for bar in bars:
        bar.set_picker(True)  # enable clicking

    ax.set_title(f"Average {dataset_var.get().split('_')[0]} for {strain}")
    ax.set_xlabel(dataset_var.get().split('_')[0])
    ax.set_ylabel("Average")

    # Progress label
    progress_var.set(f"Showing {dataset_var.get().split('_')[0]} range {start_idx+1}-{end_idx} of {n_proteins}")

    for widget in plot_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.mpl_connect("pick_event", on_bar_click)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

def move_left():
    global start_idx
    start_idx = max(0, start_idx - window_size)
    plot_averages()

def move_right():
    global start_idx
    if dataset_var is None:
        return
    
    df = datasets[dataset_var.get()]
    subset = df[df['Strain Only'] == strain_var.get()]
    averages = subset.mean(numeric_only=True)
    start_idx = min(len(averages) - window_size, start_idx + window_size)
    plot_averages()

def on_bar_click(event):
    bar = event.artist
    label = bar.get_x() + bar.get_width() / 2
    height = bar.get_height()

    protein = bar.axes.get_xticklabels()[int(bar.get_x() + 0.5)].get_text()

    messagebox.showinfo(
        "Bar Clicked",
        f"You clicked:\n{protein}\nValue: {height:.2f}"
    )

    # this is where we can have:
    # - heatmap
    # - species breakdown
    # - new window

application = tk.Tk()
strand_var = tk.StringVar()

application.title("Data Visualization Tool")
application.geometry("900x950")
application.resizable(False, False)
application.configure(bg="white")

top_frame = tk.Frame(application)
top_frame.pack(pady=10)

tk.Button(top_frame, text="Load CSV", command=load_csv).grid(row=0, column=0, padx=5)

dataset_var = tk.StringVar()
dataset_var.set("Select Dataset")  # default text
dataset_var.trace('w', lambda *args: [update_strand_dropdown(), plot_averages()])

dataset_menu = tk.OptionMenu(top_frame, dataset_var, "")
dataset_menu.config(fg="black")
dataset_menu.grid(row=0, column=2, padx=5)

strain_var = tk.StringVar()
strain_var.set("Select Strain")  # default text before CSV is loaded

tk.Label(top_frame, text="Strain:").grid(row=0, column=4, padx=5)
strain_menu = tk.OptionMenu(top_frame, strain_var, "")  # empty for now, filled after CSV load
strain_menu.config(fg="black")
strain_menu.grid(row=0, column=5, padx=5)

tk.Label(top_frame, text="Window size:").grid(row=0, column=6, padx=5)
window_entry = tk.Entry(top_frame, width=5)
window_entry.insert(0, str(window_size))
window_entry.grid(row=0, column=7, padx=5)

tk.Button(top_frame, text="Plot", command=plot_averages).grid(row=0, column=8, padx=10)

# Left / Right buttons
tk.Button(top_frame, text="←", width=3, command=move_left).grid(row=0, column=9, padx=5)
tk.Button(top_frame, text="→", width=3, command=move_right).grid(row=0, column=10, padx=5)

progress_var = tk.StringVar()
progress_label = tk.Label(application, textvariable=progress_var)
progress_label.pack()

plot_frame = tk.Frame(application, bg="white")
plot_frame.pack(fill=tk.BOTH, expand=True)

application.mainloop()