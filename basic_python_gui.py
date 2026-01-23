import tkinter as tk
from tkinter import filedialog, messagebox # https://docs.python.org/3/library/tkinter.html

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backend_bases import key_press_handler # https://matplotlib.org/stable/gallery/user_interfaces/embedding_in_tk_sgskip.html
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

data = None
strain = "A.J"
start_idx = 0
window_size = 10  # how many proteins to show at once

def load_csv():
    global data, strain_var
    file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
    if file_path:
        data = pd.read_csv(file_path)
        data['Strain Only'] = data['Strains'].str.split('_').str[0]
        strains = data['Strain Only'].unique()  # get all unique strain names
        strain_var.set(strains[0])             # default selection
        strain_menu['menu'].delete(0, 'end')   # clear old menu
        for s in strains:
            strain_menu['menu'].add_command(label=s, command=tk._setit(strain_var, s))
        messagebox.showinfo("Success", "CSV loaded successfully")

def plot_averages():
    global data, start_idx, window_size, strain
    if data is None:
        messagebox.showerror("Error", "Load a CSV first.")
        return
    
    data['Strain Only'] = data['Strains'].str.split('_').str[0]
    
    strain = strain_var.get()
    if strain not in data['Strain Only'].unique():
        return
    
    try:
        window_size = int(window_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Window size must be an integer.")
        return

    subset = data[data['Strain Only'] == strain]
    if subset.empty:
        messagebox.showerror("Error", "No rows found for that strain.")
        return

    averages = subset.mean(numeric_only=True)
    n_proteins = len(averages)
    end_idx = min(start_idx + window_size, n_proteins)
    plot_data = averages.iloc[start_idx:end_idx]

    fig = plt.Figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    plot_data.plot(kind='bar', ax=ax)

    ax.set_title(f"Average Protein Levels for {strain}")
    ax.set_xlabel("Proteins")
    ax.set_ylabel("Average")

    # Progress label
    progress_var.set(f"Showing proteins {start_idx+1}-{end_idx} of {n_proteins}")

    for widget in plot_frame.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

def move_left():
    global start_idx
    start_idx = max(0, start_idx - window_size)
    plot_averages()

def move_right():
    global start_idx
    if data is None:
        return
    subset = data[data['Strain Only'] == strain]
    averages = subset.mean(numeric_only=True)
    start_idx = min(len(averages) - window_size, start_idx + window_size)
    plot_averages()

application = tk.Tk()

application.title("Data Visualization Tool")
application.geometry("900x600")
# application.resizable(False, False)
application.configure(bg="white")

top_frame = tk.Frame(application)
top_frame.pack(pady=10)

tk.Button(top_frame, text="Load CSV", command=load_csv).grid(row=0, column=0, padx=5)

strain_var = tk.StringVar()
strain_var.set("Select Strain")  # default text before CSV is loaded

tk.Label(top_frame, text="Strain:").grid(row=0, column=2, padx=5)
strain_menu = tk.OptionMenu(top_frame, strain_var, "")  # empty for now, filled after CSV load
strain_menu.config(fg="black")
strain_menu.grid(row=0, column=3, padx=5)

tk.Label(top_frame, text="Window size:").grid(row=0, column=4, padx=5)
window_entry = tk.Entry(top_frame, width=5)
window_entry.insert(0, str(window_size))
window_entry.grid(row=0, column=5, padx=5)

tk.Button(top_frame, text="Plot", command=plot_averages).grid(row=0, column=6, padx=10)

# Left / Right buttons
tk.Button(top_frame, text="←", width=3, command=move_left).grid(row=0, column=7, padx=5)
tk.Button(top_frame, text="→", width=3, command=move_right).grid(row=0, column=8, padx=5)

progress_var = tk.StringVar()
progress_label = tk.Label(application, textvariable=progress_var)
progress_label.pack()

plot_frame = tk.Frame(application, bg="white")
plot_frame.pack(fill=tk.BOTH, expand=True)

application.mainloop()