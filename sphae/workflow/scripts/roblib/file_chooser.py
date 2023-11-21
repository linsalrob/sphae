"""
choose a file. can use this if one is not provided
"""
import tkinter as tk
from tkinter import filedialog

def choose_a_file(dialog_title="Choose a file..."):

    root = tk.Tk()

    filetypes = (
        ('Text files', '*.TXT'),
        ('All files', '*.*'),
    )

    filename = tk.filedialog.askopenfilename(
        title=dialog_title,
        filetypes=filetypes,
    )
    root.destroy()

    return filename


def write_a_file(dialog_title="Choose where to save the file..."):
    filetypes = (
        ('TSV files', '*.TSV'),
        ('XLS files', '*.XLS'),
        ('All files', '*.*'),
    )

    filename = tk.filedialog.asksaveasfilename(
        title=f'Choose where to save the file...',
        filetypes=filetypes,defaultextension=".tsv"
    )

    return filename