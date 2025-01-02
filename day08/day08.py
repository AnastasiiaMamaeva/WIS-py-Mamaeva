import tkinter as tk
from tkinter import filedialog

def count_letters(seq):
    letters_dict = dict(zip(["A", "T", "G", "C", "Unknown", "Total"], [0] * 6))
    unknown = len(seq)
    letters_dict["Total"] = len(seq)
    for let in seq:
        if let in ["A", "T", "G", "C"]:
            letters_dict[let] += 1
            unknown -= 1
    letters_dict["Unknown"] = unknown
    return letters_dict


def get_stat(dict_):
    return {k: round(dict_.get(k, 0) * 100 / max(dict_.get("Total", 1), 1), 2) for k in dict_}


def main():
    file_paths = filedialog.askopenfilenames(title="Select Text Files", filetypes=[("Text files", "*.txt")])
    if not file_paths:
        return  # Exit if no files are selected

    text_list.delete(0, tk.END)
    seq_all = {filename: open(filename, "r").read() for filename in file_paths}

    total_dict = dict(zip(["A", "T", "G", "C", "Unknown", "Total"], [0] * 6))

    for filename, seq in seq_all.items():
        text_list.insert(tk.END, f"File: {filename}")
        letters_dict = count_letters(seq)
        stat = get_stat(letters_dict)
        for key, value in stat.items():
            text_list.insert(tk.END, f"{key}: {value}%")
        text_list.insert(tk.END, f"Number of letters: {letters_dict['Total']}")
        text_list.insert(tk.END, "")


        total_dict = {k: total_dict.get(k, 0) + letters_dict.get(k, 0) for k in total_dict}


    total_stat = get_stat(total_dict)
    text_list.insert(tk.END, "Total Statistics:")
    for key, value in total_stat.items():
        text_list.insert(tk.END, f"{key}: {value}%")
    text_list.insert(tk.END, f"Total Letters: {total_dict['Total']}")


if __name__ == "__main__":
    root = tk.Tk()
    root.title("DNA Counter")

    lbl = tk.Label(root, text="Select files containing DNA sequences:")
    lbl.grid(row=0, column=0)

    open_button = tk.Button(root, text="Open Files", command=main)
    open_button.grid(row=1, column=0)

    text_list = tk.Listbox(root, width=60, height=25)
    text_list.grid(row=2, column=0)

    root.mainloop()