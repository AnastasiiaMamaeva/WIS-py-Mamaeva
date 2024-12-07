file_name = input("input path to the file")

with open(file_name, "r") as f:
    seq = f.read()


def count_letters(seq):
    letters_dict = dict(zip(["A", "T", "G", "C"], [0]*4))
    for let in seq:
        for i in letters_dict:
            if i == let:         
                letters_dict[i] += 1
    return letters_dict

count_letters(seq)
