def count_letters(seq):
    letters_dict = dict(zip(["A", "T", "G", "C","Unknown", "Total"], [0]*6))
    Unknown = len(seq)
    letters_dict["Total"] = len(seq)
    for let in seq:
        for i in letters_dict:
            if i == let:         
                letters_dict[i] += 1
                Unknown -=1 
            
    letters_dict["Unknown"] = Unknown
    return letters_dict


def get_stat(dict_):
    proc = {k:round(dict_.get(k,0)*100/dict_.get("Total", 0),2) for k in list(dict_.keys())}
    return proc

if __name__ == "__main__":
    import sys
    def get_input():
        if len(sys.argv) < 2:
            print("please tipe: python seq_counter.py <filename>")
            sys.exit(1)
        seq_all = []
        for filename in sys.argv[1:]:
            with open(filename, "r") as fh:
                seq_init = fh.read()
                seq_all.append(seq_init) 
            return seq_all
    
    seq_all = get_input()
    total_dict = dict(zip(["A", "T", "G", "C","Unknown", "Total"], [0]*6))
    i = 1
    for seq in seq_all:
        print(sys.argv[i])
        letters_dict = count_letters(seq)
        stat = get_stat(letters_dict)
        stat_keys = list(stat.keys())
        stat_values = list(stat.values())
        for t in range(5):
            print(f"{stat_keys[t]}:\t{stat_values[t]}")
        print(f"Number of letters:\t {letters_dict["Total"]}\n")
        total_dict = {k: total_dict.get(k, 0) + letters_dict.get(k, 0) 
        for k in ["A", "T", "G", "C", "Unknown", "Total"]}
        i +=1
    total_stat = get_stat(total_dict)

    if __name__ == "__main__":
        print("Total:")
        for t in range(5):
            print(f"{list(total_stat.keys())[t]}:\t{
                list(total_stat.values())[t]}")
        print(f"Number of letters:\t {total_dict["Total"]}")
