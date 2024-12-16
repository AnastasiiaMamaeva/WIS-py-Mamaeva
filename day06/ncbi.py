def get_input():
    if len(sys.argv) < 2:
        print("please tipe: python ncbi.py  `--database protein --term TERM --number NUMBER` \n or at least `python ncbi.py --term TERM`")
        sys.exit(1)
    parsed = {}
    for i in ["--database", "--term","--number"]:
        if i in sys.argv:
            parsed[i] = sys.argv[sys.argv.index(i) +1]
        if ("--database" not in parsed.keys()):
            parsed["--database"] = "protein"
        if ("--number" not in parsed.keys()):
            parsed["--number"] = 1
    return parsed

def get_id_and_rec(db, term, max):
    handle = Entrez.esearch(db=db, term=term)  # Set retmax=0 to avoid fetching IDs
    record = Entrez.read(handle)
    handle.close()
    total = int(record["Count"])
    ids  = record["IdList"][:max]
    return ids, total


def get_data(db, ids):

    handle = Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text")
    data = handle.read()
    handle.close()
    return data


def save_data(filename, data):
    with open(filename, 'w') as fh:
        fh.write(data)

if __name__ == "__main__":
    import sys
    from Bio import Entrez
    import datetime
    
    Entrez.email = "mamaeva.anastas@gmail.com"
    
    parsed = get_input()
    term = parsed["--term"]
    max = int(parsed["--number"])
    db = parsed["--database"]
    
    ids, total = get_id_and_rec(db, term, max)
    f_names = []
    for it in ids:
        f_name = it+ ".fasta"
        f_names.append(f_name)
        data = get_data(db, [it])
        save_data(f_name,data)
    print(datetime.datetime.now(),term,max,total, sep = ",")
    print(f"file names:{f_names}")