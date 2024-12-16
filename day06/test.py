from Bio import Entrez

# Set your email
Entrez.email = "mamaeva.anastas@gmail.com"

# Define the search query and database
query = "Homo sapiens"
db = "protein"
amount_to_fetch = 10  # Specify how many records to fetch

# Step 1: Search the database and fetch IDs
handle = Entrez.esearch(db=db, term=query, retmax=amount_to_fetch)
record = Entrez.read(handle)
handle.close()

# Extract the list of IDs
ids = record["IdList"]
print(f"Fetched IDs: {ids}")

# Step 2: Fetch the records for these IDs
if ids:
    handle = Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text")  # Change rettype as needed
    data = handle.read()
    handle.close()

    # Print or process the fetched data
    print(data, sep = "\n and \n")
else:
    print("No records found.")
data.split(">")
"".join([i, "fasta"])