import sys

def find_longest_duplicate(sequence):
    n = len(sequence)
    longest = ""
    
    for length in range(1, n):
        for i in range(n - length):
            substring = sequence[i:i + length]
            if sequence.count(substring) > 1 and len(substring) > len(longest):
                longest = substring
    
    return longest

def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0
    return gc_content

def find_six_letter_dna_palindromes(sequence):
    palindromes = []
    complement = str.maketrans('ATCG', 'TAGC')
    n = len(sequence)
    for i in range(n - 5):
        substring = sequence[i:i + 6]
        if substring == substring[::-1].translate(complement):
            palindromes.append(substring)
    return palindromes

def parse_file(filepath):
    sequence = ""
    with open(filepath, 'r') as file:
        for line in file:
            if not line.startswith('>') and not line.startswith('LOCUS'):
                sequence += line.strip()
    return sequence

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python analyze.py FILE [--gc-content] [--palindromes]")
        sys.exit(1)

    filepath = sys.argv[1]
    options = sys.argv[2:]
    sequence = parse_file(filepath)


    
    longest_duplicate = find_longest_duplicate(sequence)
    print(f"Longest repeated subsequence: {longest_duplicate}")

    if "--gc-content" in options:
        gc_content = calculate_gc_content(sequence)
        print(f"GC content: {gc_content:.2f}%")

    if "--palindromes" in options:
        palindromes = find_six_letter_dna_palindromes(sequence)
        print(f"6-letter DNA palindromes: {', '.join(palindromes) if palindromes else 'None found'}")


