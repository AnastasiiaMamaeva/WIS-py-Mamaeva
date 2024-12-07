from seq_counter import count_letters, get_stat
def test_count_letters():
    res = count_letters("CCCGTGGGTGCCCTTTGGGTC")
    assert res == {'A': 0, 'T': 6, 'G': 8, 'C': 7, 'Unknown': 0, 'Total': 21}

def test_stat():
    res = get_stat(count_letters("ACCTGXXCXXGGTTGGGTNNNNNNNTTACTGGGCXTTGTXX"))
    assert res == {'A': 4.88, 'T': 24.39, 'G': 24.39, 'C': 12.2, 'Unknown': 34.15, 'Total': 100.0}

def test_stat2():
    res = get_stat(count_letters("CCCGTGGGTGCCCTTTGGGTC"))
    assert res == {'A': 0.0, 'T': 28.57, 'G': 38.1, 'C': 33.33, 'Unknown': 0.0, 'Total': 100.0}

test_count_letters()
test_stat()
test_stat2()
print("all tests are passed")