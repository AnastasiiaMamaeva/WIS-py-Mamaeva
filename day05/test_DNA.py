from seq_counter import count_letters, get_stat
def test_stat():
    res = get_stat(count_letters("ACCTGXXCXXGGTTGGGTNNNNNNNTTACTGGGCXTTGTXX"))
    assert res == {'A': 4.88, 'T': 24.39, 'G': 24.39, 'C': 12.2, 'Unknown': 34.15, 'Total': 100.0}