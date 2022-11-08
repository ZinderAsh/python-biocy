from Graph import hash_kmer
import pytest

@pytest.mark.parametrize("kmer,expected", [
            ('ACTAC', 97),
            ('AAAAAAAAAA', 0),
            ('CTACC', 389),
            ('ACTGC', 109),
            ('CTGCC', 437),
            ('A', 0),
            ('C', 1),
            ('T', 2),
            ('G', 3),
            ])
def test_hash_kmer(kmer, expected):
    hashed = hash_kmer(kmer.encode('ASCII'), len(kmer))
    assert hashed == expected
