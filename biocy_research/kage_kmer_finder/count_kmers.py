import sys
import numpy as np
from graph_kmer_index import kmer_hash_to_sequence

infile1 = sys.argv[1]
infile2 = sys.argv[2]

indata1 = np.fromfile(infile1, dtype=np.uint64)
indata2 = np.fromfile(infile2, dtype=np.uint64)

counts1 = {}

print(f"Counting Kmers 1 (total: {len(indata1)})")
for i in indata1:
    if i not in counts1:
        counts1[i] = 1
    else:
        counts1[i] += 1
print(len(counts1), "unique kmers")

counts2 = {}

print(f"Counting Kmers 2 (total: {len(indata2)})")
for i in indata2:
    if i not in counts2:
        counts2[i] = 1
    else:
        counts2[i] += 1
print(len(counts2), "unique kmers")

print("Merging Keys")
keys = set(counts1.keys())
keys.update(counts2.keys())

print(len(keys), "total unique kmers")

def kmer_to_string(kmer):
    kmer_map = "ACGT"
    test_str = ""
    for i in range(31):
        test_str += kmer_map[np.bitwise_and(np.right_shift(kmer, np.uint64(i * 2)), np.uint64(3))]
    print("Test Str:", test_str)

print("Comparing Counts")
same = 0
diff = 0
not1 = 0
not2 = 0
for key in keys:
    if key not in counts1:
        not1 += 1
    elif key not in counts2:
        not2 += 1
    elif counts1[key] == counts2[key]:
        same += 1
    else:
        diff += 1

print("Same:", same)
print("Diff:", diff)
print("Missing From 1:", not1)
print("Missing From 2:", not2)

