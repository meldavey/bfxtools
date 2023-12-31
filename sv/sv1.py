import os
import numpy as np

verbose = 0

def load_fasta(filename):
    num_bytes = os.stat(filename).st_size
    print('file is %d bytes' % num_bytes)
    data = bytearray()
    cur = 0
    f = open(filename, 'rb')
    line = f.readline() # header, ignore for now
    done = False
    while not done:
        line = f.readline()
        l = len(line)
        if l > 0:
            data.extend(line[:-1])
            cur += (l-1)
        else:
            done = True
    f.close()
    print('loaded %d bases' % cur)

    return data


def get_pattern_locations(data):
    loc = []
    pattern = b'cgat'
    pattern_len = len(pattern)
    data_len = len(data)
    search_len = data_len - pattern_len + 1
    cur = 0
    while cur < search_len:
        if data[cur:cur+pattern_len] == pattern:
            loc.append(cur)
        cur += 1

    return loc

def get_dist_from_locations(pattern_locations):
    dist = []
    for i in range(len(pattern_locations)-1):
        dist.append(pattern_locations[i+1]-pattern_locations[i])
    return dist

# find_starts
# finds all candidate start positions in the distance array that match the pattern within some threshold
# initial implementation is simple error at each distance vs pattern and threshold
#
# notes - probably need to anchor on first distance, then compare each distance relative to anchor
# also can do simple edit distance to allow for an insertion/deletion of distance entry and still match

def find_starts(pattern, dist):
    candidate_positions = []

    pattern_len = len(pattern)
    dist_len = len(dist)
    search_len = dist_len - pattern_len + 1
    thresh = pattern_len * 5

    if verbose > 1:
        print('looking for pattern: %s' % str(pattern))

    for i in range(search_len):
        err = 0
        for j in range(pattern_len):
            err += abs(dist[i+j] - pattern[j])
        if err < thresh:
            candidate_positions.append(i)
    return candidate_positions


# very simple test first, generate pattern locations, then use the start & end as anchors and compare the read's mapped position to the genome
def eval_read(read, ref_dist, ref_patloc):
    read_patloc = get_pattern_locations(read)
    read_dist = get_dist_from_locations(read_patloc)
    read_starts = find_starts(read_dist[:6], ref_dist)
    read_ends = find_starts(read_dist[-6:], ref_dist)

    if verbose > 0:
        print('insertion starts: %s ends: %s' % (str(read_starts), str(read_ends)))
    read_start = read_starts[0]
    read_end = read_ends[0]
    delta_ref = ref_patloc[read_end] - ref_patloc[read_start]
    delta_read = read_patloc[-6] - read_patloc[0]
    delta_len = delta_read - delta_ref
    print('%s is roughly %dbp in length' % ('insertion' if delta_len > 0 else 'deletion', abs(delta_len)))


# ---------

# load a reference genome
genome = load_fasta('ecoli.fasta')

# generate location array where we find each pattern in the genome
pattern_locations = get_pattern_locations(genome)
print('found %d occurences of pattern' % len(pattern_locations))

# find distances between pattern starts
dist = get_dist_from_locations(pattern_locations)

# test 1 - use a known pattern, then try and find it in our distance array, as well as other similar ones
pattern_start = 50
pattern_len = 6 # this represents len+1 pattern locations since its the spacing between them
pattern = dist[pattern_start:pattern_start+pattern_len]
print('pattern: %s' % str(pattern))

starts = find_starts(pattern, dist)
print('found %d candidate starts' % len(starts))
print('first few starts: %s' % str(starts[:5]))

# detect insertions, deletions, inversions
# not clear yet on translocations, but likely this is just whether the read maps to the known location in the genome or not

# generate an insertion by starting with a portion of the genome and inserting random bases
np.random.seed(1134)
base_lookup = (b'a',b'c',b'g',b't')
read_insertion = genome[1000000:2000000]
insertion_len = 20000
for i in range(insertion_len):
    base = base_lookup[np.random.randint(0,4)]
    read_insertion.extend(base)
read_insertion.extend(genome[2000000:2500000])
print('read insertion is length: %d' % len(read_insertion))

# generate a 100kb deletion
read_deletion = genome[400000:400000+200000]
read_deletion.extend(genome[700000:1000000])

eval_read(read_insertion, dist, pattern_locations)
eval_read(read_deletion, dist, pattern_locations)

