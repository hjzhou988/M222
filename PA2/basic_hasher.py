from collections import defaultdict, Counter
import cPickle as pickle
from os.path import join, exists, splitext
import time
import os
import sys
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from CS_CM122.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref


def hash_end(end, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param end: A single end of a read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """
    key_length = len(genome_ht.keys()[0])
    end_pieces = [end[i * key_length: (i + 1) * key_length]
                  for i in range(len(end) / key_length)]  # List of fragments?

    hashed_read_locations = [genome_ht[read_piece]
                             for read_piece in end_pieces] # List of List of Index of end_pieces. 
    start_positions = [[x - i * key_length for x in hashed_read_locations[i]] # List of list of start positions
                       for i in range(len(hashed_read_locations))]  
    start_counter = Counter()

    for position_list in start_positions:
        start_counter.update(position_list)

    if not start_counter:
        return -1, 0
    else:
        best_alignment_location, best_alignment_count = \
            start_counter.most_common(1)[0]

    if best_alignment_count < 2:
        return -1, best_alignment_count
    else:
        return best_alignment_location, best_alignment_count


def hash_read(read, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param read: A single read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """

    oriented_reads = [(read[0][::i], read[1][::j]) for i, j in ((1, -1), (-1, 1))]
    ## Either one end is forward and the other end is reversed, or vice versa.

    best_score = -1
    best_alignment_locations = (-1, -1)
    best_oriented_read = ('', '')
    for oriented_read in oriented_reads:
        hash_results = [hash_end(_, genome_ht) for _ in oriented_read]
        hash_locations = [_[0] for _ in hash_results]
        hash_score = sum([_[1] for _ in hash_results]) # score is the sum of the best alignment count of the two reads]
        if hash_score > best_score:
            best_alignment_locations = hash_locations
            best_oriented_read = oriented_read
            best_score = hash_score  #insertion
    return best_oriented_read, best_alignment_locations # read pair, two locations


def make_genome_hash(reference, key_length):
    """

    :param reference: The reference as a string stored
    :param key_length: The length of keys to use.
    :return:
    """
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash


def build_hash_and_pickle(ref_fn, key_length, force_rebuild=False):
    reference_hash_pkl_fn = '{}_hash_keylength_{}.pkl'.format(splitext(ref_fn)[0], key_length)
    if exists(reference_hash_pkl_fn) and not force_rebuild:
        ref_genome_hash = pickle.load(open(reference_hash_pkl_fn, 'rb'))
        if len(ref_genome_hash.keys()[0]) == key_length:
            return ref_genome_hash
        else:
            pass
    else:
        pass
    reference = read_reference(ref_fn)
    ref_genome_hash = make_genome_hash(reference, key_length)
    pickle.dump(ref_genome_hash, open(reference_hash_pkl_fn, 'wb'))
    return ref_genome_hash


def hashing_algorithm(paired_end_reads, genome_ht):
    """

    :param paired_end_reads:
    :param genome_ht:
    :return:
    """
    alignments = []  # index
    genome_aligned_reads = []  # DNA fragments
    count = 0
    start = time.clock()

    for read in paired_end_reads:
        alignment, genome_aligned_read = hash_read(read, genome_ht) #best_oriented_read, best_alignment_locations
        alignments.append(alignment)
        genome_aligned_reads.append(genome_aligned_read)
        count += 1
        if count % 100 == 0:
            time_passed = (time.clock()-start)/60
            print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print 'Approximately {:.3} minutes remaining'.format(remaining_time)
    return alignments, genome_aligned_reads

def determine_STR_inRef(ref):
    STRlist=[]
    i=0
    potentialSTR=''
    for k in range(2,7):
        i=0
        while i<len(ref):
            potentialSTR=ref[i:i+k]
            copynumber=0
            if ref[i+k:i+6*k]==potentialSTR*5:
                if k==4 and potentialSTR[0:2]==potentialSTR[2:4]:
                    i+=1
                    continue
                if k==6 and potentialSTR[0:2]==potentialSTR[2:4]==potentialSTR[4:6]:
                    i+=1
                    continue
                if k==6 and potentialSTR[0:3]==potentialSTR[3:6]:
                    i+=1
                    continue
                copynumber=6
                j=i+6*k
                while True:
                    if ref[j:j+k]==potentialSTR:
                        copynumber+=1
                        j+=k
                    else:
                        STRlist.append(['STR',potentialSTR*copynumber,i,copynumber])
                        break
                i=i+k*copynumber+1  
            else:
                i+=1
    STRlist = sorted(STRlist, key=lambda STR: STR[-2])
    return STRlist

if __name__ == "__main__":
    genome_name = 'hw2undergrad_E_2'#'practice_W_3'##'hw1_W_2'#
    input_folder = '../data/{}'.format(genome_name)
    chr_name = '{}_chr_1'.format(genome_name)
    reads_fn_end = 'reads_{}.txt'.format(chr_name)
    reads_fn = join(input_folder, reads_fn_end)
    ref_fn_end = 'ref_{}.txt'.format(chr_name)
    ref_fn = join(input_folder, ref_fn_end)
    key_length = 7
    start = time.clock()
    reads = read_reads(reads_fn)
    # If you want to speed it up, cut down the number of reads by
    # changing the line to reads = read_reads(reads_fn)[:<x>] where <x>
    # is the number of reads you want to work with.
    genome_hash_table = build_hash_and_pickle(ref_fn, key_length)
    ref = read_reference(ref_fn)
    genome_aligned_reads, alignments = hashing_algorithm(reads, genome_hash_table)
    pickle.dump((genome_aligned_reads, alignments),open(join(input_folder,genome_name+'_aligned_reads_and_coordinates.pkl'),'wb'))
    # print genome_aligned_reads
    # print alignments
    output_str = pretty_print_aligned_reads_with_ref(genome_aligned_reads, alignments, ref)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(chr_name))
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)

    pickle.dump(determine_STR_inRef(ref),open(join(input_folder,genome_name+'_STR.pkl'),'wb'))
    print output_str[:5000]
