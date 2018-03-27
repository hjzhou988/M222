import sys
import numpy as np
from collections import defaultdict
import time
from os.path import join
from basic_hasher import build_hash_and_pickle, hashing_algorithm
import os
import zipfile
import cPickle as pickle
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from CS_CM122.helpers import *

READ_LENGTH = 50


def generate_pileup(aligned_fn):
    """
    :param aligned_fn: The filename of the saved output of the basic aligner
    :return: SNPs (the called SNPs for uploading to the herokuapp server)
             output_lines (the reference, reads, consensus string, and diff string to be printed)
    """
    line_count = 0
    lines_to_process = []
    changes = []
    start = time.clock()
    with open(aligned_fn, 'r') as input_file:
        for line in input_file:
            line_count += 1
            line = line.strip() # this step will strip off the white spaces after the read (if there are any). So in the later step, need to add white spaces for " padding"?
            if line_count <= 4 or line == '':  # The first 4 lines need to be skipped
                continue
            if len(line) > 0 and all(x == '-' for x in line):  # The different pieces of the genome are set off
                # with lines of all dashes '--------'
                new_changes = process_lines(lines_to_process)
                lines_to_process = []
                changes += new_changes
                # print time.clock() - start, 'seconds'
            else:
                lines_to_process.append(line)
    snps = [v for v in changes if v[0] == 'SNP']
    insertions = [v for v in changes if v[0] == 'INS']
    deletions = [v for v in changes if v[0] == 'DEL']
    return snps, insertions, deletions


def process_lines(genome_lines):
    """

    :param genome_lines: Lines in between dashes from the saved output of the basic_aligner
    :return: snps (the snps from this set of lines)
             output_lines (the lines to print, given this set of lines)
    """
    line_count = 0
    consensus_lines = []
    for line in genome_lines:
        line_count += 1
        if line_count == 1:  # The first line contains the position in the reference where the reads start.
            raw_index = line.split(':')[1]
            line_index = int(raw_index)
        else:
            consensus_lines.append(line[6:])
    ref = consensus_lines[0]
    aligned_reads = consensus_lines[1:]
    donor = generate_donor(ref, aligned_reads)
    print donor
    changes = alignment(ref,donor,line_index) #identify_changes(ref, donor, line_index)
    return changes


def align_to_donor(donor, read):
    """
    :param donor: Donor genome (a character string of A's, T's, C's, and G's, and spaces to represent unknown bases).
    :param read: A single read padded with spaces
    :return: The best scoring
    """

    mismatches = [1 if donor[i] != ' ' and read[i] != ' ' and
                       read[i] != donor[i] else 0 for i in range(len(donor))]
    n_mismatches = sum(mismatches)
    overlaps = [1 if donor[i] != ' ' and read[i] != ' ' else 0 for i in range(len(donor))]
    n_overlaps = sum(overlaps)
    score = n_overlaps - n_mismatches #number of matches
    if n_mismatches <= 2:
        return read, score
    else: #n_mismatches >2
        best_read = read
        best_score = score

    for shift_amount in range(-3, 0) + range(1, 4):  # This can be improved
        if shift_amount > 0:
            shifted_read = ' ' * shift_amount + read
        elif shift_amount < 0:
            shifted_read = read[-shift_amount:] + ' ' * (-shift_amount)
        mismatches = [1 if donor[i] != ' ' and shifted_read[i] != ' ' and
                           shifted_read[i] != donor[i] else 0 for i in range(len(donor))]
        n_mismatches = sum(mismatches)
        overlaps = [1 if donor[i] != ' ' and shifted_read[i] != ' ' else 0 for i in range(len(donor))]
        n_overlaps = sum(overlaps)
        score = n_overlaps - n_mismatches #- abs(shift_amount) # why 3 times?
        if score > best_score:
            best_read = shifted_read
            best_score = score
    return best_read, best_score


def generate_donor(ref, aligned_reads):
    """
    Aligns the reads against *each other* to generate a hypothesized donor genome.
    There are lots of opportunities to improve this function.
    :param aligned_reads: reads aligned to the genome (with pre-pended spaces to offset correctly)
    :return: hypothesized donor genome
    """
    cleaned_aligned_reads = [_.replace('.', ' ') for _ in aligned_reads]

    ## Start by appending spaces to the reads so they line up with the reference correctly.
    padded_reads = [aligned_read + ' ' * (len(ref) - len(aligned_read)) for aligned_read in cleaned_aligned_reads]
    consensus_string = consensus(ref, aligned_reads)

    ## Seed the donor by choosing the read that best aligns to the reference. (The best alin)
    read_scores = [sum([1 if padded_read[i] == ref[i] and padded_read[i] != ' '
                        else 0 for i in range(len(padded_read))])
                   for padded_read in padded_reads]  # for every padded read, it has a read_score. The score is the number of positions that align to the reference. It could be imperfect. 
    if not read_scores:  # What does this mean? If there is no padded read?
        return consensus_string
    longest_read = padded_reads[read_scores.index(max(read_scores))]
    donor_genome = longest_read

    # While there are reads that haven't been aligned, try to align them to the donor.
    while padded_reads:
        un_donored_reads = []
        for padded_read in padded_reads:
            re_aligned_read, score = align_to_donor(donor_genome, padded_read)
            if score < 10:  # If the alignment isn't good, throw the read back in the set of reads to be aligned.
                un_donored_reads.append(padded_read)
            else:
                donor_genome = ''.join([re_aligned_read[i] if donor_genome[i] == ' ' else donor_genome[i]
                                        for i in range(len(donor_genome))])

        if len(un_donored_reads) == len(padded_reads):
            # If we can't find good alignments for the remaining reads, quit
            break
        else:
            # Otherwise, restart the alignment with a smaller set of unaligned reads
            padded_reads = un_donored_reads

    if len(un_donored_reads) >=20:  # if there are plenty of un_donored_reads that were not aligned to donor genome, it means there is a big insertion segment that prevent this alignment.
        read_scores=[sum([1 if un_donored_read[i] == ref[i] and un_donored_read[i] != ' '
                        else 0 for i in range(len(un_donored_read))])
                   for un_donored_read in un_donored_reads]
        longest_read = un_donored_reads[read_scores.index(max(read_scores))]
        donor_genome_2 = longest_read
        padded_reads=un_donored_reads
        while padded_reads:
            un_donored_reads = []
            for padded_read in padded_reads:
                re_aligned_read, score = align_to_donor(donor_genome_2, padded_read)
                if score < 10:  # If the alignment isn't good, throw the read back in the set of reads to be aligned.
                    un_donored_reads.append(padded_read)
                else:
                    donor_genome_2 = ''.join([re_aligned_read[i] if donor_genome_2[i] == ' ' else donor_genome_2[i]
                                            for i in range(len(donor_genome_2))])

            if len(un_donored_reads) == len(padded_reads):
                # If we can't find good alignments for the remaining reads, quit
                break
            else:
                # Otherwise, restart the alignment with a smaller set of unaligned reads
                padded_reads = un_donored_reads

    #print donor_genome_2
    ## Fill in any gaps with the consensus sequence and return the donor genome.
    #donor_genome = ''.join([donor_genome[i] if donor_genome[i] != ' ' else consensus_string[i] for i
    #                       in range(len(donor_genome))]) # It seemed that donor_genome has the same length as consensus_string? 
    #print donor_genome
    #donor genome might incorporate sequencing errors, but it could be more accurate in terms of insertions and deletions. 
    #consensus sequence has the same length as reference sequnce. If length of Donar genome equals to length of consensus sequence, then use cnosensus sequence. 
    #donor_genome=donor_genome.strip()
        new_donor_genome=''.join([donor_genome_2[i] if donor_genome[i]==' ' else donor_genome[i] for i in range(len(donor_genome))])
        print new_donor_genome
        #if new_donor_genome==donor_genome:  

        return new_donor_genome

    if len(donor_genome.strip())==len(consensus_string):
        return consensus_string
    return donor_genome


def edit_distance_matrix(ref, donor):
    """
    Computes the edit distance matrix between the donor and reference

    This algorithm makes substitutions, insertions, and deletions all equal.
    Does that strike you as making biological sense? You might try changing the cost of
    deletions and insertions vs snps.
    :param ref: reference genome (as an ACTG string)
    :param donor: donor genome guess (as an ACTG string)
    :return: complete (len(ref) + 1) x (len(donor) + 1) matrix computing all changes
    """

    output_matrix = np.zeros((len(ref), len(donor)),dtype=int)
    # print len(ref), len(donor)
    # print output_matrix
    # This is a very fast and memory-efficient way to allocate a matrix
    for i in range(len(ref)):
        output_matrix[i, 0] = i
    for j in range(len(donor)):
        output_matrix[0, j] = j
    for j in range(1, len(donor)):
        for i in range(1, len(ref)):  # Big opportunities for improvement right here.
            deletion = output_matrix[i - 1, j] + 1 #1. Make it bigger or smaller?
            insertion = output_matrix[i, j - 1] + 1 # 1
            identity = output_matrix[i - 1, j - 1] if ref[i] == donor[j] else np.inf
            substitution = output_matrix[i - 1, j - 1] + 1 if ref[i] != donor[j] else np.inf
            output_matrix[i, j] = min(insertion, deletion, identity, substitution)
    return output_matrix

def alignment_matrix(ref,donor):
    import numpy as np
    output_matrix = np.zeros((len(ref), len(donor)),dtype=int)
    mismatch=-2
    match=2
    inse=-2
    dele=-2
    PM=np.zeros((len(ref)-1,len(donor)-1),dtype=int) # Pointer matrix (backtrack_matrix). right= 4, down= 3, Diagonal= 2
    Right=4
    Down=3
    Diagonal=2
    Start=1
    for i in range(len(ref)):
        output_matrix[i, 0] = dele*i
    for j in range(len(donor)):
        output_matrix[0, j] = inse*j
    for j in range(1, len(donor)):
        for i in range(1, len(ref)):  # Big opportunities for improvement right here.
            deletion = output_matrix[i - 1, j] + dele #1. Make it bigger or smaller?
            insertion = output_matrix[i, j - 1] + inse # 1
            identity = output_matrix[i - 1, j - 1] + match  if ref[i] == donor[j] else -np.inf
            substitution = output_matrix[i - 1, j - 1] + mismatch if ref[i] != donor[j] else -np.inf
            output_matrix[i, j] = max(insertion, deletion, identity, substitution,0)
            if output_matrix[i, j] ==identity or output_matrix[i, j] ==substitution:
                PM[i-1,j-1]=Diagonal
            elif output_matrix[i, j] ==deletion:
                PM[i-1,j-1]=Down
            elif output_matrix[i, j] ==insertion:#insertion
                PM[i-1,j-1]=Right
            else:#0, start position
                PM[i-1,j-1]=Start
    return output_matrix,PM              


def identify_changes(ref, donor, offset):
    """
    Performs a backtrace-based re-alignment of the donor to the reference and identifies
    SNPS, Insertions, and Deletions.
    Note that if you let either sequence get too large (more than a few thousand), you will
    run into memory issues.

    :param ref: reference sequence (ATCG string)
    :param donor: donor sequence (ATCG string)
    :param offset: The starting location in the genome.
    :return: SNPs, Inserstions, and Deletions
    """
    # print offset
    ref = '${}'.format(ref.strip())  # add a dollar sign to the front
    donor = '${}'.format(donor.strip())
    edit_matrix = edit_distance_matrix(ref=ref, donor=donor)
    #print edit_matrix
    #current_row = len(ref) - 1
    #current_column = len(donor) - 1
    last_col=edit_matrix[:,-1]
    last_row=edit_matrix[-1,:]
    if last_row.min()<=last_col.min():
        current_row = len(ref) - 1
        current_column = last_row.argmin()
    else:
        current_row = last_col.argmin()
        current_column =len(donor) - 1
    changes = []
    while current_row > 0 or current_column > 0:
        if current_row == 0:
            pvs_row = -np.inf
        else:
            pvs_row = current_row - 1

        if current_column == 0:
            pvs_column = -np.inf
        else:
            pvs_column = current_column - 1

        try:
            insertion_dist = edit_matrix[current_row, pvs_column]
        except IndexError:
            insertion_dist = np.inf

        try:
            deletion_dist = edit_matrix[pvs_row, current_column]
        except IndexError:
            deletion_dist = np.inf

        try:
            if ref[current_row] == donor[current_column]:
                identity_dist = edit_matrix[pvs_row, pvs_column]
            else:
                identity_dist = np.inf

            if ref[current_row] != donor[current_column]:
                substitution_dist = edit_matrix[pvs_row, pvs_column]
            else:
                substitution_dist = np.inf
        except (TypeError, IndexError) as e:
            identity_dist = np.inf
            substitution_dist = np.inf

        min_dist = min(insertion_dist, deletion_dist, identity_dist, substitution_dist)

        ref_index = current_row + offset - 1
        if min_dist == identity_dist:
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == substitution_dist:
            changes.append(['SNP', ref[current_row], donor[current_column], ref_index])
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == insertion_dist:
            if len(changes) > 0 and changes[-1][0] == 'INS' and changes[-1][-1] == ref_index + 1:
                changes[-1][1] = donor[current_column] + changes[-1][1]
            else:
                #print current_column
                #print donor[100-current_column]
                changes.append(['INS', donor[current_column], ref_index + 1])
            current_column = pvs_column
        elif min_dist == deletion_dist:
            if len(changes) > 0 and changes[-1][0] == 'DEL' and changes[-1][-1] == ref_index + 1:
                changes[-1] = ['DEL', ref[current_row] + changes[-1][1], ref_index]
            else:
                changes.append(['DEL', ref[current_row], ref_index])
            current_row = pvs_row
        else:
            raise ValueError
    changes = sorted(changes, key=lambda change: change[-1])

    if len(changes) >=10:  # When there are so many changes in 100 bp, it means there is a big trunk of insertion. 
        for i in range(len(changes)):
            if changes[i][0]=='INS':
                change=changes[i]
                changes=[]
                changes.append(change)
                break
    print str(changes)
    return changes

def alignment(ref, donor, offset): #new version of itendify_changes
    import numpy as np
    # print offset
    ref = '${}'.format(ref.strip())  # add a dollar sign to the front
    donor = '${}'.format(donor.strip())
    align_matrix,PM = alignment_matrix(ref=ref, donor=donor)
    last_col=align_matrix[:,-1]
    last_row=align_matrix[-1,:]
    aligned_ref=''
    aligned_donor=''
    if last_row.max()>=last_col.max():
        current_row = len(ref) - 1
        current_column = last_row.argmax()
        if last_row.argmax()<len(donor)-1:
            aligned_donor=donor[-(len(donor)-1-last_row.argmax()):]
            aligned_ref='-'*(len(donor)-1-last_row.argmax())
    else:
        current_row = last_col.argmax()
        current_column =len(donor) - 1
        aligned_donor='-'*(len(ref)-1-last_col.argmax())
        aligned_ref=ref[-(len(ref)-1-last_col.argmax()):]
        
    changes = []
    Right=4
    Down=3
    Diagonal=2
    Start=1
    while current_row >= 1 and current_column >= 1:
        #print (current_row,current_column)
        ref_index = current_row + offset - 1 
        
        if PM[current_row-1,current_column-1]==Right:#insertion
            aligned_ref='-'+aligned_ref
            aligned_donor=donor[current_column]+aligned_donor
            if len(changes) > 0 and changes[-1][0] == 'INS' and changes[-1][-1] == ref_index + 1:
                changes[-1][1] = donor[current_column] + changes[-1][1]
            else:
                changes.append(['INS', donor[current_column], ref_index + 1])
            current_column=current_column-1
            
        elif PM[current_row-1,current_column-1]==Down:#deletion
            aligned_ref=ref[current_row]+aligned_ref
            aligned_donor='-'+aligned_donor
            if len(changes) > 0 and changes[-1][0] == 'DEL' and changes[-1][-1] == ref_index + 1:
                changes[-1] = ['DEL', ref[current_row] + changes[-1][1], ref_index]
            else:
                changes.append(['DEL', ref[current_row], ref_index])
            current_row=current_row-1
            
        elif PM[current_row-1,current_column-1]==Diagonal and ref[current_row]==donor[current_column]: #match
            aligned_ref=ref[current_row]+aligned_ref
            aligned_donor=donor[current_column]+aligned_donor
            current_row -=1
            current_column -=1
        elif PM[current_row-1,current_column-1]==Diagonal and ref[current_row]!=donor[current_column]: #mismatch
            aligned_ref=ref[current_row]+aligned_ref
            aligned_donor=donor[current_column]+aligned_donor
            changes.append(['SNP', ref[current_row], donor[current_column], ref_index])
            current_row -=1
            current_column -=1
        else: #start
            break

        ####
#     print(align_matrix[current_row,current_column])
#     print(current_row)
#     print(current_column)
    if current_row-1 > current_column-1:
        aligned_donor='-'*(current_row-current_column)+donor[1:current_column+1] + aligned_donor
        aligned_ref= ref[1:current_row+1]+aligned_ref
    if current_row-1 < current_column-1:
        aligned_ref='-'*(current_column-current_row)+ref[1:current_row+1]+aligned_ref
        aligned_donor=donor[1:current_column+1]+aligned_donor
#     print(aligned_ref)
#     print(aligned_donor)

    changes = sorted(changes, key=lambda change: change[-1])
    if len(changes) >=10:  # When there are so many changes in 100 bp, it means there is a big trunk of insertion. 
        for i in range(len(changes)):
            if changes[i][0]=='INS':
                change=changes[i]
                changes=[]
                changes.append(change)
                break
    print changes
    return changes

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

def determine_STR(STRinRef, aligned_reads, alignments):
    import numpy as np
    STRlist=[]
    for STR in STRinRef:
        if len(STR[1])<44:
            reads=[]
            for i in range(len(alignments)):
                if STR[2]-20<alignments[i][0] <STR[2]-5:
                    reads.append(aligned_reads[i][0])
                if STR[2]-20<alignments[i][1] <STR[2]-5:
                    reads.append(aligned_reads[i][1])
            copynumber=[]
            STRunit=''
            for read in reads:
                L=determine_STR_inRef(read)
                if len(L)>0:
                    copynumber.append(L[0][-1])
                    unitlength=int(len(L[0][-3])/L[0][-1])
                    STRunit=L[0][-3][0:unitlength]
            if len(copynumber)>0:
            	cn=np.array(copynumber)
            	determined_copynumber=cn.max() 
            	STRlist.append(['STR',STRunit*determined_copynumber,STR[2],determined_copynumber])
        else:
            STRlist.append(STR)
    return STRlist

def consensus(ref, aligned_reads):
    """
    Identifies a consensus sequence by calling the most commmon base at each location
    in the reference.
    :param ref: reference string
    :param aligned_reads: the list of reads.
    :return: The most common base found at each position in the reads (i.e. the consensus string)
    """
    consensus_string = ''
    padded_reads = [aligned_read + ' ' * (len(ref) - len(aligned_read)) for aligned_read in aligned_reads]
    # The reads are padded with spaces so they are equal in length to the reference
    for i in range(len(ref)):
        base_count = defaultdict(float)
        ref_base = ref[i]
        base_count[ref_base] += 1.1  # If we only have a single read covering a region, we favor the reference.
        read_bases = [padded_read[i] for padded_read in padded_reads if padded_read[i] not in '. ']
        # Spaces and dots (representing the distance between paired ends) do not count as DNA bases
        for base in read_bases:
            base_count[base] += 1
        consensus_base = max(base_count.iterkeys(), key=(lambda key: base_count[key]))
        # The above line chooses (a) key with maximum value in the read_bases dictionary.
        consensus_string += consensus_base
    return consensus_string # this consensus_string has the same length as reference string



if __name__ == "__main__":

    # ### Testing code for Smith-Waterman Algorithm
    # print edit_distance_matrix('$PRETTY', '$PRTTEIN')
    # identify_changes('PRETTY', 'PRTTEIN', offset=0)
    # identify_changes(ref='ACACCC', donor='ATACCCGGG', offset=0)
    # identify_changes(ref='ATACCCGGG', donor='ACACCC', offset=0)
    # identify_changes(ref='ACACCC', donor='GGGATACCC', offset=0)
    # identify_changes(ref='ACA', donor='AGA', offset=0)
    # identify_changes(ref='ACA', donor='ACGTA', offset=0)
    # identify_changes(ref='TTACCGTGCAAGCG', donor='GCACCCAAGTTCG', offset=0)
    # ### /Testing Code
    #
    genome_name ='hw2undergrad_E_2' #'practice_W_3'###'hw1_W_2'###'practice_W_1'
    input_folder = '../data/{}'.format(genome_name)
    chr_name = '{}_chr_1'.format(genome_name)
    reads_fn_end = 'reads_{}.txt'.format(chr_name)
    reads_fn = join(input_folder, reads_fn_end)
    ref_fn_end = 'ref_{}.txt'.format(chr_name)
    ref_fn = join(input_folder, ref_fn_end)
    start = time.clock()
    input_fn = join(input_folder, 'aligned_{}.txt'.format(chr_name))
    #snps, insertions, deletions = generate_pileup(input_fn)

    STRinRef=pickle.load(open(join(input_folder,genome_name+'_STR.pkl'), 'rb'))
    (aligned_reads, alignments)=pickle.load(open(join(input_folder, genome_name+'_aligned_reads_and_coordinates.pkl'),'rb'))
    STRs=determine_STR(STRinRef, aligned_reads, alignments)

    '''
    output_fn = join(input_folder, 'changes_{}.txt'.format(chr_name))
    zip_fn = join(input_folder, 'changes_{}.zip'.format(chr_name))
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x[1:]]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x[1:]]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x[1:]]) + '\n')
        output_file.write('>STR\n')  
        for x in STRs:
        	output_file.write(x[1]+','+ str(x[2]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
    '''
    output_STR=join(input_folder, 'STR.txt'.format(chr_name))
    with open(output_STR,'w') as output_f:
    	output_f.write('>STR\n') 
    	for x in STRs:
        	output_f.write(x[1]+','+ str(x[2]) + '\n')
