import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from CS_CM122.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref
import pandas as pd


#####################################
def BWT(T):
    T=T+'$'
    s=[]
    s.append(T)
    for i in range(len(T)-1):
        T=T[1:]+T[0]
        s.append(T)
    s=sorted(s)
    transformed=''
    for t in s:
        transformed+=t[-1]
    return transformed

def invBWT(BWT): #inverse BW transformation  of sequence BWT. return the first occurrence of each character 
    d={'$':0,'A':0,'C':0,'G':0,'T':0}
    for s in BWT:
        d[s]+=1
    FstOcc=[0,1,1+d['A'],1+d['A']+d['C'],1+d['A']+d['C']+d['G']] #$,'A','C','G','T'
    return FstOcc,d

def counttable(BWT):# build the count table of BWT:
    import numpy as np
    L=len(BWT)
    table=np.zeros((L+1,5),dtype=int)
    d={'$':0,'A':1,'C':2,'G':3,'T':4}
    count=np.array((0,0,0,0,0))

    for i in range(L):
        count[d[BWT[i]]]+=1
        table[i+1]=count

    return table

def buildindex(T): #Build original index of T for BWT, which was sorted with the partial suffix array
    suffixlist=[]
    T=T+'$'
    for i in range(len(T)):
        suffixlist.append((T[i:],i))
    sortedsuffix=sorted(suffixlist)
    import pandas as pd
    ss=pd.DataFrame(sortedsuffix)
    index=ss.loc[:,1]
    return index

def BWT_algorithm(paired_end_reads, ref):
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()

    bwt=BWT(ref)
    FirstOcc,numOfNu=invBWT(bwt)
    counttb=counttable(bwt)
    IND=buildindex(ref)

    def BWT_algorithm_readmapping(read,Perfect_match_location):
        d={'$':0,'A':1,'C':2,'G':3,'T':4}
        index=d[read[-1]]
        top=FirstOcc[index]
        if index==4:
            bottom=FirstOcc[index]+numOfNu['T']-1
        else:
            bottom=FirstOcc[index+1]-1
        for i in range(1,len(read)):
            index=d[read[-1-i]]
            top=FirstOcc[index]+counttb[top,index]
            bottom=FirstOcc[index]+counttb[bottom+1,index]-1
            if bottom < top: # no match
                break
        if bottom >top: # find matches
            print('multiple matches:')
            for i in range(top, bottom+1):
                print(IND.loc[i])
                #Perfect_match_location.append(IND.loc[i]) 
        elif bottom==top:
            Perfect_match_location.append(IND.loc[top])


    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))
        for read in read_pair:
            Perfect_forward_match_location=[]
            BWT_algorithm_readmapping(read,Perfect_forward_match_location)

            reversed_read = read[::-1]
            Perfect_reverse_match_location=[]
            BWT_algorithm_readmapping(reversed_read,Perfect_reverse_match_location)

            if Perfect_forward_match_location!=[] and Perfect_reverse_match_location!=[]:
                print("find reads have both forward and reverse mapping")
            elif Perfect_forward_match_location!=[]:
                read_alignment_locations.append(Perfect_forward_match_location[0])
                output_read_pair.append(read)
            elif Perfect_reverse_match_location!=[]:
                read_alignment_locations.append(Perfect_reverse_match_location[0])
                output_read_pair.append(read)

        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)
    return all_read_alignment_locations, output_read_pairs



def hash_algorithm(paired_end_reads, ref):
    def buildhash(T):
        k=12#len(read)//4
        d={}
        for i in range(len(T)-k+1):
            kmer=T[i:i+k]
            if kmer not in d:
                d[kmer]=[i]
            else:
                d[kmer].append(i)
        return d

    def hash_read_mapping(read,T,read_alignment_locations,output_read_pair,hashtable):
        def consecutive(s,t):# check whether the difference nucleotides between two reads are consecutive (as long as two consecutive nucleotides are different)
            l=[]
            for i in range(len(s)):
                if s[i]!=t[i]:
                    l.append(i)
            if len(l)==2:
                if l[1]-l[0]==1:
                    return True
            else: # len(l)==3:
                if l[1]-l[0]==1 or l[2]-l[1]==1:
                    return True
                else:
                    return False
        def addreadsifnotconsecutive(d,s,t,i,read,read_alignment_locations,output_read_pair):# d: Hamming distance, s and t are strings, i is the index of reference genome
            if d>=2 and d<=3: #putatively consecutive
                if not consecutive(s,t):
                    read_alignment_locations.append(i)
                    output_read_pair.append(read)
                else: 
                    read_alignment_locations.append(i)
                    output_read_pair.append(read)

        def dH(s,t): # Hamming distance (if length of s and t are not equal, return -1)
            if len(s)!=len(t):
                return -1
            d=0 #initialize distance =0
            for i in range(0,len(s)):
                if s[i] != t[i]:
                    d+=1
            return d
        r1=read[0:12]
        r2=read[13:25]
        r3=read[26:38]
        r4=read[38:]
        L=len(T)
        if r1 in hashtable:
            for i in hashtable[r1]:#compare read[12:]
                d=dH(T[i+12:min(i+50,L)],read[12:])
                if d>=0 and d<=3:
                    if d>=2 and d<=3: #putatively consecutive
                        if not consecutive(T[i+12:min(i+50,L)],read[12:]):
                            read_alignment_locations.append(i)
                            output_read_pair.append(read)
                    else: 
                        read_alignment_locations.append(i)
                        output_read_pair.append(read)
        elif r4 in hashtable:
            for i in hashtable[r4]:#compare read[0:38] 
                d=dH(T[max(0,i-38):i],read[0:38])
                if d>=0 and d<=3:
                    if d>=2 and d<=3: #putatively consecutive
                        if not consecutive(T[max(0,i-38):i],read[0:38]):
                            read_alignment_locations.append(i-38)
                            output_read_pair.append(read)
                    else:
                        read_alignment_locations.append(i-38)
                        output_read_pair.append(read)
        elif r2 in hashtable:
            #compare read[0:13] and read[25:]
            for i in hashtable[r2]:
                d=dH(T[max(0,i-13):i+37],read)
                if d>=0 and d<=3:
                    if d>=2 and d<=3: #putatively consecutive
                        if not consecutive(T[max(0,i-13):i+37],read):
                            read_alignment_locations.append(i-13)
                            output_read_pair.append(read)
                    else:
                        read_alignment_locations.append(i-13)
                        output_read_pair.append(read)
        elif r3 in hashtable:
            #compare read[0:26] and read[38:]
            for i in hashtable[r3]:
                d=dH(T[max(0,i-26):i+24],read)
                if d>=0 and d<=3:
                    if d>=2 and d<=3: #putatively consecutive
                        if not consecutive(T[max(0,i-26):i+24],read):
                            read_alignment_locations.append(i-26)
                            output_read_pair.append(read)
                    else:
                        read_alignment_locations.append(i-26)
                        output_read_pair.append(read)

    hashtable=buildhash(ref)
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    unmappedcount=0
    start = time.clock()
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))
        for read in read_pair:
            hash_read_mapping(read,ref,read_alignment_locations,output_read_pair,hashtable)
            reversed_read = read[::-1]
            hash_read_mapping(reversed_read,ref,read_alignment_locations,output_read_pair,hashtable)
        if read_alignment_locations==[]:
            unmappedcount+=2
        elif len(read_alignment_locations)==1:
            unmappedcount+=1
        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)
    print('Number of reads unmapped: %s' %unmappedcount)
    return all_read_alignment_locations, output_read_pairs
#########################################################
###################################################

def trivial_algorithm(paired_end_reads, ref):
    """
    This is a functional aligner, but it's a huge simplification that
    generate a LOT of potential bugs.  It's also very slow.

    Read the spec carefully; consider how the paired-end reads are
    generated, and ideally, write your own algorithm
    instead of trying to tweak this one (which isn't very good).

    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()
    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))
        for read in read_pair:
            min_mismatches = len(read) + 1
            min_mismatch_location = -1
            for i in range(len(ref) - len(read)):
                mismatches = [1 if read[j] != ref[i + j] else 0 for j in range(len(read))]
                n_mismatches = sum(mismatches)
                # The above line should be familiar to Python users, but bears  some explanation for
                # people who are getting started with it. The "mismatches = ..." line
                # is called a "list comprehension. Basically, this is a short way of writing the loop:
                #
                # n_mismatches = 0
                # for j in range(len(read)):
                # if read[j] != ref[i+j]:
                #         n_mismatches += 1
                #
                # The first line creates a list which has a 1 for every mismatch and a 0 for every match.
                # The second line sums the list created by the first line, which counts the number of mismatches.
                if n_mismatches < min_mismatches:
                    min_mismatches = n_mismatches
                    min_mismatch_location = i

            reversed_read = read[::-1]
            for i in range(len(ref) - 50):
                mismatches = [1 if reversed_read[j] != ref[i + j] else 0 for j in range(len(read))]
                n_mismatches = sum(mismatches)
                if n_mismatches < min_mismatches:
                    min_mismatches = n_mismatches
                    min_mismatch_location = i
                    read = reversed_read
            read_alignment_locations.append(min_mismatch_location)
            output_read_pair.append(read)
            # # Note that there are some huge potential problems here.

        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)
    return all_read_alignment_locations, output_read_pairs


if __name__ == "__main__":
    data_folder = 'hw1_W_2'#'practice_W_1'
    input_folder = join('../data/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.

    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    #alignments, reads = trivial_algorithm(input_reads, reference)
    #alignments, reads = BWT_algorithm(input_reads, reference)
    alignments, reads = hash_algorithm(input_reads, reference)
    print(alignments)
    print(reads)
    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
