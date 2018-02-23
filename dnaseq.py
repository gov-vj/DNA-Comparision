#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        #raise Exception("Not implemented!")
        self.dic={}
        for k,v in pairs:
          self.put(k,v)
    # Associates the value v with the key k.
    def put(self, k, v):
        #raise Exception("Not implemented!")
        if k in self.dic:
          self.dic[k].append(v)
        else:
          self.dic[k]=[v]
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        #raise Exception("Not implemented!")
        if k in self.dic: return self.dic[k]
        return []

# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    #use this method as we iterate over sequence 'b'	
	#return: sub seq, its index position in the seq and the hash value of sub seq
    sub_seq=""
    for i in range(k):
      sub_seq+=next(seq)
    rh=RollingHash(sub_seq)
    yield (sub_seq,rh.current_hash(),0)
    index=1
    for i in seq:
      temp=rh.slide(sub_seq[0],i)
      sub_seq=sub_seq[1:]+i                      
      yield (sub_seq,temp,index)                  
      index+=1
      

# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    #same as subsequenceHashes but returns one sub seq every m nucleotides
	#use this method as we iterate over sequence 'a'
    i=0
    for ch in seq:
        if i==0:
            index=seq.pos-1
            sub_seq=ch
            for _ in range(k-1):
                sub_seq+=next(seq)
                i=(i+1)%m
            rh=RollingHash(sub_seq)
            yield sub_seq,rh.current_hash(),index
        i=(i+1)%m
    

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    a_d=Multidict()  #slightly modified dictionary
    for sub_seq,hash,index in intervalSubsequenceHashes(a,k,m):
      a_d.put(hash,(sub_seq,index))
    for sub_seq,hash,index in subsequenceHashes(b,k):
      x=a_d.get(hash)
      if len(x)>0:
        for sub_seq_a,index_a in x:
            if sub_seq_a==sub_seq: yield index_a,index

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0])
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
