# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 19:22:59 2025

@author: giera
"""

from hashlib import sha256
import os
import pickle
from functools import wraps

def disk_cache(cache_dir="cache"):
    
    if os.path.exists(cache_dir) != True:
        os.makedirs(cache_dir)
        
            
    def decorator(func):
        
        @wraps(func)
        def wrapper(arg):
            # Checking if given arg is a sequence string
            if hasattr(arg, 'sequence_string'):
                
                seq_text = arg.sequence_string
                key = sha256(seq_text.encode('utf‑8')).hexdigest()
            
            # path to pickle file in cache folder
            key_path = os.path.normpath(os.path.join(cache_dir, f'{key}.pkl'))
            
            # Loading or creating and than loading a pkl file.
            if os.path.exists(key_path) == True:
                print("Loading from cache")
                with open(key_path, 'rb') as pickle_file:  
                    kmer_db = pickle.load(pickle_file)
                
                return kmer_db
            
            else:
                # return function output
                result = func(arg)
                
                # writing result into a file
                with open(key_path, 'wb') as pickle_file:
                    pickle.dump(result, pickle_file)
                
                return result
        return wrapper
    return decorator



class Sequence:
    
    def __init__(self, header, sequence_string, quality):
        self.header = header
        self.sequence_string = sequence_string
        self.quality = quality
    
    def __repr__(self):
        return f'''Sequence header: {self.header}
Sequence: {self.sequence_string}
Length: {len(self.sequence_string)}
'''
    
    # Classmethod creating a Sequence Class from fasta file
    @classmethod
    def from_fasta(cls, text: str):
        
        # creating temp headers and seq
        header = ''
        sequence_string = ''
        
        # loop to get header and sequence from file
        
        for line in text.strip().splitlines():
            if line.startswith(">"):
                header = line.strip()[1:]
            else:
               sequence_string += line.strip()
               
        if sequence_string == '':
            raise ValueError("Sequence is empty. Please check your FASTA file.")
        
        return cls(header, sequence_string, quality = None)
    
    # Classmethod creating a Sequence Class from fastq file
    @classmethod
    def from_fastq(cls, text: str):
        
        # creating headers and sequence_string
        
        splitted_text = text.strip().splitlines()
        
        if len(splitted_text) < 4:
            raise ValueError("Your Fastq file is too short (less than 4 lines). Please check your file")
        
        header = splitted_text[0][1:]
        sequence_string = splitted_text[1]
        quality_string = splitted_text[3]
        
        if quality_string == '':
            raise ValueError("Quality string is empty. Please check your FASTQ file.")
        

        return cls(header, sequence_string, quality_string)
    
                    
    @disk_cache()
    def expensive(self):
        
        k = len(self.sequence_string)
        k_mers = {}

        for i in range(1, k + 1):
            temp_list = []
            for j in range(0, k + 1):
                if j + i <= k:
                    temp_list.append(self.sequence_string[j : j + i])
                k_mers[i] = temp_list
        
        return k_mers
        
    

seq = Sequence.from_fasta(">s1\nATGC\n")
# First call: should compute and create "cache/<hash>.pkl"
db1 = seq.expensive()
# Second call: should load from disk, not recompute
db2 = seq.expensive()

fasta = """>seq1
        AGCTTGA
        GCTTATT
        CCCAAGT"""

fastq = """@seq2
GGTTAACC
+
IIIIIIII
"""

# for line in text.strip().splitlines(): print(line.startswith(">"))

test_fasta = Sequence.from_fasta(fasta)
test_fasta

test_fastq = Sequence.from_fastq(fastq)
test_fastq


cache_patch = "C:/Users/giera/OneDrive/Dokumenty/Python_Scripts/cache"

# creating key and assiging sequence 
key = sha256(test_fasta.sequence_string.encode()).hexdigest()
db = {}
db[key] = test_fasta.sequence_string

# Checking if catche path is present 
os.path.exists(cache_patch)


# path_to_pckle = os.path.normpath(os.path.join(cache_patch, f'{key}.pkl'))
# print(path_to_pckle)
# print(os.path.exists(path_to_pckle))


# with open('test.pkl', 'rb') as f:
#     db1 = pickle.load(f)
#     print(db1)

# seq = "ATTGCC"

# def expensive(seq):
    
#     k = len(seq)
#     k_mers = {}

#     for i in range(1, k + 1):
#         temp_list = []
#         for j in range(0, k + 1):
#             if j + i <= k:
#                 temp_list.append(seq[j : j + i])
#             k_mers[i] = temp_list
    
#     return k_mers

# k_mers = expensive(test_fasta.sequence_string)

# sequence_key = sha256(test_fasta.sequence_string.encode()).hexdigest()

# with open(f'{sequence_key}.pkl', 'wb') as pickle_file:  
#     pickle_test = pickle.dump(sequence_key, pickle_file)
# pickle_file.close()
