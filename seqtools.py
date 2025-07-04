# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 09:57:22 2025

@author: giera
"""
from  OOP_learning import Sequence, disk_cache
import argparse
import os
import pickle
import re
import zlib
import zipfile
from datetime import datetime, date

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True, help="Input file path")
parser.add_argument("-o","--output", required = True, help="Output path for your pickle file")
parser.add_argument("--clear_cache", default=False ,choices=[True, False] ,help="Clears out cache folder. True = Yes, False = No")

args = parser.parse_args()

fasta_formats = ('fasta', 'fa', 'fna', 'ffn', 'faa', 'frn', 'mpfa')

file = 'test.fa'

zip_path = os.path.join(args.output, 'old_cache.zip')

def main():
    
    # checks if file is fasta or fastq format
    if str(args.input).endswith('.fastq'):
        records = Sequence.parse_fastq_file(args.input)
    elif args.input.split('.')[-1] in fasta_formats:
        records = Sequence.parse_fasta_file(args.input)
        
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Setting cache
    cached_kmer = disk_cache(cache_dir = args.output)(Sequence.expensive_kmer)
    
    # Creating one dict for whole file
    all_kmers ={}
    for seq in records:
        kmer_dict = cached_kmer(seq)
        all_kmers[seq.header] = kmer_dict
    
    # Writing out file to specyfic file
    all_kmers_path = os.path.normpath(os.path.join(args.output, 'all_kmers.pkl'))
    with open(all_kmers_path, 'wb') as out_f:
        pickle.dump(all_kmers, out_f)
        
    print(f"Wrote {len(all_kmers)} records to {all_kmers_path} (with file combining them all named: all_kmers.pkl)")
    
    if args.clear_cache:
        for file in os.listdir(args.output):
            to_remove = re.match(".*pkl", file).string
            os.remove(f'{args.output}/{to_remove}')

    # loading files and creating old_files list
    old_files = []
    for file in os.listdir(args.output):
        # defining names of files to later loop through
        filepath = os.path.normpath(os.path.join(args.output, file))
        
        # calculatnig file age
        mtime = os.path.getmtime(filepath)
        file_age = datetime.fromtimestamp(mtime).date()
        now = datetime.now().date()
        age = now - file_age
        if age.days >= 7:
            old_files.append(filepath)
    
    with zipfile.ZipFile(zip_path, 'w') as archive:
        for filename in old_files:
            archive.write(filename, arcname=os.path.basename(filename))
    
    for filepath in old_files:
        os.remove(filepath)
    
if __name__ == "__main__":
    main()