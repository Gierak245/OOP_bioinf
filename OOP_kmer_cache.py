from hashlib import sha256
import os
import pickle
from functools import wraps
import re

def disk_cache(cache_dir:str ="cache"):
    """
    Creates cache folder and decorator used for creating .pkl files from Sequence class.

    Parameters
    ----------
    cache_dir : str, optional
        Path to cache folder. The default is "cache".

    Returns
    -------
    Callable[[Callable], Callable]
        A decorator that wraps a function taking a Sequence and returns its cached result.

    """
    # checks if cache directory exists
    if os.path.exists(cache_dir) != True:
        os.makedirs(cache_dir)
        
            
    def decorator(func):
        
        @wraps(func)
        def wrapper(arg):
            # Checking if given arg (Class) has a sequence string
            if hasattr(arg, 'sequence_string'):
                
                # hetting sequence from Class
                seq_text = arg.sequence_string
                key = sha256(seq_text.encode('utf‑8')).hexdigest()  # creating unique key value for .pkl file
            
            # creating path to pickle file in cache folder
            key_path = os.path.normpath(os.path.join(cache_dir, f'{key}.pkl'))
            
            # Loading, or creating and than loading, a .pkl file.
            if os.path.exists(key_path) == True:
                print(f"Loading from cache data for {arg.header}")
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
    """
    Represents a biological sequence parsed from FASTA or FASTQ.

    Attributes
    ----------
    header : str
        The record identifier.
    sequence_string : str
        The nucleotide (or amino‑acid) string.
    quality : Optional[str]
        Quality scores (only for FASTQ), otherwise None.
    """

    def __init__(self, header, sequence_string, quality):
        """Initialize with header, sequence, and optional quality."""
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
        """
        Parse a single FASTA record into a Sequence object.

        Parameters
        ----------
        text : str
            Two‑line FASTA record, including the leading '>'.

        Raises
        ------
        ValueError
            If fasta sequence lines are empty

        Returns
        -------
        Sequence
            A new Sequence instance with header, sequence_string, and quality == None attributes.

        """
        # creating temp headers and seq
        header = None
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
        """
        Loads fastq sequence given in 'text' parameter.

        Parameters
        ----------
        text : str
            Fastq sequence with header and quality string.

        Raises
        ------
        ValueError
            If fastq file is too short (less than 4 lines) or fastq sequence is empty.

        Returns
        -------
        Sequence
            A new Sequence instance with header, sequence_string, and quality attributes.

        """
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
    
    @classmethod
    def parse_fasta_file(cls, filepath: str):
        """
        Create multiple Sequence classes from fasta file

        Parameters
        ----------
        filepath : str
            Path to fasta file

        Yields
        ------
        Sequence
            Each record in the input file as a Sequence instance.

        """
        with open(filepath, 'r') as f:
            
            header, seq_lines = None, []
            
            for line in f.readlines():
                line = line.rstrip()
                if line.startswith('>'):
                    if header:
                        yield cls(header, ''.join(seq_lines), None)                
                    header = line[1:]
                    
                    seq_lines = []
                    
                else:
                    seq_lines.append(line)
            if header and seq_lines:
                yield cls(header, ''.join(seq_lines), None)   
            
    @classmethod
    def parse_fastq_file(cls, filepath: str):
        """
        Create multiple Sequence classes from fastq file

        Parameters
        ----------
        filepath : str
            Path to fasta file

         Yields
         ------
         Sequence
             Each record in the input file as a Sequence instance.

        """
        with open(filepath, 'r') as f:
            
            file = f.readlines()
            header, sequence_string, quality = None, None, None
            
            i = 0
            
            while i <= len(file) - 4:
                
                header = file[i]
                sequence_string = file[i + 1]
                quality = file[i + 3]
                
                yield cls(header, sequence_string, quality)
                i += 4
            
                    
   # @disk_cache()
    def expensive_kmer(self):
        """
        Computationally expensive function (k-mer in this script)

        Returns
        -------
        k_mers : dict[str:list]
            
        Dictionary with k (length of substrings) as keys and list of substrings as values
        
        Examples
        --------
        >>> seq = Sequence.from_fasta(">s1\nATGC\n")
        >>> seq.expensive_kmer()
        {1: ['A','T','G','C'], 2: ['AT','TG','GC'], …}

        """
        k = len(self.sequence_string)
        k_mers = {}

        for i in range(1, k + 1):
            temp_list = []
            for j in range(0, k + 1):
                if j + i <= k:
                    temp_list.append(self.sequence_string[j : j + i])
                k_mers[i] = temp_list
        
        return k_mers
        
