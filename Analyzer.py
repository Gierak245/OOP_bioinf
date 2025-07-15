import sys
import os
import re
from abc import ABC, ABCMeta, abstractmethod
from typing import Any, List, Dict, Type, Iterable
sys.path.append(os.path.abspath(os.path.join(__file__, '..', '..', 'OOP')))
from OOP_kmer_cache import Sequence

class AnalyzerMeta(ABCMeta):
    
    registry: Dict[str, Type['Analyzer']] = {}
    
    def __init__(cls, name, base, namespace):
        super().__init__(name, base, namespace)
            
        if name != 'Analyzer':
            AnalyzerMeta.registry[name] = cls

        

class Analyzer(metaclass=AnalyzerMeta):
    
    '''
    Common interface/blueprint for all sequence analyses
    '''
    
    @abstractmethod
    def run(self, seq: Sequence) -> Any:
        '''
        Perform analysis on sequence instance.
        
        Parameters
        ----------
        
        seq: Sequence
            A sequence object to analyze
            
        Returns
        -------
        Any
            The analysis result (e.g., GC contetnt as float, dict of motif counts)
        '''
        pass

class GCContentAnalyzer(Analyzer):
    
    """
    Get GC value in seq.sequence_string.

    Returns
    -------
    float
        float percent value of GC content in sequence 
    """
    
    def run(self, seq: Sequence) -> float:
        """
        Compute the GC content of a sequence.
    
        Parameters
        ----------
        seq : Sequence
            Sequence object whose `sequence_string` will be used.
    
        Returns
        -------
        float
            Fraction of bases that are G or C (0.0–1.0).
    
        Examples
        --------
        >>> analyzer = GCContentAnalyzer()
        >>> analyzer.run(Sequence("s1","GCGT",None))
        0.75
        """
        
        nucl_string = seq.sequence_string.upper()
        
        gc_content = (nucl_string.count('C') + nucl_string.count('G'))/len(nucl_string)
        
        return gc_content 
        

class MotifSearchAnalyzer(Analyzer): 
    
    
    """
    Find all occurrences of self.motif in seq.sequence_string.

    Returns
    -------
    List[int]
        1-based start positions of each motif occurrence.
    
    Examples
    --------
    >>> analyzer = MotifSearchAnalyzer()
    >>> analyzer.run(Sequence("s1","ATGTTGGGCATG",None))
    [1, 10]
    """
    
    def __init__(self, motif: str):
        self.motif = motif
        
    def run(self, seq: Sequence) -> List[int]:
        
        all_motifs = [occ.span()[0] + 1 for occ in re.finditer(self.motif, seq.sequence_string)]
        if all_motifs:
            return all_motifs
        else:
            return 0

class PluginManager:
    
    """
    Manages discovery and execution of Analyzer plugins.

    This class allows you to register Analyzer subclasses, specify the order
    in which they should run, and execute them on `Sequence` instances.

    Attributes
    ----------
    name_map : Dict[str, Type[Analyzer]]
        Mapping of analyzer class names to their classes.
    configs : Dict[str, dict]
        Initialization parameters for each analyzer, keyed by class name.
    order : List[str]
        The ordered list of analyzer names to run.
    instances : List[Analyzer]
        Instantiated analyzer objects in current discovery order.

    Examples
    --------
    >>> from Analyzer import PluginManager, GCContentAnalyzer, MotifSearchAnalyzer
    >>> pm = PluginManager({'GCContentAnalyzer': {}, 'MotifSearchAnalyzer': {'motif': 'ATG'}})
    >>> pm.discover(['GCContentAnalyzer', 'MotifSearchAnalyzer'])
    GCContentAnalyzer discovered and added
    MotifSearchAnalyzer discovered and added
    >>> seq = Sequence.from_fasta(">s1\\nATGC\\n")
    >>> results = pm.run_instances(seq)
    >>> 'GCContentAnalyzer' in results
    True
    """     
    
    
    def __init__(self, analyzer_configs: Dict[Type[Analyzer], dict]):
        
        """
        Initialize the PluginManager with analyzer configurations.

        Parameters
        ----------
        analyzer_configs : dict
            A mapping from analyzer class name (str) to a dict of __init__ kwargs.
            Example:
            {
                'GCContentAnalyzer': {},
                'MotifSearchAnalyzer': {'motif': 'ATG'}
            }
        """
        # Build lookup of available Analyzer classes
        self.name_map = AnalyzerMeta.registry.copy()
        
        self.configs = analyzer_configs
        self.instances: List[Analyzer] = []
        self.order: List[str] = []
        
    def discover(self, action_order: List, verbose: bool = True):
        """
        Instantiate and register analyzers in the specified order.

        Parameters
        ----------
        action_order : list of str
            List of analyzer class names, in the order they should run.
        verbose : bool, optional
            If True, print “discovered and added” messages; by default True.

        Returns
        -------
        None

        Examples
        --------
        >>> pm.discover(['GCContentAnalyzer', 'MotifSearchAnalyzer'])
        GCContentAnalyzer discovered and added
        MotifSearchAnalyzer discovered and added
        """

        self.instances = []
        
        self.order = action_order.copy()
        

        for class_name in self.order:
            
            cls = self.name_map.get(class_name)
            
            if cls is None:
                if verbose:
                    print(f"{class_name} not present in order list.")
                    continue
            try:
                kwargs = self.configs[class_name]
            except KeyError:
                if verbose:
                    print(f"There are no configs for {class_name}")
                    continue
            # Instantiate and register
            self.instances.append(cls(**kwargs))
            if verbose:
                print(f"\n{class_name} discovered and added")

       
    def run_instances(self, seq: Sequence) -> Dict[str, Any]:
        """
        Run all registered analyzers on a single Sequence instance.

        Parameters
        ----------
        seq : Sequence
            The sequence object to analyze.

        Returns
        -------
        dict
            Mapping from analyzer class name to its result for this sequence.
        """
        results: Dict[str, Any] = {}
        for instance in self.instances:
            instance_name = type(instance).__name__
            try:
                results[instance_name] = instance.run(seq)
            except:
                results[instance_name] = 'Error'
                print(f"{instance_name} caused an error and won't be analyzed.")
                continue
        return results
    
    def run_all(self, sequences: Iterable[Sequence]) -> Dict[str, Dict[str, Any]]:
        """
        Run all registered analyzers on multiple sequences.

        Parameters
        ----------
        sequences : iterable of Sequence
            An iterable (list or generator) of Sequence instances.

        Returns
        -------
        dict
            Nested mapping: { sequence.header: { analyzer_name: result, ... }, ... }
        """
        multiple_sequence_results = {}
        
        for seq in sequences:
            multiple_sequence_results[seq.header] = self.run_instances(seq)
        
        return multiple_sequence_results
        
        
    def add_analyzer(self, cls, kwargs: Dict, order_position: int):
        """
        Dynamically register a new analyzer plugin.

        Parameters
        ----------
        cls : Type[Analyzer]
            The Analyzer subclass to add.
        kwargs : dict
            Initialization keyword arguments for the analyzer.
        order_position : int
            Zero-based index at which to insert the new analyzer in the order.

        Returns
        -------
        None

        Examples
        --------
        >>> pm.add_analyzer(MotifSearchAnalyzer, {'motif':'ATG'}, order_position=1)
        >>> pm.order
        ['GCContentAnalyzer', 'MotifSearchAnalyzer']
        """
        # Update lookup and configuration
        self.configs[cls.__name__] = kwargs
        # Insert into desired order, then rebuild instances silently
        self.order.insert(order_position, cls.__name__)
        print(f"\nAdded {cls.__name__}({kwargs});\nRebuilding plugin list…\n")
        self.discover(self.order, verbose=False)
        

                


pm = PluginManager({'GCContentAnalyzer':{}})
pm.discover(['GCContentAnalyzer'])
pm.add_analyzer(MotifSearchAnalyzer, {'motif':'ATG'}, order_position=1)
# Now pm.instances should already include both analyzers
print([type(i).__name__ for i in pm.instances])

pm.run_instances(Sequence('Seq1','GCGC', None))
