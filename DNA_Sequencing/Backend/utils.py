from collections import Counter
from Bio.Seq import Seq

def get_nucleotide_composition(sequence):
    """Calculate the nucleotide composition of a DNA sequence."""
    if not sequence:
        return {}
    
    # Convert to uppercase to standardize counting
    sequence = sequence.upper()
    # Count nucleotides
    composition = Counter(sequence)
    # Calculate percentages
    total = len(sequence)
    return {base: (count / total) * 100 for base, count in composition.items()}

def find_common_motifs(sequence, motif_length=6, min_occurrences=2):
    """Find common motifs in a DNA sequence."""
    if not sequence or len(sequence) < motif_length:
        return []
    
    # Convert to uppercase to standardize searching
    sequence = sequence.upper()
    motifs = {}
    
    # Slide through the sequence and count motifs
    for i in range(len(sequence) - motif_length + 1):
        motif = sequence[i:i + motif_length]
        motifs[motif] = motifs.get(motif, 0) + 1
    
    # Filter motifs that appear more than once
    common_motifs = [(motif, count) for motif, count in motifs.items() 
                     if count >= min_occurrences]
    
    # Sort by frequency, highest first
    common_motifs.sort(key=lambda x: x[1], reverse=True)
    
    return common_motifs
