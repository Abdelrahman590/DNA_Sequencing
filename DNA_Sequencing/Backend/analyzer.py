from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight  # تحديث هنا
import numpy as np
import re
from collections import Counter
import io
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DNAAnalyzer:
    def __init__(self):
        self.common_motifs = ['TATA', 'CAAT', 'GCCGCC', 'ATGC', 'CGCG', 'ATAT']
    
    def clean_sequence(self, sequence):
        """تنظيف التسلسل من الرموز غير المرغوبة"""
        cleaned = re.sub(r'[^ATGC]', '', sequence.upper())
        return cleaned
    
    def calculate_gc_content(self, sequence):
        """حساب محتوى GC"""
        if len(sequence) == 0:
            return 0
        # استخدام gc_fraction بدلاً من GC
        return gc_fraction(Seq(sequence)) * 100
    
    def get_composition(self, sequence):
        """تحليل تركيب النيوكليوتيدات"""
        total = len(sequence)
        if total == 0:
            return {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        
        composition = {
            'A': (sequence.count('A') / total) * 100,
            'T': (sequence.count('T') / total) * 100, 
            'G': (sequence.count('G') / total) * 100,
            'C': (sequence.count('C') / total) * 100
        }
        return composition
    
    def find_motifs(self, sequence):
        """البحث عن المواضع الشائعة"""
        found_motifs = {}
        for motif in self.common_motifs:
            positions = []
            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1
            
            if positions:
                found_motifs[motif] = {
                    'count': len(positions),
                    'positions': positions[:10]
                }
        return found_motifs
    
    def translate_sequence(self, sequence):
        """ترجمة DNA إلى RNA وبروتين"""
        try:
            seq_obj = Seq(sequence)
            rna = seq_obj.transcribe()
            protein = seq_obj.translate()
            orfs = self.find_orfs(sequence)
            
            return {
                'rna': str(rna)[:500],
                'protein': str(protein)[:200],
                'orfs': orfs
            }
        except Exception as e:
            return {
                'rna': 'خطأ في الترجمة',
                'protein': 'خطأ في الترجمة',
                'orfs': [],
                'error': str(e)
            }
    
    def find_orfs(self, sequence):
        """إيجاد أطر القراءة المفتوحة (ORFs)"""
        orfs = []
        start_codon = 'ATG'
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        for frame in range(3):
            for i in range(frame, len(sequence) - 2, 3):
                codon = sequence[i:i+3]
                if len(codon) == 3 and codon == start_codon:
                    for j in range(i + 3, len(sequence) - 2, 3):
                        next_codon = sequence[j:j+3]
                        if len(next_codon) == 3 and next_codon in stop_codons:
                            orf_length = j - i + 3
                            if orf_length >= 90:
                                orfs.append({
                                    'start': i,
                                    'end': j + 3,
                                    'length': orf_length,
                                    'frame': frame + 1,
                                    'sequence': sequence[i:j+3][:60] + '...'
                                })
                            break
        
        return sorted(orfs, key=lambda x: x['length'], reverse=True)[:5]
    
    def analyze_sequence(self, sequence):
        """التحليل الشامل للتسلسل"""
        cleaned_seq = self.clean_sequence(sequence)
        
        if len(cleaned_seq) == 0:
            raise ValueError("التسلسل فارغ أو يحتوي على رموز غير صالحة")
        
        results = {
            'sequence_info': {
                'length': len(cleaned_seq),
                'original_length': len(sequence),
                'cleaned_sequence': cleaned_seq[:100] + ('...' if len(cleaned_seq) > 100 else '')
            },
            'gc_content': round(self.calculate_gc_content(cleaned_seq), 2),
            'composition': {k: round(v, 2) for k, v in self.get_composition(cleaned_seq).items()},
            'motifs': self.find_motifs(cleaned_seq),
            'translation': self.translate_sequence(cleaned_seq),
            'statistics': self.get_sequence_statistics(cleaned_seq)
        }
        
        return results
    
    def get_sequence_statistics(self, sequence):
        """إحصائيات إضافية"""
        seq_len = len(sequence)
        at_content = round((sequence.count('A') + sequence.count('T')) / seq_len * 100, 2)
        gc_content = round((sequence.count('G') + sequence.count('C')) / seq_len * 100, 2)
        purine_content = round((sequence.count('A') + sequence.count('G')) / seq_len * 100, 2)
        pyrimidine_content = round((sequence.count('C') + sequence.count('T')) / seq_len * 100, 2)
        
        # Generate sequence summary
        summary = []
        
        # Length analysis
        if seq_len < 100:
            summary.append("قصير جداً")
        elif seq_len < 1000:
            summary.append("متوسط الطول")
        else:
            summary.append("طويل")
            
        # GC content analysis
        if gc_content < 30:
            summary.append("محتوى GC منخفض")
        elif gc_content > 70:
            summary.append("محتوى GC مرتفع")
        else:
            summary.append("محتوى GC متوازن")
            
        # Balance analysis
        at_gc_diff = abs(at_content - gc_content)
        if at_gc_diff < 10:
            summary.append("توازن جيد بين AT و GC")
        
        return {
            'at_content': at_content,
            'gc_content': gc_content,
            'purine_content': purine_content,
            'pyrimidine_content': pyrimidine_content,
            'molecular_weight': round(molecular_weight(Seq(sequence), seq_type='DNA'), 2) if seq_len < 1000 else 'كبير جداً',
            'summary': ' | '.join(summary)
        }
    
    def detect_file_format(self, filepath):
        """Detect the format of the input file"""
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                second_line = f.readline().strip() if first_line else ""
                
                # Check for FASTA format (starts with >)
                if first_line.startswith('>'):
                    return 'fasta'
                
                # Check for FASTQ format (4-line repeating pattern starting with @)
                elif first_line.startswith('@'):
                    return 'fastq'
                
                # Check for TSV format (has tab and 'sequence' in header)
                elif '\t' in first_line and 'sequence' in first_line.lower():
                    return 'tsv'
                
                # Check if it's a raw sequence (contains only ATGCN)
                elif re.match(r'^[ATGCN\s]+$', first_line, re.I):
                    return 'raw'
                
            return 'unknown'
        except:
            return 'unknown'

    def analyze_from_file(self, filepath):
        """تحليل التسلسل من ملف"""
        sequences = []
        file_format = self.detect_file_format(filepath)
        
        try:
            # Handle based on detected format
            if file_format == 'fasta':
                # Parse FASTA format (including .fna files)
                for record in SeqIO.parse(filepath, "fasta"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq)
                    })
            
            elif file_format == 'fastq':
                # Parse FASTQ format
                for record in SeqIO.parse(filepath, "fastq"):
                    sequences.append({
                        'id': record.id,
                        'description': record.description,
                        'sequence': str(record.seq)
                    })
            
            elif file_format == 'tsv':
                # Parse TSV format
                with open(filepath, 'r') as f:
                    content = f.read()
                    lines = [line.strip() for line in content.split('\n') if line.strip()]
                    
                    # Parse header
                    header = lines[0].split('\t')
                    if 'sequence' in [col.lower().strip() for col in header]:
                        seq_index = [col.lower().strip() for col in header].index('sequence')
                        
                        # Process each line
                        for line in lines[1:]:
                            parts = line.split('\t')
                            if len(parts) > seq_index:
                                sequence = parts[seq_index].strip().replace(' ', '')
                                sequence_class = parts[1].strip() if len(parts) > 1 else "unknown"
                                
                                if sequence:  # Only add if we have a sequence
                                    sequences.append({
                                        'id': f'Sequence_{len(sequences)+1}',
                                        'description': f'Sequence from TSV file, class: {sequence_class}',
                                        'sequence': sequence
                                    })
            else:
                # Handle raw sequence format
                with open(filepath, 'r') as f:
                    content = f.read().strip()
                    # Join lines and clean the sequence
                    sequence = ''.join(line.strip() for line in content.split('\n'))
                    cleaned = self.clean_sequence(sequence)
                    if cleaned and len(cleaned) >= 10:  # Minimum sequence length
                        sequences.append({
                            'id': 'Sequence_1',
                            'description': 'Raw DNA sequence',
                            'sequence': cleaned
                        })
                    
        except Exception as e:
            print(f"File reading error: {str(e)}")
            logger.error(f"Error reading file {filepath}: {str(e)}")
        except:
            # Try as TSV/CSV
            try:
                with open(filepath, 'r') as f:
                    content = f.read()
                    lines = [line.strip() for line in content.split('\n') if line.strip()]
                    
                    # Parse header
                    header = lines[0].split('\t')
                    if 'sequence' in [col.lower().strip() for col in header]:
                        seq_index = [col.lower().strip() for col in header].index('sequence')
                        
                        # Join multi-line sequences
                        current_sequence = ""
                        current_class = ""
                        
                        for line in lines[1:]:
                            parts = line.split('\t')
                            if len(parts) > seq_index:
                                # If we have a complete sequence, add it
                                if current_sequence:
                                    sequences.append({
                                        'id': f'Sequence_{len(sequences)+1}',
                                        'description': f'Sequence from TSV file, class: {current_class}',
                                        'sequence': current_sequence.replace(' ', '')
                                    })
                                
                                # Start new sequence
                                current_sequence = parts[seq_index].strip()
                                current_class = parts[1].strip() if len(parts) > 1 else "unknown"
                            else:
                                # Continuation of previous sequence
                                current_sequence += line.strip()
                        
                        # Add the last sequence
                        if current_sequence:
                            sequences.append({
                                'id': f'Sequence_{len(sequences)+1}',
                                'description': f'Sequence from TSV file, class: {current_class}',
                                'sequence': current_sequence.replace(' ', '')
                            })
            except:
                # Try as plain text
                with open(filepath, 'r') as f:
                    content = f.read()
                    sequences.append({
                        'id': 'Unknown',
                        'description': 'من ملف نصي',
                        'sequence': content
                    })
        
        if not sequences:
            raise ValueError("لم يتم العثور على تسلسلات صالحة في الملف")
        
        first_sequence = sequences[0]
        results = self.analyze_sequence(first_sequence['sequence'])
        results['file_info'] = {
            'sequence_id': first_sequence['id'],
            'description': first_sequence['description'],
            'total_sequences': len(sequences)
        }
        
        return results
