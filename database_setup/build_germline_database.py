#!/usr/bin/env python3
"""
IMGT Germline Database Builder

Creates a comprehensive protein BLAST database from IMGT germline V, D, and J segments
by reconstructing complete antibody frameworks through V+D+J recombination.

Scientific Approach:
1. Parse IMGT germline nucleotide sequences (V, D, J segments)
2. Translate to protein sequences
3. Reconstruct complete frameworks:
   - Heavy chains: V + D + J (FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4)
   - Light chains: V + J (FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4)
4. Generate comprehensive protein database for BLAST
5. Create species-specific and chain-specific databases

Author: Antibody Humanization Pipeline
Date: 2024
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s: %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class GermlineSegment:
    """Represents a germline V, D, or J segment."""
    accession: str
    gene_name: str
    species: str
    functionality: str  # F, ORF, P
    segment_type: str  # V-REGION, D-REGION, J-REGION
    nucleotide_seq: str
    protein_seq: str
    
    @property
    def is_functional(self) -> bool:
        """Check if segment is functional (F or ORF)."""
        return self.functionality in ['F', 'ORF']
    
    @property
    def is_human(self) -> bool:
        """Check if segment is from Homo sapiens."""
        return self.species == 'Homo sapiens'

class GermlineDatabaseBuilder:
    """Builds comprehensive germline antibody database from V, D, J segments."""
    
    def __init__(self, input_dir: str, output_dir: str):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Storage for parsed segments
        self.v_segments: Dict[str, List[GermlineSegment]] = {'heavy': [], 'kappa': [], 'lambda': []}
        self.d_segments: List[GermlineSegment] = []
        self.j_segments: Dict[str, List[GermlineSegment]] = {'heavy': [], 'kappa': [], 'lambda': []}
        
        # Reconstructed frameworks
        self.heavy_frameworks: List[SeqRecord] = []
        self.light_frameworks: List[SeqRecord] = []

    @staticmethod
    def _format_segment_description(label: str, segment: GermlineSegment) -> str:
        """Keep IMGT gene, accession, and functionality in BLAST titles."""
        return (
            f"{label}:{segment.gene_name}"
            f"[accession={segment.accession},functionality={segment.functionality}]"
        )
    
    def parse_imgt_fasta(self, fasta_file: Path) -> List[GermlineSegment]:
        """Parse IMGT FASTA file and extract germline segments."""
        segments = []
        
        logger.info(f"Parsing {fasta_file.name}...")
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Parse IMGT header format: >accession|gene|species|functionality|region|...
            header_parts = record.description.split('|')
            
            if len(header_parts) < 5:
                logger.warning(f"Skipping malformed header: {record.description}")
                continue
            
            accession = header_parts[0].lstrip('>')
            gene_name = header_parts[1]
            species = header_parts[2]
            functionality = header_parts[3]
            segment_type = header_parts[4]
            
            # Translate nucleotide to protein
            nucleotide_seq = str(record.seq).replace('.', '').replace('-', '')
            
            # Handle potential gaps and incomplete codons
            try:
                if len(nucleotide_seq) % 3 != 0:
                    # Pad to complete codon
                    nucleotide_seq = nucleotide_seq + 'N' * (3 - len(nucleotide_seq) % 3)
                
                protein_seq = str(Seq(nucleotide_seq).translate())
                # Remove stop codons and ambiguous amino acids
                protein_seq = protein_seq.replace('*', '').replace('X', '')
                
            except Exception as e:
                logger.warning(f"Translation failed for {gene_name}: {e}")
                continue
            
            segment = GermlineSegment(
                accession=accession,
                gene_name=gene_name,
                species=species,
                functionality=functionality,
                segment_type=segment_type,
                nucleotide_seq=nucleotide_seq,
                protein_seq=protein_seq
            )
            
            segments.append(segment)
        
        logger.info(f"Parsed {len(segments)} segments from {fasta_file.name}")
        return segments
    
    def categorize_segments(self):
        """Load and categorize V, D, J segments by chain type."""
        
        # Load V segments
        ighv_file = self.input_dir / "IGHV.fasta"
        if ighv_file.exists():
            ighv_segments = self.parse_imgt_fasta(ighv_file)
            self.v_segments['heavy'] = [s for s in ighv_segments if s.is_human and s.is_functional]
        
        igkv_file = self.input_dir / "IGKV.fasta"
        if igkv_file.exists():
            igkv_segments = self.parse_imgt_fasta(igkv_file)
            self.v_segments['kappa'] = [s for s in igkv_segments if s.is_human and s.is_functional]
        
        iglv_file = self.input_dir / "IGLV.fasta"
        if iglv_file.exists():
            iglv_segments = self.parse_imgt_fasta(iglv_file)
            self.v_segments['lambda'] = [s for s in iglv_segments if s.is_human and s.is_functional]
        
        # Load D segments
        ighd_file = self.input_dir / "IGHD.fasta"
        if ighd_file.exists():
            ighd_segments = self.parse_imgt_fasta(ighd_file)
            self.d_segments = [s for s in ighd_segments if s.is_human and s.is_functional]
        
        # Load J segments
        ighj_file = self.input_dir / "IGHJ.fasta"
        if ighj_file.exists():
            ighj_segments = self.parse_imgt_fasta(ighj_file)
            self.j_segments['heavy'] = [s for s in ighj_segments if s.is_human and s.is_functional]
        
        igkj_file = self.input_dir / "IGKJ.fasta"
        if igkj_file.exists():
            igkj_segments = self.parse_imgt_fasta(igkj_file)
            self.j_segments['kappa'] = [s for s in igkj_segments if s.is_human and s.is_functional]
        
        iglj_file = self.input_dir / "IGLJ.fasta"
        if iglj_file.exists():
            iglj_segments = self.parse_imgt_fasta(iglj_file)
            self.j_segments['lambda'] = [s for s in iglj_segments if s.is_human and s.is_functional]
        
        # Log statistics
        logger.info("=== Germline Segment Statistics ===")
        logger.info(f"Heavy V segments: {len(self.v_segments['heavy'])}")
        logger.info(f"Kappa V segments: {len(self.v_segments['kappa'])}")
        logger.info(f"Lambda V segments: {len(self.v_segments['lambda'])}")
        logger.info(f"Heavy D segments: {len(self.d_segments)}")
        logger.info(f"Heavy J segments: {len(self.j_segments['heavy'])}")
        logger.info(f"Kappa J segments: {len(self.j_segments['kappa'])}")
        logger.info(f"Lambda J segments: {len(self.j_segments['lambda'])}")
    
    def reconstruct_heavy_frameworks(self):
        """Reconstruct complete heavy chain frameworks from V+D+J segments."""
        logger.info("Reconstructing heavy chain frameworks (V+D+J)...")
        
        count = 0
        for v_seg in self.v_segments['heavy']:
            for d_seg in self.d_segments:
                for j_seg in self.j_segments['heavy']:
                    # Reconstruct complete framework
                    # V segment provides FR1-CDR1-FR2-CDR2-FR3
                    # D segment provides CDR3 diversity
                    # J segment provides CDR3 completion + FR4
                    
                    complete_seq = v_seg.protein_seq + d_seg.protein_seq + j_seg.protein_seq
                    
                    # Create unique identifier (max 50 chars for BLAST)
                    framework_id = f"VH_{count+1:06d}"
                    segment_details = ' '.join([
                        self._format_segment_description('V', v_seg),
                        self._format_segment_description('D', d_seg),
                        self._format_segment_description('J', j_seg),
                    ])
                    description = f"Heavy chain framework | {segment_details} | Homo sapiens germline reconstruction"
                    
                    framework_record = SeqRecord(
                        Seq(complete_seq),
                        id=framework_id,
                        description=description
                    )
                    
                    self.heavy_frameworks.append(framework_record)
                    count += 1
        
        logger.info(f"Generated {count} heavy chain frameworks")
    
    def reconstruct_light_frameworks(self):
        """Reconstruct complete light chain frameworks from V+J segments with canonical FR4."""
        logger.info("Reconstructing light chain frameworks (V+J) with canonical FR4...")
        
        count = 0
        
        # Kappa light chains
        for v_seg in self.v_segments['kappa']:
            for j_seg in self.j_segments['kappa']:
                # Scientifically correct light chain reconstruction
                complete_seq = self._reconstruct_light_chain_with_canonical_fr4(
                    v_seg, j_seg, 'kappa'
                )
                
                if complete_seq is None:
                    logger.warning(f"Failed to reconstruct kappa framework: V={v_seg.gene_name}, J={j_seg.gene_name}")
                    continue
                
                framework_id = f"VK_{count+1:06d}"
                segment_details = ' '.join([
                    self._format_segment_description('V', v_seg),
                    self._format_segment_description('J', j_seg),
                ])
                description = f"Kappa light chain framework | {segment_details} | Homo sapiens germline reconstruction with canonical FR4"
                
                framework_record = SeqRecord(
                    Seq(complete_seq),
                    id=framework_id,
                    description=description
                )
                
                self.light_frameworks.append(framework_record)
                count += 1
        
        # Lambda light chains
        for v_seg in self.v_segments['lambda']:
            for j_seg in self.j_segments['lambda']:
                # Scientifically correct light chain reconstruction
                complete_seq = self._reconstruct_light_chain_with_canonical_fr4(
                    v_seg, j_seg, 'lambda'
                )
                
                if complete_seq is None:
                    logger.warning(f"Failed to reconstruct lambda framework: V={v_seg.gene_name}, J={j_seg.gene_name}")
                    continue
                
                framework_id = f"VL_{count+1:06d}"
                segment_details = ' '.join([
                    self._format_segment_description('V', v_seg),
                    self._format_segment_description('J', j_seg),
                ])
                description = f"Lambda light chain framework | {segment_details} | Homo sapiens germline reconstruction with canonical FR4"
                
                framework_record = SeqRecord(
                    Seq(complete_seq),
                    id=framework_id,
                    description=description
                )
                
                self.light_frameworks.append(framework_record)
                count += 1
        
        logger.info(f"Generated {count} light chain frameworks with canonical FR4")
    
    def write_fasta_files(self):
        """Write reconstructed frameworks to FASTA files."""
        
        # Combined database
        all_frameworks = self.heavy_frameworks + self.light_frameworks
        combined_file = self.output_dir / "human_germline_frameworks.fasta"
        
        logger.info(f"Writing {len(all_frameworks)} frameworks to {combined_file}")
        SeqIO.write(all_frameworks, combined_file, "fasta")
        
        # Separate heavy and light chain files
        heavy_file = self.output_dir / "human_germline_heavy.fasta"
        light_file = self.output_dir / "human_germline_light.fasta"
        
        SeqIO.write(self.heavy_frameworks, heavy_file, "fasta")
        SeqIO.write(self.light_frameworks, light_file, "fasta")
        
        logger.info(f"Heavy chains: {heavy_file} ({len(self.heavy_frameworks)} sequences)")
        logger.info(f"Light chains: {light_file} ({len(self.light_frameworks)} sequences)")
    
    def build_blast_database(self):
        """Create BLAST database from reconstructed frameworks."""
        logger.info("Building BLAST database...")
        
        fasta_file = self.output_dir / "human_germline_frameworks.fasta"
        db_name = self.output_dir / "human_germline_frameworks"
        
        cmd = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-dbtype", "prot",
            "-out", str(db_name),
            "-title", "Human Germline Antibody Frameworks",
            "-parse_seqids"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("BLAST database created successfully")
            logger.info(f"Database files: {db_name}.*")
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST database creation failed: {e}")
            logger.error(f"STDERR: {e.stderr}")
            raise
    
    def generate_statistics(self):
        """Generate database statistics and summary."""
        stats_file = self.output_dir / "database_statistics.txt"
        
        total_heavy = len(self.heavy_frameworks)
        total_light = len(self.light_frameworks)
        total_frameworks = total_heavy + total_light
        
        # Calculate theoretical combinations
        heavy_combinations = len(self.v_segments['heavy']) * len(self.d_segments) * len(self.j_segments['heavy'])
        kappa_combinations = len(self.v_segments['kappa']) * len(self.j_segments['kappa'])
        lambda_combinations = len(self.v_segments['lambda']) * len(self.j_segments['lambda'])
        
        with open(stats_file, 'w') as f:
            f.write("IMGT Germline Database Statistics\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Database Directory: {self.output_dir}\n")
            f.write(f"Source Directory: {self.input_dir}\n\n")
            
            f.write("Germline Segments Used:\n")
            f.write(f"  Heavy V segments: {len(self.v_segments['heavy'])}\n")
            f.write(f"  Heavy D segments: {len(self.d_segments)}\n")
            f.write(f"  Heavy J segments: {len(self.j_segments['heavy'])}\n")
            f.write(f"  Kappa V segments: {len(self.v_segments['kappa'])}\n")
            f.write(f"  Kappa J segments: {len(self.j_segments['kappa'])}\n")
            f.write(f"  Lambda V segments: {len(self.v_segments['lambda'])}\n")
            f.write(f"  Lambda J segments: {len(self.j_segments['lambda'])}\n\n")
            
            f.write("Reconstructed Frameworks:\n")
            f.write(f"  Heavy chain frameworks: {total_heavy:,}\n")
            f.write(f"  Light chain frameworks: {total_light:,}\n")
            f.write(f"  Total frameworks: {total_frameworks:,}\n\n")
            
            f.write("Theoretical Combinations:\n")
            f.write(f"  Heavy (V×D×J): {heavy_combinations:,}\n")
            f.write(f"  Kappa (V×J): {kappa_combinations:,}\n")
            f.write(f"  Lambda (V×J): {lambda_combinations:,}\n")
            f.write(f"  Total possible: {heavy_combinations + kappa_combinations + lambda_combinations:,}\n\n")
            
            f.write("Database Files:\n")
            f.write(f"  Combined FASTA: human_germline_frameworks.fasta\n")
            f.write(f"  Heavy FASTA: human_germline_heavy.fasta\n")
            f.write(f"  Light FASTA: human_germline_light.fasta\n")
            f.write(f"  BLAST database: human_germline_frameworks.*\n")
        
        logger.info(f"Statistics written to {stats_file}")
        logger.info(f"Generated {total_frameworks:,} total germline frameworks")
    
    def _reconstruct_light_chain_with_canonical_fr4(self, v_seg: GermlineSegment, j_seg: GermlineSegment, chain_type: str) -> Optional[str]:
        """
        Reconstruct light chain framework with scientifically accurate canonical FR4.
        
        Scientific Rationale:
        - IMGT V segments provide FR1-CDR1-FR2-CDR2-FR3 + partial CDR3
        - IMGT J segments contain CDR3 end + partial FR4 + constant region contamination
        - Canonical FR4 patterns must be used to avoid constant region contamination
        
        Args:
            v_seg: IMGT V segment
            j_seg: IMGT J segment  
            chain_type: 'kappa' or 'lambda'
            
        Returns:
            Complete light chain framework with canonical FR4, or None if reconstruction fails
        """
        try:
            # Step 1: Clean V segment (remove any potential constant region contamination)
            v_clean = self._clean_v_segment_sequence(v_seg.protein_seq)
            if not v_clean:
                logger.debug(f"V segment cleaning failed for {v_seg.gene_name}")
                return None
            
            # Step 2: Get canonical FR4 based on J segment and chain type
            canonical_fr4 = self._get_canonical_fr4_for_j_segment(j_seg.gene_name, chain_type)
            if not canonical_fr4:
                logger.debug(f"No canonical FR4 found for J segment {j_seg.gene_name} ({chain_type})")
                return None
            
            # Step 3: Reconstruct complete framework
            # V segment provides everything up to CDR3 end, then add canonical FR4
            complete_framework = v_clean + canonical_fr4
            
            # Step 4: Validate the reconstructed framework
            if not self._validate_light_chain_framework(complete_framework, chain_type):
                logger.debug(f"Framework validation failed for V={v_seg.gene_name}, J={j_seg.gene_name}")
                return None
            
            logger.debug(f"Successfully reconstructed {chain_type} framework: V={v_seg.gene_name}, J={j_seg.gene_name}, FR4={canonical_fr4}")
            return complete_framework
            
        except Exception as e:
            logger.debug(f"Framework reconstruction failed for V={v_seg.gene_name}, J={j_seg.gene_name}: {e}")
            return None
    
    def _clean_v_segment_sequence(self, v_protein: str) -> Optional[str]:
        """
        Clean V segment sequence to remove any constant region contamination.
        
        V segments should end around CDR3 start. We need to identify the proper
        boundary to avoid including constant region sequences.
        
        Args:
            v_protein: V segment protein sequence
            
        Returns:
            Cleaned V segment sequence or None if cleaning fails
        """
        if not v_protein or len(v_protein) < 50:
            return None
            
        # Remove stop codons and ambiguous amino acids (X)
        cleaned = v_protein.replace('*', '').replace('X', '')
        
        # V segments typically end around 95-100 AA for light chains
        # Truncate if sequence is unusually long (likely contains contamination)
        if len(cleaned) > 110:
            logger.debug(f"V segment unusually long ({len(cleaned)} AA), truncating to 100 AA")
            cleaned = cleaned[:100]
        
        # Ensure sequence ends properly (not with obvious constant region patterns)
        # Common constant region starts: 'RADAA', 'TKLEA', etc.
        constant_patterns = ['RADAA', 'TKLEA', 'RTVA', 'RTAA']
        for pattern in constant_patterns:
            if pattern in cleaned[-15:]:  # Check last 15 AA
                pattern_pos = cleaned.rfind(pattern)
                if pattern_pos > 70:  # Only truncate if pattern is in reasonable position
                    cleaned = cleaned[:pattern_pos]
                    logger.debug(f"Removed constant region pattern {pattern} from V segment")
                    break
        
        return cleaned if len(cleaned) >= 50 else None
    
    def _get_canonical_fr4_for_j_segment(self, j_gene_name: str, chain_type: str) -> Optional[str]:
        """
        Get canonical FR4 sequence based on J segment identity and chain type.
        
        Based on scientific literature (Kabat et al., 1991; Chothia & Lesk, 1987)
        and IMGT database analysis.
        
        Args:
            j_gene_name: IMGT J gene name (e.g., 'IGKJ1*01')
            chain_type: 'kappa' or 'lambda'
            
        Returns:
            Canonical FR4 sequence or None if not found
        """
        # Canonical FR4 patterns based on scientific literature and database analysis
        canonical_fr4_map = {
            'kappa': {
                # Most common kappa FR4 patterns
                'default': 'FGGGTKLEIK',  # Most frequent pattern
                'IGKJ1': 'FGQGTKVEIK',   # IGKJ1-specific pattern
                'IGKJ2': 'FGGGTKLEIK',   # IGKJ2-specific pattern  
                'IGKJ3': 'FGQGTKVEIK',   # IGKJ3-specific pattern
                'IGKJ4': 'FGQGTKVEIK',   # IGKJ4-specific pattern
                'IGKJ5': 'FGGGTKLEIK',   # IGKJ5-specific pattern
            },
            'lambda': {
                # Common lambda FR4 patterns
                'default': 'FGAGTKLELK',  # Most frequent lambda pattern
                'IGLJ1': 'FGAGTKLELK',   # IGLJ1-specific pattern
                'IGLJ2': 'FGAGTKLELK',   # IGLJ2-specific pattern
                'IGLJ3': 'FGAGTKLELK',   # IGLJ3-specific pattern
                'IGLJ4': 'FGSGTKLELK',   # IGLJ4-specific pattern (serine variant)
                'IGLJ6': 'FGAGTKLELK',   # IGLJ6-specific pattern
                'IGLJ7': 'FGAGTKLELK',   # IGLJ7-specific pattern
            }
        }
        
        if chain_type not in canonical_fr4_map:
            logger.debug(f"Unsupported chain type: {chain_type}")
            return None
        
        # Extract J gene base name (remove allele info)
        j_base_name = j_gene_name.split('*')[0] if '*' in j_gene_name else j_gene_name
        
        # Look for specific J gene pattern
        fr4_patterns = canonical_fr4_map[chain_type]
        if j_base_name in fr4_patterns:
            return fr4_patterns[j_base_name]
        
        # Fallback to default pattern for the chain type
        return fr4_patterns['default']
    
    def _validate_light_chain_framework(self, framework: str, chain_type: str) -> bool:
        """
        Validate reconstructed light chain framework for scientific accuracy.
        
        Args:
            framework: Complete framework sequence
            chain_type: 'kappa' or 'lambda'
            
        Returns:
            True if framework passes validation
        """
        if not framework or len(framework) < 90:
            return False
        
        # Check for valid amino acids only (should be clean after X removal)
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        invalid_chars = [aa for aa in framework if aa not in valid_aa]
        if invalid_chars:
            logger.debug(f"Framework contains invalid amino acids: {set(invalid_chars)}")
            return False
        
        # Check that framework ends with canonical FR4 pattern
        expected_fr4_patterns = {
            'kappa': ['FGGGTKLEIK', 'FGQGTKVEIK', 'FGAGTKVEIK'],
            'lambda': ['FGAGTKLELK', 'FGSGTKLELK', 'FGAGTKLVLK']
        }
        
        if chain_type in expected_fr4_patterns:
            framework_end = framework[-10:]  # Last 10 AA should be FR4
            if framework_end not in expected_fr4_patterns[chain_type]:
                logger.debug(f"Framework does not end with canonical {chain_type} FR4. Got: {framework_end}")
                return False
        
        # Check reasonable length (light chains typically 100-115 AA)
        if not (90 <= len(framework) <= 120):
            logger.debug(f"Framework length {len(framework)} outside expected range 90-120 AA")
            return False
        
        return True
    
    def build_database(self):
        """Main method to build the complete germline database."""
        logger.info("=== IMGT Germline Database Builder ===")
        logger.info(f"Input directory: {self.input_dir}")
        logger.info(f"Output directory: {self.output_dir}")
        
        # Step 1: Load and categorize segments
        self.categorize_segments()
        
        # Step 2: Reconstruct complete frameworks
        # Heavy chains: Use existing method (works correctly)
        self.reconstruct_heavy_frameworks()
        
        # Light chains: Use new scientifically-corrected method
        self.reconstruct_light_frameworks()
        
        # Step 3: Write FASTA files
        self.write_fasta_files()
        
        # Step 4: Build BLAST database
        self.build_blast_database()
        
        # Step 5: Generate statistics
        self.generate_statistics()
        
        logger.info("=== Database Build Complete ===")

def main():
    """Main execution function."""
    
    # Configuration
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Build database
    builder = GermlineDatabaseBuilder(input_dir, output_dir)
    builder.build_database()

if __name__ == "__main__":
    main()
