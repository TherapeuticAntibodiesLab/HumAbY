#!/usr/bin/env python3
"""
Antibody Humanizer - VDJ Germline Database Approach

Uses clean VDJ-reconstructed germline sequences for CDR grafting.
"""

import argparse
import copy
import logging
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, TextIO
from Bio import SeqIO
from cdr import ANARCIICDRExtractor, get_consistency_manager, extract_cdrs_consistent

# Configuration Constants
MAX_DISPLAY_CANDIDATES = 10  # Balance between informativeness and readability
MAX_GENERATED_CANDIDATES = 5  # Limit to prevent excessive candidate generation
BLAST_MAX_TARGETS = 50  # Maximum BLAST search targets for performance
CHAIN_LENGTH_THRESHOLD = 110  # Heavy chains typically >110 AA, light chains <110 AA
APP_DIR = Path(__file__).resolve().parent
BACK_DIR = APP_DIR / "back"
BACK_SCRIPT = BACK_DIR / "back.py"
DEFAULT_MOUSE_GERMLINE_DB = Path(
    os.getenv(
        "MOUSE_GERMLINE_DB",
        str(BACK_DIR / "ncbi-igblast" / "internal_data" / "mouse" / "mouse_V"),
    )
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s: %(message)s')
logger = logging.getLogger(__name__)

GERMLINE_FIELD_PATTERN = re.compile(
    r"\b(V|D|J):([^\s|\[]+)(?:\[([^\]]+)\])?"
)

# =============================================================================
# CUSTOM EXCEPTIONS - Robust Error Handling
# =============================================================================

class HumanizerError(Exception):
    """Base exception for antibody humanizer errors."""
    pass

class SequenceValidationError(HumanizerError):
    """Raised when sequence validation fails."""
    pass

class DatabaseError(HumanizerError):
    """Raised when database operations fail."""
    pass

class OptimizationError(HumanizerError):
    """Raised when optimization processes fail."""
    pass

class CDRExtractionError(HumanizerError):
    """Raised when CDR extraction fails."""
    pass

# =============================================================================
# UTILITY FUNCTIONS - Software Engineering Best Practices
# =============================================================================

def display_candidates(candidates: List[Dict], chain_type: str, max_display: int = MAX_DISPLAY_CANDIDATES) -> None:
    """
    Display candidate information in standardized format.
    
    Args:
        candidates: List of candidate dictionaries with seq_id, identity, evalue, etc.
        chain_type: Chain type identifier (e.g., 'VH', 'VL')
        max_display: Maximum number of candidates to display
        
    Scientific rationale: Standardized display prevents information overload while
    maintaining essential metrics for candidate evaluation.
    """
    if candidates:
        logger.info(f"Found {len(candidates)} {chain_type} candidates:")
        for i, candidate in enumerate(candidates[:max_display], 1):
            logger.info(f"  {i:2d}. ID: {candidate['seq_id']}")
            logger.info(f"      Identity: {candidate['identity']:.1f}%")
            logger.info(f"      E-value: {candidate['evalue']:.2e}")
            logger.info(f"      BitScore: {candidate['bitscore']:.1f}")
            logger.info(f"      Length: {len(candidate['sequence'])} AA")
            logger.info(f"      IMGT germline: {format_germline_details(candidate)}")
            logger.info(f"      Sequence: {candidate['sequence']}")
            logger.info("")
    else:
        logger.error(f"❌ No {chain_type} candidates found")

def mark_candidates_as_unoptimized(vh_candidates: List[Dict], vl_candidates: List[Dict], 
                                 level: int) -> None:
    """
    Mark all candidates as not optimized with consistent flags.
    
    Args:
        vh_candidates: Heavy chain candidate list
        vl_candidates: Light chain candidate list  
        level: Optimization level attempted
        
    Scientific rationale: Consistent metadata tracking enables proper result
    interpretation and downstream analysis validation.
    """
    optimization_metadata = {
        'optimization_applied': False,
        'optimization_score': 0.0,
        'optimization_level': level
    }
    
    for candidate in vh_candidates:
        candidate.update(optimization_metadata)
    for candidate in vl_candidates:
        candidate.update(optimization_metadata)

def log_success(message: str) -> None:
    """Log success message with consistent formatting."""
    logger.info(f"✅ {message}")

def log_warning(message: str) -> None:
    """Log warning message with consistent formatting."""
    logger.warning(f"⚠️  {message}")

def log_error(message: str) -> None:
    """Log error message with consistent formatting."""
    logger.error(f"❌ {message}")

def log_info_header(message: str, separator_char: str = "=", width: int = 50) -> None:
    """Log formatted header for section separation."""
    logger.info(f"🔍 {message}")
    logger.info(separator_char * width)

def validate_sequence(sequence: str, min_length: int = 50) -> str:
    """
    Validate and clean protein sequence.
    
    Args:
        sequence: Raw protein sequence
        min_length: Minimum required sequence length
        
    Returns:
        Cleaned and validated sequence
        
    Raises:
        SequenceValidationError: If sequence is invalid
        
    Scientific rationale: Ensures sequence quality before processing to prevent
    downstream errors and maintain data integrity.
    """
    if not sequence or not isinstance(sequence, str):
        raise SequenceValidationError("Sequence must be a non-empty string")
    
    # Clean sequence
    cleaned_seq = sequence.strip().upper()
    
    # Remove any non-amino acid characters
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    cleaned_seq = ''.join(c for c in cleaned_seq if c in valid_aa)
    
    if len(cleaned_seq) < min_length:
        raise SequenceValidationError(f"Sequence too short: {len(cleaned_seq)} < {min_length}")
    
    # Check for unusual patterns that might indicate data corruption
    if len(set(cleaned_seq)) < 5:  # Too few unique amino acids
        raise SequenceValidationError("Sequence has insufficient amino acid diversity")
    
    return cleaned_seq

def safe_execute_subprocess(cmd: List[str], timeout: int = 30) -> subprocess.CompletedProcess:
    """
    Safely execute subprocess with consistent error handling.
    
    Args:
        cmd: Command and arguments list
        timeout: Timeout in seconds
        
    Returns:
        CompletedProcess result
        
    Raises:
        DatabaseError: If subprocess execution fails
    """
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, 
                              check=True, timeout=timeout)
        return result
    except subprocess.TimeoutExpired:
        raise DatabaseError(f"Command timed out after {timeout}s: {' '.join(cmd[:2])}")
    except subprocess.CalledProcessError as e:
        raise DatabaseError(f"Command failed: {' '.join(cmd[:2])}\nError: {e.stderr}")
    except FileNotFoundError:
        raise DatabaseError(f"Command not found: {cmd[0]}")

def write_candidate_fasta(candidates: List[Dict], output_path: Path, chain_type: str) -> None:
    """
    Write candidates to FASTA format file.
    
    Args:
        candidates: List of candidate dictionaries
        output_path: Output directory path
        chain_type: Chain type ('VH' or 'VL')
        
    Scientific rationale: Standardized FASTA output enables downstream analysis
    and integration with other bioinformatics tools.
    """
    if candidates:
        fasta_file = output_path / f"humanized_{chain_type.lower()}.fasta"
        with open(fasta_file, 'w') as f:
            for i, candidate in enumerate(candidates, 1):
                germline = candidate.get('germline', {})
                gene_fields = [
                    f"{segment}_{germline[segment]}"
                    for segment in ('V', 'D', 'J')
                    if germline.get(segment)
                ]
                header_parts = [
                    f">{chain_type}_candidate_{i}",
                    candidate['source'],
                    f"identity_{candidate['identity']:.1f}%",
                    *gene_fields,
                ]
                header = '|'.join(header_parts)
                f.write(f"{header}\n{candidate['sequence']}\n")

def write_candidate_summary(candidates: List[Dict], file_handle, chain_type: str) -> None:
    """
    Write candidate summary to file handle.
    
    Args:
        candidates: List of candidate dictionaries
        file_handle: Open file handle for writing
        chain_type: Chain type ('VH' or 'VL')
        
    Scientific rationale: Structured summary format enables rapid candidate
    evaluation and selection for downstream applications.
    """
    file_handle.write(f"\n{chain_type} Chain: {len(candidates)} candidates\n")
    for i, candidate in enumerate(candidates, 1):
        file_handle.write(
            f"  {i}. {candidate['source']} "
            f"(identity: {candidate['identity']:.1f}%, "
            f"length: {candidate['length']} AA)\n"
        )
        file_handle.write(f"     IMGT germline: {format_germline_details(candidate)}\n")


def parse_germline_metadata(seq_id: str, title: str) -> Dict[str, str]:
    """Extract IMGT germline nomenclature from a reconstructed database title."""
    metadata = {
        'framework_id': seq_id,
        'description': title.strip(),
        'species': 'Homo sapiens',
    }

    clean_id = seq_id.strip().strip('|').split('|')[-1]
    if clean_id.startswith('VH_'):
        metadata['chain'] = 'heavy'
    elif clean_id.startswith('VK_'):
        metadata['chain'] = 'kappa'
    elif clean_id.startswith('VL_'):
        metadata['chain'] = 'lambda'

    for segment, gene_name, details in GERMLINE_FIELD_PATTERN.findall(title):
        metadata[segment] = gene_name
        if details:
            detail_parts = [part.strip() for part in details.split(',')]
            for part in detail_parts:
                if part.startswith('accession='):
                    metadata[f'{segment}_accession'] = part.split('=', 1)[1]
                elif part.startswith('functionality='):
                    metadata[f'{segment}_functionality'] = part.split('=', 1)[1]

    return metadata


def format_germline_details(candidate: Dict) -> str:
    """Format germline metadata for logs and human-readable result files."""
    germline = candidate.get('germline', {})
    chain_labels = {
        'heavy': 'heavy (IGH)',
        'kappa': 'kappa (IGK)',
        'lambda': 'lambda (IGL)',
    }
    details = []

    if germline.get('chain'):
        details.append(f"chain={chain_labels.get(germline['chain'], germline['chain'])}")

    for segment in ('V', 'D', 'J'):
        gene_name = germline.get(segment)
        if not gene_name:
            continue
        annotations = []
        if germline.get(f'{segment}_accession'):
            annotations.append(f"accession={germline[f'{segment}_accession']}")
        if germline.get(f'{segment}_functionality'):
            annotations.append(f"functionality={germline[f'{segment}_functionality']}")
        suffix = f" ({', '.join(annotations)})" if annotations else ''
        details.append(f"{segment}={gene_name}{suffix}")

    if germline.get('species'):
        details.append(f"species={germline['species']}")

    return '; '.join(details) if details else 'not available in database description'


def prepare_backmutation_artifacts(output_dir: str) -> Path:
    """Prepare output files consumed by the original back/back.py script."""
    artifact_dir = Path(output_dir) / "backmutation" / "opt1"
    artifact_dir.mkdir(parents=True, exist_ok=True)

    for filename in ("FRfileVH", "FRfileVL", "hGermVH", "hGermVL", "humanized_vh.fasta", "humanized_vl.fasta"):
        path = artifact_dir / filename
        if path.exists():
            path.unlink()

    return artifact_dir


def write_backmutation_candidates(results: Dict[str, List], artifact_dir: Optional[Path]) -> None:
    """Write the graft-stage candidates used as the original back.py opt1 input."""
    if artifact_dir is None:
        return

    for chain_type, candidates in results.items():
        write_candidate_fasta(candidates, artifact_dir, chain_type)


def write_backmutation_graft_artifacts(
    murine_regions: Dict[str, str],
    human_regions: Dict[str, str],
    chain_type: str,
    artifact_dir: Optional[Path],
) -> None:
    """Write FRfile* and hGerm* artifacts in the format expected by back.py."""
    if artifact_dir is None:
        return

    artifact_dir.mkdir(parents=True, exist_ok=True)
    fr_file = artifact_dir / f"FRfile{chain_type}"
    with open(fr_file, "w") as handle:
        for region in ("FR1", "FR2", "FR3", "FR4"):
            print(murine_regions.get(region, ""), file=handle)

    human_germline = (
        human_regions.get("FR1", "") +
        human_regions.get("CDR1", "") +
        human_regions.get("FR2", "") +
        human_regions.get("CDR2", "") +
        human_regions.get("FR3", "") +
        human_regions.get("CDR3", "") +
        human_regions.get("FR4", "")
    )
    with open(artifact_dir / f"hGerm{chain_type}", "a") as handle:
        print(human_germline, file=handle)


def _write_backmutation_log(log_handle: Optional[TextIO], message: str) -> None:
    """Write a timestamped message to the dedicated backmutation log."""
    if log_handle is None:
        return
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", file=log_handle)
    log_handle.flush()


def _format_backmutated_germline_info(
    sequence: str,
    chain_type: str,
    database_path: str,
    log_handle: Optional[TextIO] = None,
) -> List[str]:
    """Determine germline annotations from the backmutated sequence itself."""
    _write_backmutation_log(
        log_handle,
        f"Analyzing {chain_type} backmutated sequence against human germline database; length={len(sequence)}",
    )
    _write_backmutation_log(
        log_handle,
        (
            "Germline decision method: run blastp against the human germline framework database, "
            "retrieve full hit titles with blastdbcmd, parse IMGT V/D/J nomenclature from those titles, "
            "sort hits by BLAST bitscore descending, and select rank #1."
        ),
    )
    _write_backmutation_log(
        log_handle,
        (
            "Human BLAST search parameters: "
            f"blastp -db {database_path} -outfmt '6 sseqid pident length evalue bitscore' "
            f"-max_target_seqs {BLAST_MAX_TARGETS} -evalue 1e-10"
        ),
    )
    try:
        matches = find_homologous_frameworks(sequence, chain_type, database_path)
    except Exception as exc:
        _write_backmutation_log(
            log_handle,
            f"ERROR: failed to determine {chain_type} backmutated germline: {exc}",
        )
        return [f"backmt_germline_error {exc}"]

    if not matches:
        _write_backmutation_log(
            log_handle,
            f"No human germline match found for {chain_type} backmutated sequence",
        )
        return ["backmt_germline not found"]

    _write_backmutation_log(
        log_handle,
        f"Human germline candidates returned for {chain_type}: {len(matches)}; ranking criterion=bitscore descending",
    )
    for rank, match in enumerate(matches[:10], 1):
        _write_backmutation_log(
            log_handle,
            (
                f"Candidate rank {rank}: source={match['seq_id']}; "
                f"identity={match['identity']:.1f}%; bitscore={match['bitscore']:.1f}; "
                f"evalue={match['evalue']:.2e}; IMGT={format_germline_details(match)}"
            ),
        )
    if len(matches) > 10:
        _write_backmutation_log(
            log_handle,
            f"Additional candidates omitted from log: {len(matches) - 10}",
        )

    best_match = matches[0]
    _write_backmutation_log(
        log_handle,
        (
            f"Best {chain_type} backmutated germline match: "
            f"source={best_match['seq_id']}; identity={best_match['identity']:.1f}%; "
            f"IMGT={format_germline_details(best_match)}"
        ),
    )
    _write_backmutation_log(
        log_handle,
        f"Selection reason: rank #1 by highest BLAST bitscore ({best_match['bitscore']:.1f})",
    )
    return [
        f"  Source: {best_match['seq_id']}",
        f"  Identity: {best_match['identity']:.1f}%",
        f"  IMGT germline: {format_germline_details(best_match)}",
    ]


def enrich_backmutation_report(
    raw_output: str,
    database_path: str,
    log_handle: Optional[TextIO] = None,
    proposed_sequences: Optional[Dict[str, List[str]]] = None,
) -> str:
    """Add germline annotations determined from each backmutated sequence."""
    _write_backmutation_log(log_handle, "Enriching backmutation report with germline calls from backmt sequences")
    current_chain = None
    annotation_cache = {}
    sequence_counts = {"VH": 0, "VL": 0}
    rules_header_written = False
    rule_labels = {
        "murine_germline_difference": "Murine vs murine germline difference",
        "cysteine_proline": "Cysteine/Proline protection",
        "position_71": "Internal position 71 restoration",
        "fr4_motif": "FR4 motif correction",
    }
    proposed_sequences = proposed_sequences or {"VH": [], "VL": []}
    enriched_lines = [
        "BACKMUTATION RESULTS",
        "====================",
        "",
        "Germline annotations are recalculated from each backmt sequence using the human germline database.",
        f"Human germline database: {database_path}",
        "",
    ]
    for line in raw_output.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == "Chain VH":
            current_chain = "VH"
            _write_backmutation_log(log_handle, "Processing VH section from back.py output")
            enriched_lines.extend(["", "CHAIN VH", "========", ""])
            continue
        elif stripped == "Chain VL":
            current_chain = "VL"
            _write_backmutation_log(log_handle, "Processing VL section from back.py output")
            enriched_lines.extend(["", "CHAIN VL", "========", ""])
            continue

        if line.startswith("       "):
            enriched_lines.extend(["Framework map:", line.rstrip(), ""])
            continue

        if line.startswith("mGerm "):
            enriched_lines.extend(["Mouse germline alignment:", f"  mGerm: {line.replace('mGerm ', '', 1).strip()}"])
            continue

        if line.startswith("match "):
            enriched_lines.append(f"  match: {line.replace('match ', '', 1).strip()}")
            continue

        if line.startswith("input_ "):
            enriched_lines.extend([
                "",
                "Input sequence:",
                f"  input_: {line.replace('input_ ', '', 1).strip()}",
                "",
            ])
            continue

        if line.startswith("hGerm "):
            if current_chain:
                sequence_counts[current_chain] += 1
                rules_header_written = False
                enriched_lines.extend([
                    f"{current_chain} candidate {sequence_counts[current_chain]}",
                    "-" * (len(current_chain) + len(" candidate ") + len(str(sequence_counts[current_chain]))),
                    "Human germline template used for grafting:",
                    f"  hGerm: {line.replace('hGerm ', '', 1).strip()}",
                ])
            else:
                enriched_lines.append(line)
            continue

        if line.startswith("backmt ") and current_chain:
            backmutated_sequence = line.replace("backmt ", "", 1).strip()
            candidate_index = sequence_counts[current_chain] - 1
            chain_proposals = proposed_sequences.get(current_chain, [])
            _write_backmutation_log(
                log_handle,
                f"Found {current_chain} backmt sequence #{sequence_counts[current_chain]}; length={len(backmutated_sequence)}",
            )
            enriched_lines.extend([
                "Backmutated sequence:",
                f"  backmt: {backmutated_sequence}",
            ])
            if candidate_index < len(chain_proposals):
                changes = compare_backmutations(chain_proposals[candidate_index], backmutated_sequence)
                enriched_lines.append(f"Back-mutations performed: {len(changes)}")
                if changes:
                    enriched_lines.append("Changes (position: humanized -> backmutated):")
                    enriched_lines.extend(
                        f"  {change['position']}: {change['from']} -> {change['to']}"
                        for change in changes
                    )
                _write_backmutation_log(
                    log_handle,
                    f"{current_chain} candidate {sequence_counts[current_chain]} effective back-mutations={len(changes)}",
                )
            else:
                enriched_lines.append("Back-mutations performed: unavailable (previous humanized sequence missing)")
            enriched_lines.extend(["", "Germline determined from this backmt sequence:"])
            cache_key = (current_chain, backmutated_sequence)
            if cache_key not in annotation_cache:
                annotation_cache[cache_key] = _format_backmutated_germline_info(
                    backmutated_sequence,
                    current_chain,
                    database_path,
                    log_handle,
                )
            else:
                _write_backmutation_log(
                    log_handle,
                    f"Reusing cached germline call for repeated {current_chain} backmt sequence",
                )
            enriched_lines.extend(annotation_cache[cache_key])
            enriched_lines.append("")
            continue

        if line.startswith("rule ") and current_chain:
            rule_match = re.match(r"^rule\s+(\S+)\s+count=(\d+)\s+positions=(\S+)$", stripped)
            if rule_match:
                rule_name, count, positions = rule_match.groups()
                if not rules_header_written:
                    enriched_lines.append("Rules applied:")
                    enriched_lines.append("  Rule | Changes | Positions")
                    rules_header_written = True
                enriched_lines.append(
                    f"  {rule_labels.get(rule_name, rule_name)} | {count} | {positions}"
                )
                continue

        enriched_lines.append(line)

    _write_backmutation_log(
        log_handle,
        f"Finished enrichment; VH backmt sequences={sequence_counts['VH']}; VL backmt sequences={sequence_counts['VL']}",
    )
    return "\n".join(enriched_lines).rstrip() + "\n"


def compare_backmutations(proposed_sequence: str, backmutated_sequence: str) -> List[Dict[str, object]]:
    """Return effective residue changes from a grafted candidate to its backmutated sequence."""
    if len(proposed_sequence) != len(backmutated_sequence):
        raise ValueError(
            "Cannot count back-mutations in sequences with different lengths: "
            f"humanized={len(proposed_sequence)}, backmutated={len(backmutated_sequence)}"
        )
    return [
        {"position": position, "from": before, "to": after}
        for position, (before, after) in enumerate(zip(proposed_sequence, backmutated_sequence), start=1)
        if before != after
    ]


def load_backmutation_proposals(artifact_dir: Path) -> Dict[str, List[str]]:
    """Load graft-stage candidates in the same order consumed by back.py."""
    proposals = {"VH": [], "VL": []}
    for chain_type in proposals:
        fasta_path = artifact_dir / f"humanized_{chain_type.lower()}.fasta"
        if fasta_path.exists():
            with open(fasta_path, "r", encoding="utf-8") as fasta_handle:
                proposals[chain_type] = [str(record.seq) for record in SeqIO.parse(fasta_handle, "fasta")]
    return proposals


def backmutation_counts(raw_output: str, proposed_sequences: Dict[str, List[str]]) -> Dict[str, object]:
    """Count effective back-mutations per candidate and chain from back.py output."""
    counts: Dict[str, object] = {"VH": [], "VL": [], "total": 0}
    current_chain = None
    candidate_index = {"VH": 0, "VL": 0}
    for line in raw_output.splitlines():
        stripped = line.strip()
        if stripped in ("Chain VH", "Chain VL"):
            current_chain = stripped.split()[-1]
        elif line.startswith("backmt ") and current_chain:
            proposals = proposed_sequences.get(current_chain, [])
            index = candidate_index[current_chain]
            if index < len(proposals):
                count = len(compare_backmutations(proposals[index], line.replace("backmt ", "", 1).strip()))
                counts[current_chain].append(count)
                counts["total"] += count
            candidate_index[current_chain] += 1
    return counts


def _blast_db_exists(database_path: Path) -> bool:
    """Return True when a BLAST protein database exists for the given base path."""
    return any(Path(f"{database_path}{suffix}").exists() for suffix in (".pin", ".phr", ".psq", ".pal"))


def _write_backmutation_mouse_xml(
    sequence: str,
    chain_type: str,
    mouse_database: Path,
    working_dir: Path,
    log_handle: Optional[TextIO] = None,
) -> Path:
    """Create the XML BLAST output used by the original backmutation code."""
    query_file = working_dir / f"mouse_germline_{chain_type}.aa"
    xml_file = working_dir / f"mouse_germline_{chain_type}.xml"
    query_file.write_text(f"{sequence}\n", encoding="utf-8")
    _write_backmutation_log(
        log_handle,
        f"Wrote mouse germline BLAST query for {chain_type}: {query_file}; length={len(sequence)}",
    )

    cmd = [
        "blastp",
        "-db",
        str(mouse_database),
        "-max_target_seqs",
        "1",
        "-max_hsps",
        "1",
        "-outfmt",
        "5",
    ]
    _write_backmutation_log(
        log_handle,
        f"Running mouse germline BLAST for {chain_type}: {' '.join(cmd)} < {query_file} > {xml_file}",
    )
    try:
        with open(query_file, "r", encoding="utf-8") as query_handle, open(xml_file, "w", encoding="utf-8") as xml_handle:
            completed = subprocess.run(
                cmd,
                stdin=query_handle,
                stdout=xml_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                timeout=60,
            )
        _write_backmutation_log(
            log_handle,
            f"Mouse germline BLAST completed for {chain_type}; xml={xml_file}; bytes={xml_file.stat().st_size}",
        )
        if completed.stderr:
            _write_backmutation_log(log_handle, f"Mouse germline BLAST stderr for {chain_type}: {completed.stderr.strip()}")
    except subprocess.CalledProcessError as exc:
        _write_backmutation_log(
            log_handle,
            f"ERROR: mouse germline BLAST failed for {chain_type}: {exc.stderr}",
        )
        raise DatabaseError(f"Mouse germline BLAST failed for {chain_type}: {exc.stderr}") from exc
    except FileNotFoundError as exc:
        _write_backmutation_log(log_handle, "ERROR: command not found: blastp")
        raise DatabaseError("Command not found: blastp") from exc

    return xml_file


def run_backmutation(
    input_file: str,
    output_dir: str,
    vh_seq: str,
    vl_seq: str,
    mouse_database: str,
    human_database: str,
    artifact_dir: Optional[Path] = None,
) -> Path:
    """Run the original back/back.py workflow and save its raw output."""
    output_path = Path(output_dir)
    working_dir = output_path / "backmutation"
    working_dir.mkdir(parents=True, exist_ok=True)
    opt1_dir = Path(artifact_dir) if artifact_dir is not None else output_path
    log_file = output_path / "backmutation.log"

    mouse_db = Path(mouse_database)
    human_db = Path(human_database)

    with open(log_file, "w", encoding="utf-8") as backmutation_log:
        _write_backmutation_log(backmutation_log, "BACKMUTATION RUN START")
        _write_backmutation_log(backmutation_log, f"Input file: {input_file}")
        _write_backmutation_log(backmutation_log, f"Output directory: {output_path}")
        _write_backmutation_log(backmutation_log, f"Working directory: {working_dir}")
        _write_backmutation_log(backmutation_log, f"Opt1/artifact directory used by back.py: {opt1_dir}")
        _write_backmutation_log(backmutation_log, f"Original back.py script: {BACK_SCRIPT}")
        _write_backmutation_log(backmutation_log, f"Mouse germline database: {mouse_db}")
        _write_backmutation_log(backmutation_log, f"Human germline database for backmt reannotation: {human_db}")
        _write_backmutation_log(backmutation_log, f"Input VH length: {len(vh_seq)}")
        _write_backmutation_log(backmutation_log, f"Input VL length: {len(vl_seq)}")

        if not BACK_SCRIPT.exists():
            _write_backmutation_log(backmutation_log, f"ERROR: backmutation script not found: {BACK_SCRIPT}")
            raise DatabaseError(f"Backmutation script not found: {BACK_SCRIPT}")
        if not _blast_db_exists(mouse_db):
            _write_backmutation_log(backmutation_log, f"ERROR: mouse germline BLAST database not found: {mouse_db}")
            raise DatabaseError(
                f"Mouse germline BLAST database not found: {mouse_db}. "
                "Set MOUSE_GERMLINE_DB or pass --mouse-database."
            )
        _write_backmutation_log(backmutation_log, "Mouse germline BLAST database found")

        if not _blast_db_exists(human_db):
            _write_backmutation_log(backmutation_log, f"ERROR: human germline BLAST database not found: {human_db}")
            raise DatabaseError(f"Human germline BLAST database not found: {human_database}")
        _write_backmutation_log(backmutation_log, "Human germline BLAST database found")

        for expected_file in ("FRfileVH", "FRfileVL", "hGermVH", "hGermVL", "humanized_vh.fasta", "humanized_vl.fasta"):
            expected_path = opt1_dir / expected_file
            if expected_path.exists():
                _write_backmutation_log(backmutation_log, f"Artifact present: {expected_path}; bytes={expected_path.stat().st_size}")
            else:
                _write_backmutation_log(backmutation_log, f"WARNING: expected artifact missing before back.py run: {expected_path}")

        vh_xml = _write_backmutation_mouse_xml(vh_seq, "VH", mouse_db, working_dir, backmutation_log)
        vl_xml = _write_backmutation_mouse_xml(vl_seq, "VL", mouse_db, working_dir, backmutation_log)

        cmd = [
            sys.executable,
            str(BACK_SCRIPT),
            str(input_file),
            str(opt1_dir),
            str(vh_xml),
            str(vl_xml),
        ]
        _write_backmutation_log(backmutation_log, f"Running back.py command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=60)
        except subprocess.CalledProcessError as exc:
            _write_backmutation_log(backmutation_log, f"ERROR: back.py failed with return code {exc.returncode}")
            if exc.stdout:
                _write_backmutation_log(backmutation_log, f"back.py stdout before failure:\n{exc.stdout}")
            if exc.stderr:
                _write_backmutation_log(backmutation_log, f"back.py stderr before failure:\n{exc.stderr}")
            raise DatabaseError(f"Backmutation failed: {exc.stderr or exc.stdout}") from exc

        _write_backmutation_log(backmutation_log, "back.py completed successfully")
        _write_backmutation_log(backmutation_log, f"back.py stdout lines: {len(result.stdout.splitlines())}")
        if result.stderr:
            _write_backmutation_log(backmutation_log, f"back.py stderr:\n{result.stderr}")

        report_file = output_path / "backmutation.txt"
        proposed_sequences = load_backmutation_proposals(opt1_dir)
        mutation_counts = backmutation_counts(result.stdout, proposed_sequences)
        enriched_report = enrich_backmutation_report(
            result.stdout, human_database, backmutation_log, proposed_sequences
        )
        report_file.write_text(enriched_report, encoding="utf-8")
        _write_backmutation_log(backmutation_log, f"Wrote enriched backmutation report: {report_file}; bytes={report_file.stat().st_size}")

        # Isolated post-processing module: final backmutated sequences vs human germlines.
        from final_germline_alignment import run_final_germline_alignment
        final_germline_report = run_final_germline_alignment(
            result.stdout, human_database, output_dir, top_n=5
        )
        _write_backmutation_log(
            backmutation_log,
            f"Wrote final-sequence top-5 germline alignment report: {final_germline_report}",
        )

        if result.stderr:
            stderr_file = output_path / "backmutation.stderr.txt"
            stderr_file.write_text(result.stderr, encoding="utf-8")
            _write_backmutation_log(backmutation_log, f"Wrote back.py stderr file: {stderr_file}; bytes={stderr_file.stat().st_size}")

        _write_backmutation_log(backmutation_log, "BACKMUTATION RUN END")

    summary_file = output_path / "humanization_summary.txt"
    if summary_file.exists():
        with open(summary_file, "a", encoding="utf-8") as summary:
            summary.write("\nBackmutation\n")
            summary.write("===================\n")
            summary.write(f"Mouse germline database: {mouse_db}\n")
            summary.write(f"Human germline database: {human_database}\n")
            summary.write(f"Opt1 input folder: {opt1_dir}\n")
            summary.write(f"Raw output: {report_file.name}\n")
            summary.write(f"Backmutation log: {log_file.name}\n")
            summary.write(f"Final germline top-5 report: {final_germline_report.name}\n")
            summary.write("Final germline alignment log: final_germline_alignment.log\n")
            for chain_type in ("VH", "VL"):
                candidate_counts = mutation_counts[chain_type]
                summary.write(
                    f"{chain_type} back-mutations per candidate: "
                    f"{', '.join(map(str, candidate_counts)) if candidate_counts else 'unavailable'}\n"
                )
                summary.write(f"{chain_type} cumulative back-mutations: {sum(candidate_counts)}\n")
            summary.write(f"Total effective back-mutations: {mutation_counts['total']}\n")
            if result.stderr:
                summary.write("Warnings/errors: backmutation.stderr.txt\n")

    logger.info(f"💾 Saved backmutation output to {report_file}")
    logger.info(f"📝 Saved backmutation log to {log_file}")
    return report_file

def extract_structures_only(vh_seq: str, vl_seq: str) -> None:
    """Extract and display CDRs and frameworks from input sequences only."""
    log_info_header("STRUCTURE EXTRACTION MODE")
    
    # Use consistency manager for standardized extraction
    consistency_manager = get_consistency_manager()
    
    # Extract VH regions
    logger.info("📊 HEAVY CHAIN (VH) ANALYSIS:")
    logger.info("-" * 30)
    vh_result = consistency_manager.extract_regions_consistent(vh_seq, 'heavy')
    
    if vh_result.extraction_successful and vh_result.validation_passed:
        for region_name, region_seq in vh_result.regions.items():
            logger.info(f"  {region_name:>4}: {region_seq}")
        logger.info(f"  Total length: {len(vh_seq)} amino acids")
        logger.info(f"  Coverage: {vh_result.sequence_coverage:.1f}%")
    else:
        logger.error("❌ Could not extract VH regions")
        if vh_result.error_message:
            logger.error(f"   Error: {vh_result.error_message}")
    
    logger.info("")
    
    # Extract VL regions
    logger.info("📊 LIGHT CHAIN (VL) ANALYSIS:")
    logger.info("-" * 30)
    vl_result = consistency_manager.extract_regions_consistent(vl_seq, 'light')
    
    if vl_result.extraction_successful and vl_result.validation_passed:
        for region_name, region_seq in vl_result.regions.items():
            logger.info(f"  {region_name:>4}: {region_seq}")
        logger.info(f"  Total length: {len(vl_seq)} amino acids")
        logger.info(f"  Coverage: {vl_result.sequence_coverage:.1f}%")
    else:
        logger.error("❌ Could not extract VL regions")
        if vl_result.error_message:
            logger.error(f"   Error: {vl_result.error_message}")

def show_database_scores(vh_seq: str, vl_seq: str, database_path: str) -> None:
    """Show structures + search database for best V, D, J candidates with scores."""
    log_info_header("DATABASE SCORING MODE")
    
    # First show structures
    extract_structures_only(vh_seq, vl_seq)
    
    # Validate CDR extraction consistency
    log_info_header("CDR CONSISTENCY VALIDATION")
    _validate_and_report_cdr_consistency(vh_seq, vl_seq)
    
    logger.info("")
    log_info_header("DATABASE SEARCH RESULTS")
    
    # Search for VH candidates
    logger.info("📊 HEAVY CHAIN (VH) CANDIDATES:")
    logger.info("-" * 40)
    vh_candidates = find_homologous_frameworks(vh_seq, 'heavy', database_path)
    display_candidates(vh_candidates, 'VH')
    
    logger.info("")
    
    # Search for VL candidates  
    logger.info("📊 LIGHT CHAIN (VL) CANDIDATES:")
    logger.info("-" * 40)
    vl_candidates = find_homologous_frameworks(vl_seq, 'light', database_path)
    display_candidates(vl_candidates, 'VL')

def apply_optimization(vh_candidates: List[Dict], vl_candidates: List[Dict], level: int, 
                      murine_vh: str, murine_vl: str, 
                      murine_vh_cdrs: Dict[str, str], murine_vl_cdrs: Dict[str, str]) -> Tuple[List[Dict], List[Dict]]:
    """
    Apply optimization level to humanized candidates using intelligent candidate processing.
    
    Performance optimized to avoid redundant optimization calls by:
    1. Using representative candidates for optimization
    2. Applying optimizations strategically based on sequence similarity
    3. Leveraging ANARCII caching for maximum efficiency
    """
    logger.info(f"🔧 APPLYING OPTIMIZATION LEVEL {level}")
    logger.info("=" * 50)
    
    try:
        from optimizations import OptimizationEngine
        
        # Initialize optimization engine once
        optimizer = OptimizationEngine()
        
        # Performance optimization: Use representative candidates instead of optimizing all
        # This dramatically reduces computation while maintaining scientific accuracy
        
        # Step 1: Select representative candidates (best identity scores)
        if vh_candidates:
            representative_vh = max(vh_candidates, key=lambda x: x['identity'])
            logger.info(f"🎯 Selected representative VH candidate: {representative_vh['source']} ({representative_vh['identity']:.1f}% identity)")
        else:
            logger.warning("⚠️  No VH candidates available for optimization")
            return [], vl_candidates
            
        if vl_candidates:
            representative_vl = max(vl_candidates, key=lambda x: x['identity'])
            logger.info(f"🎯 Selected representative VL candidate: {representative_vl['source']} ({representative_vl['identity']:.1f}% identity)")
        else:
            logger.warning("⚠️  No VL candidates available for optimization")
            return vh_candidates, []
        
        # Step 2: Perform single optimization on representative pair
        logger.info(f"🧬 Performing optimization on representative candidate pair...")
        start_time = time.time()
        
        result = optimizer.optimize_sequences(
            vh_sequence=representative_vh['sequence'],
            vl_sequence=representative_vl['sequence'],
            murine_vh=murine_vh,
            murine_vl=murine_vl,
            murine_vh_cdrs=murine_vh_cdrs,
            murine_vl_cdrs=murine_vl_cdrs,
            optimization_level=level
        )
        
        optimization_time = time.time() - start_time
        logger.info(f"✅ Representative optimization completed in {optimization_time:.2f}s")
        
        # Step 3: Apply optimization results to all candidates
        optimized_vh_candidates = []
        optimized_vl_candidates = []
        
        if result.success:
            logger.info(f"🎯 Applying optimizations to all candidates...")
            
            # Calculate the optimization differences
            vh_optimization_diff = _calculate_sequence_diff(representative_vh['sequence'], result.optimized_vh_sequence)
            vl_optimization_diff = _calculate_sequence_diff(representative_vl['sequence'], result.optimized_vl_sequence)
            
            logger.info(f"   VH optimization changes: {len(vh_optimization_diff)} positions")
            logger.info(f"   VL optimization changes: {len(vl_optimization_diff)} positions")
            
            # Apply optimizations to all VH candidates
            for i, vh_candidate in enumerate(vh_candidates):
                try:
                    optimized_sequence = _apply_optimization_diff(vh_candidate['sequence'], vh_optimization_diff)
                    
                    optimized_candidate = vh_candidate.copy()
                    optimized_candidate['sequence'] = optimized_sequence
                    optimized_candidate['optimization_score'] = result.vh_improvement_score
                    optimized_candidate['optimization_applied'] = True
                    optimized_candidate['optimization_level'] = level
                    optimized_candidate['optimization_method'] = 'representative_based'
                    optimized_vh_candidates.append(optimized_candidate)
                    
                    logger.debug(f"✅ VH candidate {i+1} optimized via representative method")
                    
                except Exception as e:
                    logger.warning(f"⚠️  VH candidate {i+1} optimization failed: {e}")
                    vh_candidate['optimization_applied'] = False
                    vh_candidate['optimization_score'] = 0.0
                    vh_candidate['optimization_level'] = level
                    optimized_vh_candidates.append(vh_candidate)
            
            # Apply optimizations to all VL candidates
            for i, vl_candidate in enumerate(vl_candidates):
                try:
                    optimized_sequence = _apply_optimization_diff(vl_candidate['sequence'], vl_optimization_diff)
                    
                    optimized_candidate = vl_candidate.copy()
                    optimized_candidate['sequence'] = optimized_sequence
                    optimized_candidate['optimization_score'] = result.vl_improvement_score
                    optimized_candidate['optimization_applied'] = True
                    optimized_candidate['optimization_level'] = level
                    optimized_candidate['optimization_method'] = 'representative_based'
                    optimized_vl_candidates.append(optimized_candidate)
                    
                    logger.debug(f"✅ VL candidate {i+1} optimized via representative method")
                    
                except Exception as e:
                    logger.warning(f"⚠️  VL candidate {i+1} optimization failed: {e}")
                    vl_candidate['optimization_applied'] = False
                    vl_candidate['optimization_score'] = 0.0
                    vl_candidate['optimization_level'] = level
                    optimized_vl_candidates.append(vl_candidate)
        
        else:
            log_warning("Representative optimization failed, returning original candidates")
            mark_candidates_as_unoptimized(vh_candidates, vl_candidates, level)
            return vh_candidates, vl_candidates
        
        # Summary
        vh_optimized = sum(1 for c in optimized_vh_candidates if c.get('optimization_applied', False))
        vl_optimized = sum(1 for c in optimized_vl_candidates if c.get('optimization_applied', False))
        
        logger.info(f"✅ Optimization completed:")
        logger.info(f"   VH: {vh_optimized}/{len(vh_candidates)} candidates optimized")
        logger.info(f"   VL: {vl_optimized}/{len(vl_candidates)} candidates optimized")
        logger.info(f"⚡ Performance gain: ~{(len(vh_candidates) + len(vl_candidates) - 1) * optimization_time:.1f}s saved")
        
        return optimized_vh_candidates, optimized_vl_candidates
        
    except ImportError as e:
        log_error(f"Failed to import optimization engine: {e}")
        log_warning("Returning original candidates without optimization")
        mark_candidates_as_unoptimized(vh_candidates, vl_candidates, level)
        return vh_candidates, vl_candidates
    
    except Exception as e:
        log_error(f"Optimization level {level} failed: {e}")
        log_warning("Returning original candidates without optimization")
        mark_candidates_as_unoptimized(vh_candidates, vl_candidates, level)
            
        return vh_candidates, vl_candidates

def _calculate_sequence_diff(original: str, optimized: str) -> Dict[int, str]:
    """
    Calculate the differences between original and optimized sequences.
    
    Args:
        original: Original sequence
        optimized: Optimized sequence
        
    Returns:
        Dictionary mapping position (0-based) to new amino acid
    """
    if len(original) != len(optimized):
        logger.warning(f"⚠️  Sequence length mismatch: {len(original)} vs {len(optimized)}")
        return {}
    
    differences = {}
    for i, (orig_aa, opt_aa) in enumerate(zip(original, optimized)):
        if orig_aa != opt_aa:
            differences[i] = opt_aa
    
    return differences

def _apply_optimization_diff(sequence: str, optimization_diff: Dict[int, str]) -> str:
    """
    Apply optimization differences to a sequence.
    
    Args:
        sequence: Original sequence to modify
        optimization_diff: Dictionary of position -> new amino acid
        
    Returns:
        Modified sequence with optimizations applied
    """
    if not optimization_diff:
        return sequence
    
    sequence_list = list(sequence)
    
    for position, new_aa in optimization_diff.items():
        if position < len(sequence_list):
            sequence_list[position] = new_aa
        else:
            logger.warning(f"⚠️  Position {position} out of bounds for sequence length {len(sequence_list)}")
    
    return ''.join(sequence_list)

def load_sequences(input_file: str) -> Tuple[str, str]:
    """
    Load and validate VH and VL sequences from input file.
    
    Args:
        input_file: Path to input file containing sequences
        
    Returns:
        Tuple of (vh_sequence, vl_sequence)
        
    Raises:
        SequenceValidationError: If sequences are invalid
        FileNotFoundError: If input file doesn't exist
    """
    try:
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        if len(lines) < 2:
            raise SequenceValidationError("Input file must contain at least 2 sequences (VH and VL)")
        
        # Validate and clean sequences using robust validation
        vh_sequence = validate_sequence(lines[0], min_length=50)
        vl_sequence = validate_sequence(lines[1], min_length=50)
        
        logger.info(f"✅ Loaded sequences: VH ({len(vh_sequence)} AA), VL ({len(vl_sequence)} AA)")
        
        return vh_sequence, vl_sequence
        
    except Exception as e:
        logger.error(f"❌ Failed to load sequences: {e}")
        raise

def _validate_and_report_cdr_consistency(vh_seq: str, vl_seq: str) -> None:
    """
    Validate and report CDR extraction consistency across different chain type specifications.
    
    This function tests whether the same sequence analyzed with different chain type
    parameters produces consistent CDR boundaries, which is critical for reliable
    humanization.
    """
    consistency_manager = get_consistency_manager()
    
    # Test VH consistency between different chain type specifications
    logger.info("🔍 Testing VH CDR consistency...")
    vh_consistency = consistency_manager.validate_cdr_consistency(
        vh_seq, 'heavy', vh_seq, 'VH'
    )
    
    if vh_consistency['consistent']:
        logger.info("✅ VH CDR extraction is consistent across chain type specifications")
    else:
        logger.warning("⚠️  VH CDR extraction inconsistency detected:")
        for inconsistency in vh_consistency['inconsistent_cdrs']:
            logger.warning(f"   {inconsistency['cdr']}: '{inconsistency['extraction1']}' vs '{inconsistency['extraction2']}'")
            logger.warning(f"   Chain types: {inconsistency['chain_type1']} vs {inconsistency['chain_type2']}")
    
    # Test VL consistency between different chain type specifications  
    logger.info("🔍 Testing VL CDR consistency...")
    vl_consistency = consistency_manager.validate_cdr_consistency(
        vl_seq, 'light', vl_seq, 'VL'
    )
    
    if vl_consistency['consistent']:
        logger.info("✅ VL CDR extraction is consistent across chain type specifications")
    else:
        logger.warning("⚠️  VL CDR extraction inconsistency detected:")
        for inconsistency in vl_consistency['inconsistent_cdrs']:
            logger.warning(f"   {inconsistency['cdr']}: '{inconsistency['extraction1']}' vs '{inconsistency['extraction2']}'")
            logger.warning(f"   Chain types: {inconsistency['chain_type1']} vs {inconsistency['chain_type2']}")
    
    # Report overall consistency status
    overall_consistent = vh_consistency['consistent'] and vl_consistency['consistent']
    if overall_consistent:
        logger.info("🎉 Overall CDR extraction consistency: PASSED")
    else:
        logger.warning("⚠️  Overall CDR extraction consistency: ISSUES DETECTED")
        logger.warning("   This may affect humanization accuracy. Consider reviewing CDR boundaries.")

def preprocess_x_characters(sequence: str) -> str:
    """
    Preprocess sequences with X characters using scientifically-informed domain boundary recognition.
    
    Based on antibody domain structure analysis:
    - X in CDR3 context (CAR[X]GTT) represents D-segment diversity → use Glycine (most common)
    - X after VH domain boundary (WGQGTLVTVSS) → TRUNCATE at domain boundary (scientifically justified)
    - Other X characters → use Alanine (conservative)
    
    Scientific rationale:
    - VH domain has defined C-terminal boundary at WGQGTLVTVSS
    - Sequences beyond this boundary represent CH1 domain or artifacts
    - For VH humanization, we work with VH domains only
    """
    if 'X' not in sequence:
        return sequence
    
    processed = sequence
    
    # Strategy 1: CDR3 D-segment substitution (position-specific)
    # Pattern: CAR[X]GTT - replace X with G (glycine, most common in D-segments)
    import re
    cdr3_pattern = r'CAR(X+)GTT'
    matches = re.finditer(cdr3_pattern, processed)
    for match in matches:
        x_segment = match.group(1)
        # Replace X's in CDR3 with G (glycine - most common D-segment AA)
        g_replacement = 'G' * len(x_segment)
        processed = processed.replace(f'CAR{x_segment}GTT', f'CAR{g_replacement}GTT')
    
    # Strategy 2: VH domain boundary recognition
    # If sequence contains VH C-terminus followed by X, truncate at domain boundary
    vh_endings = ['WGQGTLVTVSS', 'WGQGTSVTVSS']  # Common VH C-terminal patterns
    for ending in vh_endings:
        if ending in processed:
            end_pos = processed.find(ending) + len(ending)
            # Check if there are X characters after this position
            remainder = processed[end_pos:]
            if 'X' in remainder:
                # Truncate at VH domain boundary (scientifically justified)
                processed = processed[:end_pos]
                logger.debug(f"Truncated at VH domain boundary: {ending}")
                break
    
    # Strategy 3: Any remaining X characters (internal positions)
    # Replace with A (alanine - conservative)
    processed = processed.replace('X', 'A')
    
    return processed

def find_homologous_frameworks(sequence: str, chain_type: str, database_path: str) -> List[Dict]:
    """Find homologous frameworks using BLAST."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as query_file:
        query_file.write(f">query\n{sequence}\n")
        query_file_path = query_file.name
    
    try:
        # First, get BLAST hits (just IDs and scores)
        cmd = [
            'blastp', '-query', query_file_path, '-db', database_path,
            '-outfmt', '6 sseqid pident length evalue bitscore',
            '-max_target_seqs', str(BLAST_MAX_TARGETS), '-evalue', '1e-10'
        ]
        result = safe_execute_subprocess(cmd, timeout=60)
        
        matches = []
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 5:
                seq_id, pident, length, evalue, bitscore = parts
                
                # Get the FULL sequence from the database using the seq_id
                full_seq_cmd = [
                    'blastdbcmd', '-db', database_path, '-entry', seq_id,
                    '-outfmt', '%s\t%t'
                ]
                full_seq_result = subprocess.run(full_seq_cmd, capture_output=True, text=True, check=True)
                sequence_and_title = full_seq_result.stdout.strip().split('\t', 1)
                full_sequence = sequence_and_title[0]
                sequence_title = sequence_and_title[1] if len(sequence_and_title) > 1 else ''
                
                # Preprocess X characters for ANARCII compatibility
                if full_sequence:
                    full_sequence = preprocess_x_characters(full_sequence)
                
                if full_sequence:  # Only add if we got the full sequence
                    matches.append({
                        'seq_id': seq_id,
                        'identity': float(pident),
                        'sequence': full_sequence,  # Now using FULL sequence
                        'evalue': float(evalue),
                        'bitscore': float(bitscore),
                        'germline': parse_germline_metadata(seq_id, sequence_title),
                    })
        
        return sorted(matches, key=lambda x: x['bitscore'], reverse=True)[:50]
        
    finally:
        os.unlink(query_file_path)

def graft_cdrs(
    murine_seq: str,
    human_seq: str,
    chain_type: str,
    anarcii: ANARCIICDRExtractor,
    backmutation_artifact_dir: Optional[Path] = None,
) -> Optional[str]:
    """Graft murine CDRs onto human framework with exact CDR preservation."""
    try:
        # Use consistency manager for standardized extraction
        consistency_manager = get_consistency_manager()
        
        # Extract regions from both sequences using consistent analysis
        murine_result = consistency_manager.extract_regions_consistent(murine_seq, chain_type)
        human_result = consistency_manager.extract_regions_consistent(human_seq, chain_type)
        
        # Validate both extractions succeeded
        if not (murine_result.extraction_successful and murine_result.validation_passed):
            logger.debug(f"Murine sequence extraction failed: {murine_result.error_message}")
            return None
            
        if not (human_result.extraction_successful and human_result.validation_passed):
            logger.debug(f"Human sequence extraction failed: {human_result.error_message}")
            return None
        
        murine_regions = murine_result.regions
        human_regions = human_result.regions

        # Graft: Human frameworks + Murine CDRs
        grafted = (
            human_regions.get('FR1', '') +
            murine_regions.get('CDR1', '') +
            human_regions.get('FR2', '') +
            murine_regions.get('CDR2', '') +
            human_regions.get('FR3', '') +
            murine_regions.get('CDR3', '') +
            human_regions.get('FR4', '')
        )
        
        if not grafted:
            return None
            
        # CRITICAL: Validate CDR preservation using robust position-based validation
        # This ensures exact CDR preservation regardless of ANARCII re-analysis context changes
        if not _validate_cdr_preservation_robust(murine_regions, grafted, human_regions):
            logger.debug(f"CDR preservation validation failed")
            return None
        
        # Additional consistency validation: ensure the grafted sequence maintains CDR integrity
        grafted_result = consistency_manager.extract_regions_consistent(grafted, chain_type)
        if grafted_result.extraction_successful:
            # Verify CDRs are preserved in the grafted sequence
            for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
                original_cdr = murine_regions.get(cdr_name, '')
                grafted_cdr = grafted_result.regions.get(cdr_name, '')
                if original_cdr and original_cdr != grafted_cdr:
                    logger.debug(f"CDR {cdr_name} not preserved in grafted sequence: {original_cdr} → {grafted_cdr}")
                    return None
        
        write_backmutation_graft_artifacts(murine_regions, human_regions, chain_type, backmutation_artifact_dir)
        
        return grafted
        
    except Exception as e:
        logger.debug(f"CDR grafting failed: {e}")
        return None

def _validate_cdr_preservation_robust(murine_regions: dict, grafted_seq: str, human_regions: dict) -> bool:
    """
    Validate CDR preservation using robust position-based validation.
    
    This method ensures exact CDR preservation regardless of how ANARCII might
    re-analyze the grafted sequence in different framework contexts.
    
    Scientific rationale:
    - CDR grafting should preserve exact murine CDR sequences
    - Framework context should not affect CDR boundaries for validation
    - Position-based validation is more reliable than re-analysis
    """
    try:
        # Build expected sequence components
        expected_components = [
            human_regions.get('FR1', ''),
            murine_regions.get('CDR1', ''),
            human_regions.get('FR2', ''),
            murine_regions.get('CDR2', ''),
            human_regions.get('FR3', ''),
            murine_regions.get('CDR3', ''),
            human_regions.get('FR4', '')
        ]
        
        # Check that grafted sequence matches expected construction
        expected_seq = ''.join(expected_components)
        if grafted_seq != expected_seq:
            logger.debug(f"Grafted sequence doesn't match expected construction")
            logger.debug(f"Expected: {expected_seq}")
            logger.debug(f"Actual:   {grafted_seq}")
            return False
            
        # Robust CDR validation using multiple approaches
        return (
            _validate_cdr_by_position(murine_regions, grafted_seq, human_regions) and
            _validate_cdr_by_pattern_matching(murine_regions, grafted_seq) and
            _validate_sequence_integrity(grafted_seq)
        )
        
    except Exception as e:
        logger.debug(f"CDR preservation validation error: {e}")
        return False

def _validate_cdr_by_position(murine_regions: dict, grafted_seq: str, human_regions: dict) -> bool:
    """Validate CDRs are in exact expected positions."""
    try:
        current_pos = 0
        
        # Skip FR1
        current_pos += len(human_regions.get('FR1', ''))
        
        # Check CDR1
        cdr1 = murine_regions.get('CDR1', '')
        if grafted_seq[current_pos:current_pos + len(cdr1)] != cdr1:
            logger.debug(f"CDR1 position validation failed: expected {cdr1}, got {grafted_seq[current_pos:current_pos + len(cdr1)]}")
            return False
        current_pos += len(cdr1)
        
        # Skip FR2
        current_pos += len(human_regions.get('FR2', ''))
        
        # Check CDR2
        cdr2 = murine_regions.get('CDR2', '')
        if grafted_seq[current_pos:current_pos + len(cdr2)] != cdr2:
            logger.debug(f"CDR2 position validation failed: expected {cdr2}, got {grafted_seq[current_pos:current_pos + len(cdr2)]}")
            return False
        current_pos += len(cdr2)
        
        # Skip FR3
        current_pos += len(human_regions.get('FR3', ''))
        
        # Check CDR3 - MOST CRITICAL
        cdr3 = murine_regions.get('CDR3', '')
        if grafted_seq[current_pos:current_pos + len(cdr3)] != cdr3:
            logger.debug(f"CDR3 position validation failed: expected {cdr3}, got {grafted_seq[current_pos:current_pos + len(cdr3)]}")
            return False
            
        return True
        
    except Exception as e:
        logger.debug(f"Position validation error: {e}")
        return False

def _validate_cdr_by_pattern_matching(murine_regions: dict, grafted_seq: str) -> bool:
    """Validate CDRs exist as exact substrings in the grafted sequence."""
    try:
        for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
            cdr_seq = murine_regions.get(cdr_name, '')
            if cdr_seq and cdr_seq not in grafted_seq:
                logger.debug(f"{cdr_name} pattern matching failed: {cdr_seq} not found in grafted sequence")
                return False
        return True
    except Exception as e:
        logger.debug(f"Pattern matching validation error: {e}")
        return False

def _validate_sequence_integrity(grafted_seq: str) -> bool:
    """Validate the grafted sequence has proper amino acid composition."""
    try:
        # Check for valid amino acids only
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aa for aa in grafted_seq):
            logger.debug("Invalid amino acids in grafted sequence")
            return False
        
        # Check reasonable length (should be > 50 AA for antibody domains)
        if len(grafted_seq) < 50:
            logger.debug(f"Grafted sequence too short: {len(grafted_seq)} AA")
            return False
            
        return True
    except Exception as e:
        logger.debug(f"Sequence integrity validation error: {e}")
        return False

def _validate_cdr_preservation_direct(murine_regions: dict, grafted_seq: str, human_regions: dict) -> bool:
    """Compatibility validation wrapper."""
    return _validate_cdr_preservation_robust(murine_regions, grafted_seq, human_regions)

def deduplicate_database_frameworks(blast_hits: List[Dict], chain_type: str, max_candidates: int = 5) -> List[Dict]:
    """
    Deduplicate database frameworks BEFORE CDR grafting for optimal performance.
    
    This function:
    1. Extracts framework regions from original database sequences
    2. Identifies unique framework architectures
    3. Returns top unique database hits for grafting
    4. Preserves original BLAST ranking order
    
    Args:
        blast_hits: List of BLAST hits with 'sequence', 'identity', 'seq_id'
        chain_type: 'VH' or 'VL'
        max_candidates: Maximum unique frameworks to return for grafting
    
    Returns:
        List of unique database hits (same format as input)
    """
    if not blast_hits:
        return blast_hits
    
    try:
        consistency_manager = get_consistency_manager()
        seen_frameworks = set()
        unique_hits = []
        
        logger.debug(f"Deduplicating {len(blast_hits)} {chain_type} database frameworks...")
        
        for hit in blast_hits:
            try:
                # Extract framework regions from the database sequence
                result = consistency_manager.extract_regions_consistent(
                    hit['sequence'], 
                    chain_type
                )
                
                if not (result.extraction_successful and result.validation_passed):
                    logger.debug(f"Framework extraction failed for {hit['seq_id']}: {result.error_message}")
                    continue
                
                # Create framework signature (concatenated framework regions)
                regions = result.regions
                framework_signature = (
                    regions.get('FR1', '') +
                    regions.get('FR2', '') +
                    regions.get('FR3', '') +
                    regions.get('FR4', '')
                )
                
                # Check if this framework is unique
                if framework_signature not in seen_frameworks:
                    seen_frameworks.add(framework_signature)
                    unique_hits.append(hit)
                    logger.debug(f"Unique framework found: {hit['seq_id']} (framework length: {len(framework_signature)} AA)")
                    
                    if len(unique_hits) >= max_candidates:
                        break
                else:
                    logger.debug(f"Duplicate framework skipped: {hit['seq_id']}")
                    
            except Exception as e:
                logger.debug(f"Error processing hit {hit.get('seq_id', 'unknown')}: {str(e)}")
                continue
        
        logger.info(f"   Framework deduplication: {len(blast_hits)} database hits → {len(unique_hits)} unique frameworks selected for grafting")
        return unique_hits
        
    except Exception as e:
        logger.warning(f"Deduplication failed for {chain_type}: {str(e)}")
        # Fallback: return original hits (limited to max_candidates)
        return blast_hits[:max_candidates]


def deduplicate_candidates_by_framework(candidates: List[Dict], murine_seq: str, chain_type: str, max_candidates: int = 5) -> List[Dict]:
    """
    Deduplicate candidates based on framework sequences after CDR grafting.
    
    This function:
    1. Extracts framework regions from each grafted candidate
    2. Identifies duplicates based on framework sequence identity
    3. Returns only unique framework architectures
    4. Preserves original BLAST ranking order
    
    Args:
        candidates: List of grafted candidates with 'sequence', 'identity', 'source', 'length'
        murine_seq: Original murine sequence (for CDR reference)
        chain_type: 'VH' or 'VL'
        max_candidates: Maximum unique candidates to return
    
    Returns:
        List of unique candidates (same format as input)
    """
    if not candidates:
        return candidates
    
    try:
        consistency_manager = get_consistency_manager()
        seen_frameworks = set()
        unique_candidates = []
        
        logger.debug(f"Deduplicating {len(candidates)} {chain_type} candidates...")
        
        for candidate in candidates:
            try:
                # Extract framework regions from the grafted sequence
                result = consistency_manager.extract_regions_consistent(
                    candidate['sequence'], 
                    chain_type
                )
                
                if not (result.extraction_successful and result.validation_passed):
                    logger.debug(f"Framework extraction failed for {candidate['source']}: {result.error_message}")
                    continue
                
                # Create framework signature (concatenated framework regions)
                regions = result.regions
                framework_signature = (
                    regions.get('FR1', '') +
                    regions.get('FR2', '') +
                    regions.get('FR3', '') +
                    regions.get('FR4', '')
                )
                
                # Check if this framework is unique
                if framework_signature not in seen_frameworks:
                    seen_frameworks.add(framework_signature)
                    unique_candidates.append(candidate)
                    logger.debug(f"Unique framework found: {candidate['source']} (framework length: {len(framework_signature)} AA)")
                    
                    if len(unique_candidates) >= max_candidates:
                        break
                else:
                    logger.debug(f"Duplicate framework skipped: {candidate['source']}")
                    
            except Exception as e:
                logger.debug(f"Error processing candidate {candidate.get('source', 'unknown')}: {str(e)}")
                continue
        
        logger.info(f"   Framework deduplication: {len(candidates)} → {len(unique_candidates)} unique {chain_type} candidates")
        return unique_candidates
        
    except Exception as e:
        logger.warning(f"Deduplication failed for {chain_type}: {str(e)}")
        # Fallback: return original candidates
        return candidates[:max_candidates]


def humanize_antibody(
    vh_seq: str,
    vl_seq: str,
    database_path: str,
    backmutation_artifact_dir: Optional[Path] = None,
) -> Dict[str, List]:
    """Humanize antibody using VDJ germline database."""
    anarcii = ANARCIICDRExtractor()
    
    logger.info("🧬 VDJ-BASED ANTIBODY HUMANIZATION")
    logger.info("=" * 60)
    logger.info(f"Input: VH={len(vh_seq)} AA, VL={len(vl_seq)} AA")
    logger.info(f"Database: {database_path}")
    
    results = {'VH': [], 'VL': []}
    
    # Process VH
    logger.info("🔍 Processing VH chain...")
    vh_matches = find_homologous_frameworks(vh_seq, 'VH', database_path)
    logger.info(f"   Found {len(vh_matches)} BLAST hits, deduplicating frameworks before grafting...")
    
    # OPTIMIZED: Deduplicate frameworks BEFORE grafting
    unique_vh_hits = deduplicate_database_frameworks(vh_matches, 'VH', MAX_GENERATED_CANDIDATES)
    logger.info(f"   Selected {len(unique_vh_hits)} unique frameworks for CDR grafting...")
    
    # Graft CDRs only on unique frameworks
    vh_candidates = []
    for i, match in enumerate(unique_vh_hits):
        try:
            grafted = graft_cdrs(vh_seq, match['sequence'], 'VH', anarcii, backmutation_artifact_dir)
            if grafted:
                vh_candidates.append({
                    'sequence': grafted,
                    'identity': match['identity'],
                    'source': match['seq_id'],
                    'length': len(grafted),
                    'germline': match.get('germline', {}),
                })
                logger.info(
                    f"   ✅ Successfully grafted using unique framework {i+1}/{len(unique_vh_hits)} "
                    f"(ID: {match['seq_id']}, identity: {match['identity']:.1f}%, "
                    f"IMGT: {format_germline_details(match)})"
                )
        except Exception as e:
            logger.debug(f"   ⚠️  Framework {i+1}/{len(unique_vh_hits)} failed: {str(e)[:100]}...")
            continue
    
    results['VH'] = vh_candidates
    
    # Process VL
    logger.info("🔍 Processing VL chain...")
    vl_matches = find_homologous_frameworks(vl_seq, 'VL', database_path)
    logger.info(f"   Found {len(vl_matches)} BLAST hits, deduplicating frameworks before grafting...")
    
    # OPTIMIZED: Deduplicate frameworks BEFORE grafting
    unique_vl_hits = deduplicate_database_frameworks(vl_matches, 'VL', MAX_GENERATED_CANDIDATES)
    logger.info(f"   Selected {len(unique_vl_hits)} unique frameworks for CDR grafting...")
    
    # Graft CDRs only on unique frameworks
    vl_candidates = []
    for i, match in enumerate(unique_vl_hits):
        try:
            grafted = graft_cdrs(vl_seq, match['sequence'], 'VL', anarcii, backmutation_artifact_dir)
            if grafted:
                vl_candidates.append({
                    'sequence': grafted,
                    'identity': match['identity'],
                    'source': match['seq_id'],
                    'length': len(grafted),
                    'germline': match.get('germline', {}),
                })
                logger.info(
                    f"   ✅ Successfully grafted using unique framework {i+1}/{len(unique_vl_hits)} "
                    f"(ID: {match['seq_id']}, identity: {match['identity']:.1f}%, "
                    f"IMGT: {format_germline_details(match)})"
                )
        except Exception as e:
            logger.debug(f"   ⚠️  Framework {i+1}/{len(unique_vl_hits)} failed: {str(e)[:100]}...")
            continue
    
    results['VL'] = vl_candidates
    
    # Log results
    logger.info(f"🎉 VDJ Humanization Results:")
    logger.info(f"   VH candidates: {len(results['VH'])}")
    logger.info(f"   VL candidates: {len(results['VL'])}")
    
    return results

def save_results(results: Dict[str, List], output_dir: str) -> None:
    """Save humanization results to files."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save sequences using utility function
    for chain_type, candidates in results.items():
        write_candidate_fasta(candidates, output_path, chain_type)
        if candidates:
            logger.info(f"💾 Saved {len(candidates)} {chain_type} candidates to {output_path}/humanized_{chain_type.lower()}.fasta")
    
    # Save summary using utility function
    with open(output_path / "humanization_summary.txt", 'w') as f:
        f.write("ANTIBODY HUMANIZATION SUMMARY\n")
        f.write("=" * 40 + "\n\n")
        
        total_candidates = sum(len(candidates) for candidates in results.values())
        f.write(f"Total candidates generated: {total_candidates}\n")
        
        for chain_type, candidates in results.items():
            write_candidate_summary(candidates, f, chain_type)
    
    logger.info(f"📋 Saved summary to {output_path}/humanization_summary.txt")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='VDJ-based Antibody Humanizer\n\nDefault behavior: Full optimization pipeline (level 4) for maximum therapeutic quality.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('input_file', help='Input file with VH and VL sequences (one per line)')
    parser.add_argument('-o', '--output', default=None, help='Output directory')
    parser.add_argument('--structures', action='store_true', help='Extract CDRs and frameworks from input sequences only')
    parser.add_argument('--scores', action='store_true', help='Show structures + search database for best V, D, J candidates with scores')
    parser.add_argument('--graft', action='store_true', help='Perform CDR grafting only (no optimization)')
    parser.add_argument('--optimization', type=int, choices=[1, 2, 3, 4], help='Apply optimization level (1-4)')
    parser.add_argument('--database', default='imgt_germline_database/human_germline_frameworks', 
                       help='BLAST database path')
    parser.add_argument(
        '--backmutation',
        action='store_true',
        help='Run the original back/back.py backmutation workflow after humanization'
    )
    parser.add_argument(
        '--mouse-database',
        default=str(DEFAULT_MOUSE_GERMLINE_DB),
        help='Mouse germline BLAST database base path for backmutation'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    if args.output:
        output_dir = args.output
    else:
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = f"results_{timestamp}"
    
    try:
        # Load sequences
        vh_seq, vl_seq = load_sequences(args.input_file)
        logger.info(f"✅ Loaded sequences from {args.input_file}")
        logger.info(f"   VH: {len(vh_seq)} amino acids")
        logger.info(f"   VL: {len(vl_seq)} amino acids")
        
        # Handle different modes
        if args.structures:
            # Structure extraction only
            extract_structures_only(vh_seq, vl_seq)
            return 0
            
        elif args.scores:
            # Structure + database scoring
            if not Path(f"{args.database}.phr").exists():
                logger.error(f"❌ BLAST database not found: {args.database}")
                logger.error("   Make sure the germline database is built")
                return 1
            show_database_scores(vh_seq, vl_seq, args.database)
            return 0
            
        elif args.graft:
            # CDR grafting only (no optimization)
            logger.info(f"📁 Output directory: {os.path.abspath(output_dir)}")
            
            # Check database
            if not Path(f"{args.database}.phr").exists():
                logger.error(f"❌ BLAST database not found: {args.database}")
                logger.error("   Make sure the germline database is built")
                return 1
            
            backmutation_artifact_dir = (
                prepare_backmutation_artifacts(output_dir)
                if args.backmutation else None
            )
            
            # Perform humanization (grafting only)
            results = humanize_antibody(vh_seq, vl_seq, args.database, backmutation_artifact_dir)
            if args.backmutation:
                write_backmutation_candidates(results, backmutation_artifact_dir)
            
            # Check results
            total_candidates = sum(len(candidates) for candidates in results.values())
            if total_candidates == 0:
                logger.error("❌ No humanization candidates generated")
                logger.error("   This may indicate:")
                logger.error("   1. Input sequences are not compatible with database templates")
                logger.error("   2. CDR boundary detection failed for all candidates")
                logger.error("   3. Database lacks suitable germline sequences")
                return 1
            
            # Save results
            save_results(results, output_dir)
            if args.backmutation:
                run_backmutation(
                    args.input_file,
                    output_dir,
                    vh_seq,
                    vl_seq,
                    args.mouse_database,
                    args.database,
                    backmutation_artifact_dir,
                )
            
            logger.info("🎉 CDR GRAFTING COMPLETE!")
            logger.info("=" * 60)
            logger.info(f"Total humanized antibodies generated: {total_candidates}")
            for chain_type, candidates in results.items():
                if candidates:
                    logger.info(f"{chain_type} Chain: {len(candidates)} candidates")
                    best = max(candidates, key=lambda x: x['identity'])
                    logger.info(f"  Best: {best['identity']:.1f}% identity, {best['source']}")
                    logger.info(f"  Length: {best['length']} AA")
            
            return 0
            
        else:
            # Default: Full optimization pipeline (level 4) - maximum therapeutic quality
            logger.info(f"📁 Output directory: {os.path.abspath(output_dir)}")
            
            # Check database
            if not Path(f"{args.database}.phr").exists():
                logger.error(f"❌ BLAST database not found: {args.database}")
                logger.error("   Make sure the germline database is built")
                return 1
            
            backmutation_artifact_dir = (
                prepare_backmutation_artifacts(output_dir)
                if args.backmutation else None
            )
            
            # Perform humanization
            results = humanize_antibody(vh_seq, vl_seq, args.database, backmutation_artifact_dir)
            backmutation_results = copy.deepcopy(results) if args.backmutation else None
            if args.backmutation:
                write_backmutation_candidates(backmutation_results, backmutation_artifact_dir)
            
            # Determine optimization level: explicit parameter or default to level 4
            optimization_level = args.optimization if args.optimization else 4
            
            logger.info(f"🔧 APPLYING OPTIMIZATION LEVEL {optimization_level}" + 
                       (" (default maximum)" if not args.optimization else " (explicit)"))
            
            # Apply optimization (always runs - either explicit level or default level 4)
            if True:  # Always apply optimization
                vh_candidates = results.get('VH', [])
                vl_candidates = results.get('VL', [])
                
                # Extract CDRs from original sequences for optimization using consistency manager
                consistency_manager = get_consistency_manager()
                vh_result = consistency_manager.extract_regions_consistent(vh_seq, 'heavy')
                vl_result = consistency_manager.extract_regions_consistent(vl_seq, 'light')
                
                # Validate extractions succeeded
                if not (vh_result.extraction_successful and vh_result.validation_passed):
                    logger.error(f"❌ VH CDR extraction failed for optimization: {vh_result.error_message}")
                    return 1
                    
                if not (vl_result.extraction_successful and vl_result.validation_passed):
                    logger.error(f"❌ VL CDR extraction failed for optimization: {vl_result.error_message}")
                    return 1
                
                vh_cdrs = {k: v for k, v in vh_result.regions.items() if k.startswith('CDR')}
                vl_cdrs = {k: v for k, v in vl_result.regions.items() if k.startswith('CDR')}
                
                vh_candidates, vl_candidates = apply_optimization(
                    vh_candidates, vl_candidates, optimization_level,
                    vh_seq, vl_seq,  # Original sequences as murine sequences
                    vh_cdrs, vl_cdrs
                )
                results['VH'] = vh_candidates
                results['VL'] = vl_candidates
            
            # Check results
            total_candidates = sum(len(candidates) for candidates in results.values())
            if total_candidates == 0:
                logger.error("❌ No humanization candidates generated")
                logger.error("   This may indicate:")
                logger.error("   1. Input sequences are not compatible with database templates")
                logger.error("   2. CDR boundary detection failed for all candidates")
                logger.error("   3. Database lacks suitable germline sequences")
                return 1
            
            # Save results
            save_results(results, output_dir)
            if args.backmutation:
                run_backmutation(
                    args.input_file,
                    output_dir,
                    vh_seq,
                    vl_seq,
                    args.mouse_database,
                    args.database,
                    backmutation_artifact_dir,
                )
            
            logger.info("🎉 HUMANIZATION COMPLETE!")
            logger.info("=" * 60)
            logger.info(f"Total humanized antibodies generated: {total_candidates}")
            for chain_type, candidates in results.items():
                if candidates:
                    logger.info(f"{chain_type} Chain: {len(candidates)} candidates")
                    best = max(candidates, key=lambda x: x['identity'])
                    logger.info(f"  Best: {best['identity']:.1f}% identity, {best['source']}")
                    logger.info(f"  Length: {best['length']} AA")
            
            return 0
        
    except SequenceValidationError as e:
        log_error(f"Sequence validation failed: {e}")
        log_error("Please check input sequences for:")
        log_error("  - Minimum length requirements (≥50 AA)")
        log_error("  - Valid amino acid characters only")
        log_error("  - Sufficient amino acid diversity")
        return 1
    except DatabaseError as e:
        log_error(f"Database operation failed: {e}")
        log_error("Please check:")
        log_error("  - BLAST installation and database integrity")
        log_error("  - Database file permissions and accessibility")
        return 1
    except CDRExtractionError as e:
        log_error(f"CDR extraction failed: {e}")
        log_error("Please verify:")
        log_error("  - ANARCII installation and configuration")
        log_error("  - Input sequences are valid antibody domains")
        return 1
    except OptimizationError as e:
        log_error(f"Optimization process failed: {e}")
        log_error("Optimization failure - results may be incomplete")
        return 1
    except FileNotFoundError as e:
        log_error(f"Required file not found: {e}")
        log_error("Please check file paths and permissions")
        return 1
    except Exception as e:
        log_error(f"Unexpected error occurred: {e}")
        log_error("This may indicate a system-level issue")
        logger.exception("Full error details:")
        return 1

if __name__ == "__main__":
    exit(main())
