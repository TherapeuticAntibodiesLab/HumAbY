"""Add reference identifiers to auxiliary reports, without changing the scientific files.

Identifiers come from saved FASTA/XML. Gene descriptions are retrieved by ID
from the existing local database, never by running a new sequence alignment.
Origins are verified against preserved reference catalogs, not inferred from
the current setting or a folder named ncbi-igblast.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import subprocess
from functools import lru_cache
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIXML

from integration.paths import PROJECT, CATALOGS, resolve_path
HUMAN_IMGT_FASTA_SHA256 = "147e07c50634c0d7054d1d75bd2d3cf6ffd4f653584ae9105d5e801bd4908c5f"


@lru_cache(maxsize=64)
def file_digest_cached(path: str, mtime: int, size: int) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def file_digest(path: Path) -> str:
    info = path.stat()
    return file_digest_cached(str(path.resolve()), info.st_mtime_ns, info.st_size)


@lru_cache(maxsize=256)
def catalog_entry_cached(prefix: str, identifier: str, signatures: tuple) -> dict | None:
    try:
        result = subprocess.run(["blastdbcmd", "-db", prefix, "-entry", identifier,
                                 "-outfmt", "%f"], capture_output=True, text=True,
                                check=True, timeout=10)
        records = list(SeqIO.parse(io.StringIO(result.stdout), "fasta"))
        if len(records) != 1:
            return None
        record = records[0]
        return {"id": record.id, "description": record.description, "sequence": str(record.seq)}
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def catalog_entry(prefix: Path, identifier: str) -> dict | None:
    files = [Path(str(prefix) + extension) for extension in (".psq", ".pin", ".phr")]
    if not all(path.is_file() for path in files):
        return None
    signatures = tuple((path.stat().st_mtime_ns, path.stat().st_size) for path in files)
    return catalog_entry_cached(str(prefix.resolve()), identifier, signatures)


def _mouse_catalogs() -> list:
    provenance = CATALOGS
    if not provenance.is_file():
        return []
    data = json.loads(provenance.read_text(encoding="utf-8"))
    catalogs = []
    for key, origin in (("active", "NCBI"), ("preserved_inactive", "IMGT")):
        reference = data[key]
        prefix = PROJECT / reference["database_prefix"]
        if not all((prefix.parent / filename).is_file()
                   and file_digest(prefix.parent / filename) == digest
                   for filename, digest in reference["sha256"].items()):
            continue
        catalogs.append({"origin": origin, "prefix": prefix,
                         "sequences": reference["sequences"], "residues": reference["residues"]})
    return catalogs


def _human_reference(candidate: dict, prefix: Path | None) -> dict:
    saved = candidate["proposal_id"]
    match = re.search(r"\|ref\|((?:VH|VK|VL)_\d+)\|", saved)
    if not match:
        match = re.fullmatch(r"(?:ref\|)?((?:VH|VK|VL)_\d+)\|?", saved)
    identifier = match[1] if match else None
    reference = {"origin": None, "identifier": identifier,
                 "identifier_type": "Internal reconstructed-framework ID; not a GenBank accession",
                 "saved_fasta_id": saved, "genes": {}, "description": None,
                 "sequence_verified": False,
                 "note": "Gene description unavailable; no new alignment was performed."}
    if identifier is None or prefix is None:
        return reference
    entry = catalog_entry(prefix, identifier)
    if entry is None:
        return reference
    if entry["sequence"] != candidate["human_germline"]:
        reference["note"] = "Local database entry does not match saved hGerm; gene names and origin were not assigned."
        return reference
    reference.update(sequence_verified=True, description=entry["description"],
                     database_path=str(prefix),
                     genes=dict(re.findall(r"\b([VDJ]):([^\s|]+)", entry["description"])),
                     note="Description retrieved by saved ID from the local catalog; exact match to saved hGerm. This is the grafting template, not a new annotation of backmt.")
    fasta = Path(str(prefix) + ".fasta")
    if fasta.is_file():
        reference["database_fasta_sha256"] = file_digest(fasta)
        if reference["database_fasta_sha256"] == HUMAN_IMGT_FASTA_SHA256:
            reference["origin"] = "IMGT"
    return reference


def _mouse_reference(xml: Path, chain: dict, catalogs: list) -> dict:
    reference = {"origin": None, "identifier": None, "accession": None,
                 "identifier_type": "Reference record code; not necessarily a GenBank accession",
                 "description": None, "xml_path": None, "database_path_recorded_in_xml": None,
                 "sequence_verified": False,
                 "note": "Reference identifier unavailable: no saved BLAST XML."}
    if not xml.is_file():
        return reference
    with xml.open(encoding="utf-8") as handle:
        record = NCBIXML.read(handle)
    if len(record.alignments) != 1 or len(record.alignments[0].hsps) != 1:
        raise ValueError(f"Ambiguous reference alignment in {xml.name}.")
    hit = record.alignments[0]
    hsp = hit.hsps[0]
    if chain["mouse_germline"] not in (None, hsp.sbjct):
        raise ValueError(f"{xml.name} and raw report contain different murine alignments.")
    if chain["input_sequence"] is not None:
        saved_query = chain["input_sequence"][hsp.query_start - 1:hsp.query_end]
        if saved_query != hsp.query.replace("-", ""):
            raise ValueError(f"{xml.name} and raw report belong to different input sequences.")
    reference.update(identifier=hit.hit_id, accession=hit.accession,
                     description=hit.hit_def, xml_path="output/" + xml.name,
                     database_path_recorded_in_xml=record.database,
                     saved_database_sequences=record.database_sequences,
                     saved_database_residues=record.database_length,
                     note="Identifier and aligned sequence read from the original BLAST XML. Database origin not verified.")
    verified = []
    for catalog in catalogs:
        if (record.database_sequences != catalog["sequences"]
                or record.database_length != catalog["residues"]):
            continue
        entry = catalog_entry(catalog["prefix"], hit.hit_id)
        if entry and entry["sequence"][hsp.sbjct_start - 1:hsp.sbjct_end] == hsp.sbjct.replace("-", ""):
            verified.append(catalog)
    if len(verified) == 1:
        reference.update(origin=verified[0]["origin"], sequence_verified=True,
                         note="Saved hit ID, aligned residues and database statistics match the checksum-verified preserved catalog. This identifies the reference; it does not prove an entire historical database version.")
    return reference


def add_reference_info(document: dict, directory: Path) -> dict:
    catalogs = _mouse_catalogs()
    human_path = document["run"]["database_path"]
    prefix = resolve_path(human_path) if human_path else None
    for name, chain in document["chains"].items():
        xml = directory / (name + ".xml")
        chain["mouse_reference"] = _mouse_reference(xml, chain, catalogs)
        if xml.is_file():
            document["sources"].append({"path": "output/" + xml.name, "sha256": file_digest(xml)})
        for candidate in chain["candidates"]:
            candidate["human_reference"] = _human_reference(candidate, prefix)
    document["reference_schema_version"] = 1
    return document
