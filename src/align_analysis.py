"""Analyze TTG1 mRNA, its annotated coding sequence, and protein homologs."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

from .fetch_seq import HOMOLOGS, RAW, REFERENCE_UNIPROT, ROOT, TRANSCRIPTS

OUTPUT = ROOT / "results"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def gc_percent(sequence: str | Seq) -> float:
    """Calculate percent G+C using Biopython; undefined positions are ignored."""
    if not sequence:
        raise ValueError("Cannot calculate GC% of an empty sequence")
    return round(100 * gc_fraction(sequence, ambiguous="remove"), 3)


def sliding_gc(sequence: str | Seq, window: int = 100, step: int = 10):
    """Return the midpoint (1-based) and GC% of each complete sliding window."""
    if window <= 0 or step <= 0:
        raise ValueError("Window and step must be positive")
    sequence = str(sequence).upper().replace("U", "T")
    if len(sequence) < window:
        return []
    return [
        {"position_nt": start + 1 + (window - 1) / 2,
         "gc_percent": gc_percent(sequence[start:start + window])}
        for start in range(0, len(sequence) - window + 1, step)
    ]


def scan_orfs(sequence: str | Seq, min_aa: int = 30):
    """Find complete ATG-to-stop ORFs on the transcribed (+) mRNA strand.

    These are candidate ORFs, not experimentally validated coding regions.
    Nested ATGs are retained; positions are 1-based inclusive on the mRNA.
    """
    sequence = str(sequence).upper().replace("U", "T")
    result = []
    for start in range(len(sequence) - 2):
        if sequence[start:start + 3] != "ATG":
            continue
        for stop_start in range(start + 3, len(sequence) - 2, 3):
            if sequence[stop_start:stop_start + 3] in STOP_CODONS:
                aa_length = (stop_start - start) // 3
                if aa_length >= min_aa:
                    result.append({
                        "start_nt": start + 1,
                        "end_nt": stop_start + 3,
                        "frame": start % 3 + 1,
                        "length_aa_excluding_stop": aa_length,
                        "stop_codon": sequence[stop_start:stop_start + 3],
                    })
                break
    return sorted(result, key=lambda item: (-item["length_aa_excluding_stop"], item["start_nt"]))


def annotated_cds(record):
    """Extract the deposited CDS feature and translation from a GenBank mRNA."""
    features = [feature for feature in record.features if feature.type == "CDS"]
    if len(features) != 1:
        raise ValueError(f"Expected one CDS in {record.id}; found {len(features)}")
    cds = features[0]
    nucleotide = cds.extract(record.seq)
    translated = str(nucleotide.translate(to_stop=True))
    expected = cds.qualifiers.get("translation", [None])[0]
    if expected and translated != expected:
        raise ValueError(f"Translation mismatch in NCBI annotated CDS of {record.id}")
    return nucleotide, translated, cds


def protein_aligner():
    """Create a reproducible end-to-end global alignment with BLOSUM62."""
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def identity_percent(alignment, include_gaps: bool = True) -> float:
    """Calculate the percentage of identical amino acids in an alignment.

    If include_gaps is True, all alignment columns are counted.
    Otherwise, columns containing a gap in either sequence are excluded.

    This measures exact identity, not BLOSUM similarity.
    """
    indices = alignment.indices
    matches = sum(
        a >= 0 and b >= 0 and alignment.target[a] == alignment.query[b]
        for a, b in zip(indices[0], indices[1])
    )
    aligned_columns = sum(
        a >= 0 and b >= 0
        for a, b in zip(indices[0], indices[1])
    )

    denominator = indices.shape[1] if include_gaps else aligned_columns

    return round(100.0 * matches / denominator, 2) if denominator else 0.0

def identity_in_reference_region(alignment, start_aa: int, end_aa: int) -> float:
    """Calculate identity for a region of the Arabidopsis reference protein."""
    if not (1 <= start_aa <= end_aa <= len(alignment.target)):
        raise ValueError("Invalid reference region coordinates")

    matches = 0

    for ref_pos, other_pos in zip(*alignment.indices):
        if start_aa - 1 <= ref_pos < end_aa:
            if other_pos >= 0 and alignment.target[ref_pos] == alignment.query[other_pos]:
                matches += 1

    region_length = end_aa - start_aa + 1
    return round(100.0 * matches / region_length, 2)

def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mrna_rows, gc_windows, orf_rows = [], [], []
    proteins = {}
    uniprot_path = RAW / f"uniprot_{REFERENCE_UNIPROT}.json"
    if not uniprot_path.exists():
        raise FileNotFoundError("Reference files missing. Run: python -m src.fetch_seq")
    uniprot = json.loads(uniprot_path.read_text(encoding="utf-8"))
    reference_protein = uniprot["sequence"]["value"]
    proteins["Arabidopsis_thaliana"] = reference_protein
    wd40 = []
    for feature in uniprot.get("features", []):
        if feature.get("type") == "Repeat" and "WD" in feature.get("description", "").upper():
            wd40.append({
                "label": feature.get("description", "WD40 repeat"),
                "start_aa": int(feature["location"]["start"]["value"]),
                "end_aa": int(feature["location"]["end"]["value"]),
                "source": REFERENCE_UNIPROT,
            })
    if not wd40:
        raise ValueError("No WD repeat annotations found in the UniProt record")
    (OUTPUT / "wd40_annotations.json").write_text(json.dumps(wd40, indent=2) + "\n")

    reference_translated = None
    for accession in TRANSCRIPTS:
        gb_path = RAW / "transcripts" / f"{accession}.gb"
        record = SeqIO.read(gb_path, "genbank")
        mrna = str(record.seq).upper()
        cds_seq, peptide, cds_feature = annotated_cds(record)
        if accession == "NM_122360.2":
            reference_translated = peptide
        mrna_rows.append({
            "accession": accession,
            "mrna_length_nt": len(mrna),
            "mrna_gc_percent": gc_percent(mrna),
            "annotated_cds_length_nt": len(cds_seq),
            "cds_gc_percent": gc_percent(cds_seq),
            "annotated_protein_length_aa": len(peptide),
            "protein_exact_match_to_UniProt_Q9XGN1": peptide == reference_protein,
        })
        for point in sliding_gc(mrna):
            gc_windows.append({"accession": accession, **point})
        for orf in scan_orfs(mrna):
            orf_rows.append({"accession": accession, **orf})

    write_csv(OUTPUT / "transcript_summary.csv", mrna_rows, list(mrna_rows[0]))
    write_csv(OUTPUT / "sliding_gc.csv", gc_windows,
              ["accession", "position_nt", "gc_percent"])
    write_csv(OUTPUT / "orf_candidates.csv", orf_rows,
              ["accession", "start_nt", "end_nt", "frame", "length_aa_excluding_stop", "stop_codon"])

    for species, accession in HOMOLOGS.items():
        matching = list((RAW / "proteins").glob(f"{species}_{accession}*.fasta"))
        if len(matching) != 1:
            raise ValueError(f"Expected one protein FASTA for {species}, found {len(matching)}")
        record = SeqIO.read(matching[0], "fasta")
        protein = str(record.seq).replace("*", "").upper()
        if not protein or set(protein) - set("ACDEFGHIKLMNPQRSTVWYUXOBZJ"):
            raise ValueError(f"Invalid amino acid sequence from {record.id}")
        proteins[species] = protein

    # Include complete input protein sequences to make this analysis easy to
    # reproduce offline without fetching again.
    SeqIO.write(
        [SeqRecord(Seq(seq), id=name, description="TTG1 homologue") for name, seq in proteins.items()],
        OUTPUT / "homologs_used.fasta", "fasta",
    )

    aligner = protein_aligner()
    names = list(proteins)
    matrix = []
    with (OUTPUT / "pairwise_alignments.txt").open("w", encoding="utf-8") as handle:
        for name in names[1:]:
            alignment = aligner.align(proteins[names[0]], proteins[name])[0]
            handle.write(f"{names[0]} vs {name}\n")
            handle.write(f"Global BLOSUM62; open gap -10, extend gap -0.5\n")
            handle.write(f"Identity (including gaps): {identity_percent(alignment):.2f}%\n")
            handle.write(f"Identity (excluding gaps): {identity_percent(alignment, include_gaps=False):.2f}%\n\n")
            handle.write(str(alignment) + "\n\n")
    for first in names:
        row = {"species": first}
        for second in names:
            row[second] = 100.0 if first == second else identity_percent(
                aligner.align(proteins[first], proteins[second])[0]
            )
        matrix.append(row)
    write_csv(OUTPUT / "pairwise_identity_matrix.csv", matrix, ["species"] + names)

    wd40_rows = []

    for species in names[1:]:
        alignment = aligner.align(reference_protein, proteins[species])[0]

        full_identity = identity_in_reference_region(
            alignment, 1, len(reference_protein)
        )

        for number, repeat in enumerate(wd40, start=1):
            wd40_rows.append({
                "species": species,
                "region": f"WD{number}",
                "start_aa": repeat["start_aa"],
                "end_aa": repeat["end_aa"],
                "wd_identity_percent": identity_in_reference_region(
                    alignment, repeat["start_aa"], repeat["end_aa"]
                ),
                "full_reference_identity_percent": full_identity,
            })

    write_csv(
        OUTPUT / "wd40_identity.csv",
        wd40_rows,
        [
            "species",
            "region",
            "start_aa",
            "end_aa",
            "wd_identity_percent",
            "full_reference_identity_percent",
        ],
    )
    summary = {
        "gene": "TTG1 / AT5G24520",
        "reference_mrna": "NM_122360.2",
        "reference_protein": REFERENCE_UNIPROT,
        "reference_protein_length_aa": len(reference_protein),
        "reference_annotated_CDS_matches_UniProt": reference_translated == reference_protein,
        "transcript_count": len(TRANSCRIPTS),
        "homolog_species_count_including_reference": len(proteins),
        "wd40_repeat_count": len(wd40),
        "pairwise_alignment_method": "global BLOSUM62; gap open -10; gap extend -0.5",
        "identity_definition": "identical residues divided by all alignment columns (including gaps)",
        "scope_note": "Homology supported by source literature; one-to-one orthology and species phylogeny not inferred.",
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
