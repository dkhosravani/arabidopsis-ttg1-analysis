"""Download accession-verified TTG1 records and UniProt feature annotations.

Sources: NCBI RefSeq/GenBank E-utilities and UniProt REST. Raw source records
are retained; an SHA-256 manifest makes it possible to audit the inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
TRANSCRIPTS = (
    "NM_122360.2",
    "NM_180738.3",
    "NM_180739.2",
    "NM_001343878.1",
)
# Published TTG1/WD40 homologues; their *homology* does not independently prove
# a one-to-one orthology relationship. See the references in README.md.
HOMOLOGS = {
    "Brassica_rapa": "ABQ10570",
    "Prunus_persica": "ACQ65867",
    "Solanum_lycopersicum": "XP_004235332",
    "Oryza_sativa": "BAF09665",
}
REFERENCE_UNIPROT = "Q9XGN1"


def download(url: str) -> bytes:
    """Fetch a small public source record with bounded retries and timeout."""
    last_error = None
    for attempt in range(3):
        try:
            request = Request(url, headers={"User-Agent": "ttg1-comparative-sequence-study/1.0"})
            with urlopen(request, timeout=45) as response:
                return response.read()
        except (URLError, HTTPError, TimeoutError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Could not retrieve {url}: {last_error}")


def ncbi_record(accession: str, database: str, format_name: str):
    """Fetch and validate precisely one versioned NCBI record."""
    params = urlencode({
        "db": database,
        "id": accession,
        "rettype": format_name,
        "retmode": "text",
        "tool": "TTG1-sequence-study",
    })
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{params}"
    payload = download(url)
    text = payload.decode("utf-8")
    file_format = "genbank" if format_name == "gb" else "fasta"
    records = list(SeqIO.parse(io.StringIO(text), file_format))
    if len(records) != 1:
        raise ValueError(f"Expected one record for {accession}; found {len(records)}")
    actual = records[0].id
    if "." in accession and actual != accession:
        raise ValueError(f"Accession/version mismatch: requested {accession}; got {actual}")
    if "." not in accession and actual.split(".")[0] != accession:
        raise ValueError(f"Accession mismatch: requested {accession}; got {actual}")
    if not records[0].seq:
        raise ValueError(f"Empty sequence for {accession}")
    return payload, records[0], url


def save_bytes(path: Path, content: bytes, refresh: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not refresh:
        if path.read_bytes() != content:
            raise FileExistsError(f"{path} already exists and differs. Use --refresh.")
    else:
        path.write_bytes(content)


def main(refresh: bool = False) -> None:
    entries = []
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # Use the GenBank-format transcript record to retain the *annotated* CDS,
    # rather than assuming the first ATG in the mRNA is the translation start.
    for accession in TRANSCRIPTS:
        content, record, url = ncbi_record(accession, "nuccore", "gb")
        path = RAW / "transcripts" / f"{accession}.gb"
        save_bytes(path, content, refresh)
        entries.append((path, accession, "NCBI RefSeq mRNA", url))
        print(f"Downloaded mRNA: {record.id} ({len(record.seq)} nt)")
        time.sleep(0.35)  # NCBI's no-API-key rate limit is 3 requests/second.

    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{REFERENCE_UNIPROT}.json"
    payload = download(uniprot_url)
    uniprot = json.loads(payload)
    if uniprot.get("primaryAccession") != REFERENCE_UNIPROT:
        raise ValueError("UniProt accession mismatch")
    if len(uniprot["sequence"]["value"]) != uniprot["sequence"]["length"]:
        raise ValueError("UniProt sequence length mismatch")
    path = RAW / f"uniprot_{REFERENCE_UNIPROT}.json"
    save_bytes(path, payload, refresh)
    entries.append((path, REFERENCE_UNIPROT, "UniProtKB protein and features", uniprot_url))
    print(f"Downloaded reference protein: {REFERENCE_UNIPROT} ({uniprot['sequence']['length']} aa)")

    for species, accession in HOMOLOGS.items():
        content, record, url = ncbi_record(accession, "protein", "fasta")
        path = RAW / "proteins" / f"{species}_{record.id}.fasta"
        save_bytes(path, content, refresh)
        entries.append((path, record.id, f"NCBI protein, {species}", url))
        print(f"Downloaded homolog: {species} / {record.id} ({len(record.seq)} aa)")
        time.sleep(0.35)

    manifest = ROOT / "data" / "source_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "accession", "source", "retrieved_utc", "sha256", "source_url"])
        for path, accession, source, url in entries:
            writer.writerow([
                path.relative_to(ROOT).as_posix(), accession, source, timestamp,
                hashlib.sha256(path.read_bytes()).hexdigest(), url,
            ])
    print(f"Source manifest: {manifest.relative_to(ROOT)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Replace previously downloaded source records")
    args = parser.parse_args()
    main(refresh=args.refresh)
