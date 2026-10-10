"""Minimal FASTA reading and writing.

Keeping this in-house avoids a dependency on Biopython for what amounts to
twenty lines, and it keeps the CLI usable in a bare environment.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

__all__ = ["Record", "parse_fasta", "read_fasta", "write_fasta", "format_fasta"]


@dataclass(frozen=True)
class Record:
    """One FASTA record.

    Attributes:
        name: The identifier, i.e. the header up to the first space.
        description: The full header line without the leading ``>``.
        sequence: The sequence with whitespace removed.
    """

    name: str
    description: str
    sequence: str


def parse_fasta(lines: Iterable[str]) -> Iterator[Record]:
    """Yield records from FASTA-formatted ``lines``.

    A file with no header line is treated as a single unnamed sequence, which
    makes it easy to pass a bare sequence on the command line.

    Raises:
        ValueError: If a record has a header but no sequence.
    """
    header: str | None = None
    chunks: list[str] = []
    saw_header = False

    for line in lines:
        line = line.strip()
        if not line or line.startswith((";", "#")):
            continue
        if line.startswith(">"):
            if header is not None or chunks:
                yield _build(header, chunks)
            header = line[1:].strip()
            saw_header = True
            chunks = []
        else:
            chunks.append(line.replace(" ", ""))
    if header is not None or chunks:
        yield _build(header, chunks)
    elif saw_header:  # pragma: no cover - defensive
        raise ValueError("FASTA input ended after a header with no sequence")


def _build(header: str | None, chunks: list[str]) -> Record:
    sequence = "".join(chunks)
    if not sequence:
        raise ValueError(
            f"FASTA record {header or '(unnamed)'!r} has a header but no sequence"
        )
    description = header or ""
    name = description.split()[0] if description else "sequence"
    return Record(name=name, description=description, sequence=sequence)


def read_fasta(path: str | Path) -> list[Record]:
    """Read every record from ``path``.

    Raises:
        ValueError: If the file holds no records.
    """
    path = Path(path)
    with path.open() as handle:
        records = list(parse_fasta(handle))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def format_fasta(record: Record, *, width: int = 60) -> str:
    """Return ``record`` as FASTA text, wrapped at ``width`` columns."""
    if width < 1:
        raise ValueError("width must be at least 1")
    body = "\n".join(
        record.sequence[i : i + width] for i in range(0, len(record.sequence), width)
    )
    header = record.description or record.name
    return f">{header}\n{body}\n"


def write_fasta(
    path: str | Path, records: Iterable[Record], *, width: int = 60
) -> None:
    """Write ``records`` to ``path`` in FASTA format."""
    text = "".join(format_fasta(record, width=width) for record in records)
    Path(path).write_text(text)
