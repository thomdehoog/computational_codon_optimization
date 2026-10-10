import type { Kind } from "./types";

/** What the page can tell about the pasted text before the server sees it. */
export interface Detection {
  state: "empty" | "ok" | "error";
  message: string;
  /** How the text will be read, once any "auto" setting is resolved. */
  kind: "protein" | "dna" | null;
  /** Residues, or bases for DNA. */
  length: number;
}

const PROTEIN = new Set("ACDEFGHIKLMNPQRSTVWY*");
const DNA = new Set("ACGTU");

/**
 * Describe the pasted sequence. This mirrors the server's reading of the text
 * so the badge matches what will happen, but it only informs: the server stays
 * the authority and reports anything this misses.
 */
export function describe(text: string, requested: Kind, maxResidues: number): Detection {
  const trimmed = text.trim();
  if (!trimmed) return { state: "empty", message: "Nothing pasted yet.", kind: null, length: 0 };

  let body = trimmed;
  if (trimmed.startsWith(">")) {
    const lines = trimmed.split(/\r?\n/);
    const headers = lines.filter((line) => line.startsWith(">")).length;
    if (headers > 1) {
      return error(`${headers} FASTA records found. Paste one sequence at a time.`);
    }
    body = lines.slice(1).join("");
  }
  const sequence = body.replace(/\s+/g, "").toUpperCase();
  if (!sequence) return error("The FASTA header has no sequence under it.");

  const dnaLike = sequence.length >= 9 && [...sequence].every((c) => DNA.has(c));
  const kind = requested === "auto" ? (dnaLike ? "dna" : "protein") : requested;

  if (kind === "dna") {
    const bad = unique([...sequence].filter((c) => !DNA.has(c)));
    if (bad.length) return error(`Not DNA: ${bad.join(", ")}`);
    if (sequence.length % 3) {
      return error(
        `DNA · ${sequence.length} bases is not a whole number of codons (${sequence.length % 3} extra).`,
      );
    }
    const codons = sequence.length / 3;
    return ok("dna", sequence.length, `DNA · ${sequence.length.toLocaleString()} bases, ${codons} codons`);
  }

  const bad = unique([...sequence].filter((c) => !PROTEIN.has(c)));
  if (bad.length) return error(`Not a protein letter: ${bad.join(", ")}`);
  const stop = sequence.indexOf("*");
  if (stop !== -1 && stop !== sequence.length - 1) {
    return error(`Protein · stop symbol (*) at position ${stop + 1}, before the end.`);
  }
  if (sequence.length > maxResidues) {
    return error(
      `Protein · ${sequence.length.toLocaleString()} residues; this page takes up to ${maxResidues.toLocaleString()}.`,
    );
  }
  const noun = sequence.length === 1 ? "residue" : "residues";
  return ok("protein", sequence.length, `Protein · ${sequence.length.toLocaleString()} ${noun}`);
}

function ok(kind: "protein" | "dna", length: number, message: string): Detection {
  return { state: "ok", message, kind, length };
}

function error(message: string): Detection {
  return { state: "error", message, kind: null, length: 0 };
}

function unique(items: string[]): string[] {
  return [...new Set(items)].sort();
}
