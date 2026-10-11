/** Shapes of the JSON the codonopt server sends. Mirrors codonopt/server.py. */

export interface Defaults {
  gc_target: number;
  gc_range: [number, number];
  gc_window: number;
  gc_step: number;
  repeat_k: number;
  hairpin_stem: number;
  hairpin_min_loop: number;
  hairpin_max_loop: number;
  max_homopolymer: number;
  seconds: number;
}

export interface EffortLevel {
  max_sweeps: number;
  restarts: number;
  anneal_steps: number;
}

export interface Meta {
  version: string;
  tables: { name: string; description: string }[];
  default_table: string;
  motif_groups: Record<string, string[]>;
  restriction_sites: Record<string, string>;
  regulatory_motifs: Record<string, string>;
  splice_backends: string[];
  maxent_available: boolean;
  weights: Record<string, number>;
  defaults: Defaults;
  effort: Record<string, EffortLevel>;
  limits: { max_residues: number; max_seconds: number };
  example: { name: string; protein: string };
}

export interface Measures {
  context_length: number;
  cds_start: number;
  length: number;
  gc_fraction: number;
  gc_windows_below_range: number;
  gc_windows_above_range: number;
  cai: number;
  repeat_k: number;
  redundant_kmers: number;
  hairpin_pairs: number;
  longest_hairpin_stem: number;
  longest_homopolymer: number;
  homopolymer_excess_bases: number;
  forbidden_motif_hits: number;
  forbidden_motifs_present: string[];
  splice_backend?: string;
  splice_donor_sites?: number;
  splice_acceptor_sites?: number;
}

export interface GcProfile {
  window: number;
  step: number;
  range: [number, number];
  target: number;
  starts: number[];
  initial: number[];
  final: number[];
}

export interface MotifHit {
  name: string;
  length: number;
  positions: number[];
}

export interface Result {
  context_sequence: string;
  valid: boolean;
  unresolved_constraints: { criterion: string; count: number }[];
  stop_reason: string;
  name: string;
  kind: "protein" | "dna";
  protein: string;
  sequence: string;
  initial_sequence: string;
  cancelled: boolean;
  codon_usage_table: string;
  score: {
    initial: number;
    final: number;
    improvement: number;
    initial_terms: Record<string, number>;
    final_terms: Record<string, number>;
  };
  search: {
    restarts_run: number;
    greedy_sweeps: number;
    moves_evaluated: number;
    moves_accepted: number;
    seconds: number;
    hit_time_limit: boolean;
  };
  initial: Measures;
  final: Measures;
  gc: GcProfile;
  motif_hits: MotifHit[];
  fasta: string;
}

export type JobStatus = "running" | "done" | "cancelled" | "error";

export interface Job {
  id: string;
  status: JobStatus;
  elapsed: number;
  max_seconds: number;
  restarts_total: number;
  restarts_done: number;
  messages: string[];
  error?: string;
  result?: Result;
}

export type Kind = "auto" | "protein" | "dna";
