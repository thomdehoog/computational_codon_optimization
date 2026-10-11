import { gcChart } from "./chart";
import { fmt, h, int, replace } from "./dom";
import type { Child } from "./dom";
import { sequenceView } from "./sequence";
import type { Job, Measures, Meta, Result } from "./types";

/* ---------- small helpers ---------- */

async function copyText(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    // Older embedded web views have no async clipboard: fall back to a hidden textarea.
    const area = h("textarea", { attrs: { "aria-hidden": "true", tabindex: "-1" } });
    area.value = text;
    document.body.append(area);
    area.select();
    let done = false;
    try {
      done = document.execCommand("copy");
    } catch {
      done = false;
    }
    area.remove();
    return done;
  }
}

function download(filename: string, text: string, type: string): void {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = h("a", { attrs: { href: url, download: filename } });
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 10_000);
}

function action(label: string, onClick: (button: HTMLButtonElement) => void, primary = false): HTMLButtonElement {
  const button: HTMLButtonElement = h("button", {
    type: "button",
    class: primary ? "chip primary-chip" : "chip",
    text: label,
    on: { click: () => onClick(button) },
  });
  return button;
}

const flashing = new WeakMap<HTMLButtonElement, number>();

/** Briefly replace a button's label, surviving repeated clicks. */
function flash(button: HTMLButtonElement, message: string): void {
  const pending = flashing.get(button);
  if (pending !== undefined) window.clearTimeout(pending);
  // Capture the real label once; a second click must not capture "Copied".
  if (button.dataset.label === undefined) button.dataset.label = button.textContent ?? "";
  button.textContent = message;
  flashing.set(
    button,
    window.setTimeout(() => {
      button.textContent = button.dataset.label ?? "";
      flashing.delete(button);
    }, 1600),
  );
}

const safeName = (name: string): string => name.replace(/[^\w.-]+/g, "_") || "sequence";

/* ---------- empty, running and error states ---------- */

export function showError(output: HTMLElement, message: string): void {
  replace(
    output,
    h(
      "div",
      { class: "callout error", role: "alert" },
      h("strong", { text: "That did not work" }),
      h("p", { text: message }),
    ),
  );
  output.focus({ preventScroll: true });
}

/** A progress card that is built once and updated in place on every poll. */
export class Progress {
  private readonly bar = h("div", { class: "bar-fill" });
  private readonly elapsed = h("span", { class: "mono" });
  private readonly detail = h("span", { class: "hint" });
  private readonly log = h("ol", { class: "log mono", attrs: { "aria-label": "Search log" } });
  private readonly trouble = h("p", { class: "note warn", role: "alert" });
  private lastLog = "";

  constructor(output: HTMLElement) {
    this.trouble.hidden = true;
    replace(
      output,
      h(
        "div",
        { class: "card progress" },
        h("h2", { text: "Optimizing" }),
        h("p", { class: "hint", text: "Resolving mandatory requirements, then improving synonymous codon choices." }),
        h("div", { class: "bar", role: "progressbar", attrs: { "aria-label": "Optimization progress" } }, this.bar),
        h("div", { class: "progress-meta" }, this.elapsed, this.detail),
        this.trouble,
        this.log,
      ),
    );
    output.focus({ preventScroll: true });
  }

  /** Report something that went wrong without throwing the run's state away. */
  note(message: string): void {
    this.trouble.textContent = message;
    this.trouble.hidden = false;
  }

  update(job: Job): void {
    this.elapsed.textContent = `${fmt(job.elapsed, 1)} s of at most ${int(job.max_seconds)} s`;
    const track = this.bar.parentElement;
    if (job.restarts_total > 1) {
      const done = Math.min(job.restarts_done, job.restarts_total);
      this.detail.textContent = `Restart ${Math.min(done + 1, job.restarts_total)} of ${job.restarts_total}`;
      this.bar.style.width = `${(done / job.restarts_total) * 100}%`;
      track?.classList.remove("indeterminate");
      track?.setAttribute("aria-valuemin", "0");
      track?.setAttribute("aria-valuemax", String(job.restarts_total));
      track?.setAttribute("aria-valuenow", String(done));
    } else {
      this.detail.textContent = "Stops as soon as nothing improves.";
      track?.classList.add("indeterminate");
      track?.removeAttribute("aria-valuenow");
    }
    const key = job.messages.join("\n");
    if (key !== this.lastLog) {
      this.lastLog = key;
      replace(this.log, ...job.messages.map((message) => h("li", { text: message })));
    }
  }
}

/* ---------- the result ---------- */

interface Row {
  label: string;
  help: string;
  get: (m: Measures) => number;
  format?: (value: number) => string;
  /** How to judge the change: fewer is better, more is better, or nearer the target. */
  better: "lower" | "higher" | "target";
  show?: (m: Measures) => boolean;
}

const percent = (v: number): string => `${fmt(v * 100, 1)}%`;

function rows(after: Measures): Row[] {
  const hasSplice = after.splice_backend !== undefined;
  return [
    {
      label: "GC content",
      help: "Share of G and C bases over the whole sequence.",
      get: (m) => m.gc_fraction,
      format: percent,
      better: "target",
    },
    {
      label: "GC windows out of range",
      help: "Windows whose local GC falls outside the allowed range.",
      get: (m) => m.gc_windows_below_range + m.gc_windows_above_range,
      better: "lower",
    },
    {
      label: "Codon adaptation index",
      help: "How closely the codons follow the host's preferences (1 is ideal).",
      get: (m) => m.cai,
      format: (v) => fmt(v, 3),
      better: "higher",
    },
    {
      label: `Repeated ${after.repeat_k}-mers`,
      help: "Extra copies of any stretch that appears more than once. Repeats misassemble during synthesis.",
      get: (m) => m.redundant_kmers,
      better: "lower",
    },
    {
      label: "Hairpin pairs",
      help: "Pairs of inverted repeats that can fold into stem-loops.",
      get: (m) => m.hairpin_pairs,
      better: "lower",
    },
    {
      label: "Longest single-base run",
      help: "The longest stretch of one repeated base.",
      get: (m) => m.longest_homopolymer,
      better: "lower",
    },
    {
      label: "Forbidden sites",
      help: "Matches to the enzymes and patterns you asked to keep out.",
      get: (m) => m.forbidden_motif_hits,
      better: "lower",
    },
    {
      label: "Splice donor sites",
      help: "Windows that look like a 5′ splice site.",
      get: (m) => m.splice_donor_sites ?? 0,
      better: "lower",
      show: () => hasSplice,
    },
    {
      label: "Splice acceptor sites",
      help: "Windows that look like a 3′ splice site.",
      get: (m) => m.splice_acceptor_sites ?? 0,
      better: "lower",
      show: () => hasSplice,
    },
  ];
}

/** Everything the optimizer was asked to remove, as one headline number. */
function issues(m: Measures): number {
  return (
    m.gc_windows_below_range +
    m.gc_windows_above_range +
    m.redundant_kmers +
    m.hairpin_pairs +
    m.homopolymer_excess_bases +
    m.forbidden_motif_hits +
    (m.splice_donor_sites ?? 0) +
    (m.splice_acceptor_sites ?? 0)
  );
}

function comparison(result: Result): HTMLElement {
  const target = result.gc.target;
  const body = rows(result.final)
    .filter((row) => !row.show || row.show(result.final))
    .map((row) => {
      const before = row.get(result.initial);
      const after = row.get(result.final);
      const format = row.format ?? int;
      let verdict: "good" | "bad" | "same";
      if (row.better === "target") {
        const distance = Math.abs(after - target) - Math.abs(before - target);
        verdict = Math.abs(distance) < 5e-4 ? "same" : distance < 0 ? "good" : "bad";
      } else {
        const delta = after - before;
        verdict = Math.abs(delta) < 1e-9 ? "same" : (delta < 0) === (row.better === "lower") ? "good" : "bad";
      }
      const words = { good: "better", bad: "worse", same: "no change" }[verdict];
      const mark = { good: "✓", bad: "!", same: "–" }[verdict];
      return h(
        "tr",
        null,
        h("th", { attrs: { scope: "row", title: row.help } }, row.label),
        h("td", { class: "num", text: format(before) }),
        h("td", { class: "num strong", text: format(after) }),
        h("td", { class: `verdict ${verdict}` }, h("span", { attrs: { "aria-hidden": "true" } }, `${mark} `), words),
      );
    });
  return h(
    "table",
    { class: "compare" },
    h(
      "thead",
      null,
      h(
        "tr",
        null,
        h("th", { attrs: { scope: "col" } }, "Property"),
        h("th", { class: "num", attrs: { scope: "col" } }, "Before"),
        h("th", { class: "num", attrs: { scope: "col" } }, "After"),
        h("th", { attrs: { scope: "col" } }, "Change"),
      ),
    ),
    h("tbody", null, ...body),
  );
}

const TERMS: Record<string, [string, string]> = {
  gc_global: ["Overall GC", "Distance of the whole-sequence GC from the target."],
  gc_window: ["Local GC", "How far windows stray outside the allowed range."],
  repeat: ["Repeats", "Redundant copies of repeated stretches."],
  hairpin: ["Hairpins", "Inverted repeats that can fold back on themselves."],
  homopolymer: ["Single-base runs", "Bases beyond the longest run you allow."],
  motif: ["Forbidden sites", "Matches to enzymes and patterns you keep out."],
  splice: ["Splice sites", "How far splice-site scores exceed their cutoffs."],
  usage: ["Codon usage", "Cost of using codons the host dislikes."],
};

function breakdown(result: Result): HTMLElement {
  const before = result.score.initial_terms;
  const after = result.score.final_terms;
  const names = Object.keys(TERMS).filter((key) => (before[key] ?? 0) > 1e-9 || (after[key] ?? 0) > 1e-9);
  if (!names.length) {
    return h("p", { class: "hint", text: "No active soft penalties. Mandatory requirements are checked separately." });
  }
  const scale = Math.max(...names.flatMap((key) => [before[key] ?? 0, after[key] ?? 0]), 1e-9);
  const bar = (value: number, className: string): HTMLElement => {
    const fill = h("span", { class: `fill ${className}` });
    fill.style.width = `${Math.max(value > 0 ? 1.5 : 0, (value / scale) * 100)}%`;
    return h("span", { class: "track" }, fill);
  };
  return h(
    "div",
    { class: "terms" },
    ...names.map((key) =>
      h(
        "div",
        { class: "term", title: TERMS[key][1] },
        h("span", { class: "term-name", text: TERMS[key][0] }),
        h(
          "span",
          { class: "term-bars" },
          bar(before[key] ?? 0, "before"),
          bar(after[key] ?? 0, "after"),
        ),
        h("span", { class: "term-values mono", text: `${fmt(before[key] ?? 0, 2)} → ${fmt(after[key] ?? 0, 2)}` }),
      ),
    ),
  );
}

function card(title: string, subtitle: string | null, ...content: Child[]): HTMLElement {
  return h(
    "section",
    { class: "card" },
    h("h3", { text: title }),
    subtitle ? h("p", { class: "hint", text: subtitle }) : null,
    ...content,
  );
}

function stat(label: string, before: string, after: string, note: string): HTMLElement {
  return h(
    "div",
    { class: "stat" },
    h("span", { class: "stat-label", text: label }),
    h("span", { class: "stat-value" }, h("span", { class: "was", text: before }), " → ", h("strong", { text: after })),
    h("span", { class: "stat-note", text: note }),
  );
}

export function showResult(output: HTMLElement, result: Result, _meta: Meta): void {
  const issuesBefore = issues(result.initial);
  const issuesAfter = issues(result.final);
  const reduction =
    result.score.initial > 0
      ? `${Math.round((1 - result.score.final / result.score.initial) * 100)}% lower score`
      : "score unchanged";

  const name = safeName(result.name);
  const banner = result.cancelled
    ? h("div", { class: "callout warn", role: "status" }, "Stopped early. This is the best sequence found so far.")
    : result.search.hit_time_limit
      ? h("div", { class: "callout warn", role: "status" }, "The time limit ended the search before it settled. Raise the limit or use a thorough search for more.")
      : null;
  const unresolved = result.valid ? null : h(
    "div", { class: "callout warn", role: "alert" },
    "Mandatory requirements remain unresolved. This candidate is not a valid design: " +
    result.unresolved_constraints.map((v) => `${v.criterion} (${v.count})`).join(", ") + ".",
  );

  const toolbar = h(
    "div",
    { class: "toolbar" },
    action("Copy DNA", (b) => void copyText(result.sequence).then((ok) => flash(b, ok ? "Copied" : "Copy failed")), true),
    action("Copy FASTA", (b) => void copyText(result.fasta).then((ok) => flash(b, ok ? "Copied" : "Copy failed"))),
    action("Download FASTA", () => download(`${name}.codonopt.fna`, result.fasta, "text/plain")),
    action("Download report", () =>
      download(`${name}.codonopt.json`, JSON.stringify(result, null, 2) + "\n", "application/json"),
    ),
  );

  replace(
    output,
    h(
      "div",
      { class: "result-head" },
      h("div", null, h("h2", { text: result.valid ? "Optimized sequence" : "Unresolved candidate" }), h("p", { class: "hint", text: `${result.name} · ${result.protein.length.toLocaleString()} residues · ${result.codon_usage_table} codon usage` })),
      toolbar,
    ),
    banner,
    unresolved,
    h(
      "div",
      { class: "stats" },
      stat("Issues found", int(issuesBefore), int(issuesAfter), issuesAfter === 0 ? "none left" : "still to resolve"),
      stat("Score", fmt(result.score.initial, 2), fmt(result.score.final, 2), reduction),
      stat("GC content", percent(result.initial.gc_fraction), percent(result.final.gc_fraction), `target ${percent(result.gc.target)}`),
      stat("Codon adaptation", fmt(result.initial.cai, 2), fmt(result.final.cai, 2), "1.0 is ideal"),
    ),
    card("DNA sequence", "The protein is identical to your input; only synonymous codons changed.", sequenceView(result)),
    card("Before and after", null, comparison(result)),
    card("GC content along the construct", `Each point is a ${result.gc.window} base window, including any fixed flanks.`, gcChart(result.gc, result.context_sequence.length)),
    card("Where the score comes from", "Lower is better. Each term is weighted as in your settings.", breakdown(result)),
    h(
      "p",
      { class: "run-info hint" },
      `${result.search.restarts_run} ${result.search.restarts_run === 1 ? "restart" : "restarts"} · ${int(result.search.greedy_sweeps)} sweeps · ${int(result.search.moves_evaluated)} codon swaps tried · ${fmt(result.search.seconds, 1)} s`,
    ),
  );
  output.focus({ preventScroll: true });
}
