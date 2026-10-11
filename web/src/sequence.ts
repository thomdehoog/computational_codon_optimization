import { h, int } from "./dom";
import type { Result } from "./types";

const CODONS_PER_ROW = 20;

/**
 * The optimized sequence, codon by codon, with the protein above it.
 *
 * Codons that differ from the starting sequence can be highlighted, and any
 * codon that still overlaps a forbidden motif is marked so it cannot be missed.
 */
export function sequenceView(result: Result): HTMLElement {
  const final = result.sequence;
  const initial = result.initial_sequence;
  const codons = final.length / 3;

  const flagged = new Set<number>();
  for (const hit of result.motif_hits) {
    for (const position of hit.positions) {
      const start = position - result.final.cds_start;
      for (let i = Math.max(0, Math.floor(start / 3)); i <= Math.min(codons - 1, Math.floor((start + hit.length - 1) / 3)); i++) {
        flagged.add(i);
      }
    }
  }

  let changed = 0;
  const rows: HTMLElement[] = [];
  const labelWidth = String(final.length).length;
  for (let first = 0; first < codons; first += CODONS_PER_ROW) {
    const last = Math.min(first + CODONS_PER_ROW, codons);
    const aminoAcids = h("div", { class: "line aa", attrs: { "aria-hidden": "true" } });
    const bases = h("div", { class: "line dna" });
    for (let i = first; i < last; i++) {
      const codon = final.slice(i * 3, i * 3 + 3);
      const before = initial.slice(i * 3, i * 3 + 3);
      const differs = codon !== before;
      if (differs) changed++;
      const classes = ["cell", differs ? "chg" : "", flagged.has(i) ? "hit" : ""].filter(Boolean).join(" ");
      aminoAcids.append(h("span", { class: "cell" }, result.protein[i] ?? ""));
      bases.append(
        h("span", differs ? { class: classes, title: `${before} → ${codon}` } : { class: classes }, codon),
      );
    }
    rows.push(
      h(
        "div",
        { class: "block" },
        h("span", { class: "pos", text: String(first * 3 + 1).padStart(labelWidth, " ") }),
        h("div", { class: "lines" }, aminoAcids, bases),
      ),
    );
  }

  const view = h("div", { class: "seq show-changes show-aa", attrs: { tabindex: "0", "aria-label": "Optimized DNA sequence" } }, ...rows);

  const toggle = (label: string, className: string, checked: boolean): HTMLElement =>
    h(
      "label",
      { class: "check" },
      h("input", {
        type: "checkbox",
        checked,
        on: { change: (event) => view.classList.toggle(className, (event.target as HTMLInputElement).checked) },
      }),
      ` ${label}`,
    );
  if (!result.motif_hits.length && !changed) view.classList.remove("show-changes");

  const summary =
    changed === 0
      ? "No codon needed to change."
      : `${int(changed)} of ${int(codons)} codons changed (${Math.round((changed / codons) * 100)}%).`;
  const warn = result.motif_hits.length
    ? h(
        "p",
        { class: "note warn", role: "status" },
        `Still contains ${result.motif_hits.map((hit) => `${hit.name} (${hit.positions.length})`).join(", ")}; ` +
          "matches are in the scored construct; overlapping CDS codons are marked.",
      )
    : null;

  return h(
    "div",
    { class: "seq-wrap" },
    h(
      "div",
      { class: "seq-tools" },
      h("p", { class: "hint", text: summary }),
      h(
        "div",
        { class: "toggles" },
        toggle("Highlight changes", "show-changes", view.classList.contains("show-changes")),
        toggle("Show protein", "show-aa", true),
      ),
    ),
    warn,
    view,
  );
}
