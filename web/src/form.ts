import { $, h, replace } from "./dom";
import type { Kind, Meta } from "./types";

/** Display names for the score terms, in the order they are listed. */
const WEIGHT_LABELS: Record<string, string> = {
  gc_global: "Overall GC",
  gc_window: "Local GC",
  repeat: "Repeats",
  hairpin: "Hairpins",
  homopolymer: "Single-base runs",
  motif: "Forbidden sites",
  splice: "Splice sites",
  usage: "Codon usage",
};

const GROUP_LABELS: Record<string, string> = {
  golden_gate: "Golden Gate",
  common_cloning: "Common cloning sites",
  regulatory: "Poly(A) and TATA signals",
};

const SPLICE_LABELS: Record<string, string> = {
  auto: "Automatic",
  maxent: "MaxEntScan",
  motif: "Consensus screen",
  none: "Off",
};

const num = (id: string): HTMLInputElement => $<HTMLInputElement>(id);

/** Fill the dynamic parts of the form from what the server reports. */
export function populate(meta: Meta): void {
  $("version").textContent = `v${meta.version}`;

  const table = $<HTMLSelectElement>("table");
  replace(
    table,
    ...meta.tables.map((t) => h("option", { value: t.name, text: t.name })),
    h("option", { value: "custom", text: "Custom table…" }),
  );

  replace(
    $("avoid-groups"),
    ...Object.entries(meta.motif_groups)
      .filter(([, members]) => members.length)
      .map(([name, members]) =>
        h(
          "label",
          { class: "chip-check", title: members.join(", ") },
          h("input", { type: "checkbox", value: name, attrs: { name: "avoid" } }),
          h("span", { text: GROUP_LABELS[name] ?? name }),
        ),
      ),
  );

  replace(
    $("splice"),
    ...meta.splice_backends.map((name) => h("option", { value: name, text: SPLICE_LABELS[name] ?? name })),
  );

  replace(
    $("weights"),
    ...Object.keys(meta.weights).map((name) =>
      h(
        "label",
        { class: "field" },
        h("span", { text: WEIGHT_LABELS[name] ?? name }),
        h("input", {
          type: "number",
          id: `weight-${name}`,
          attrs: { min: "0", step: "any", inputmode: "decimal", "data-weight": name },
        }),
      ),
    ),
  );

  reset(meta);
}

/** Put every setting back to the server's defaults. The sequence is left alone. */
export function reset(meta: Meta): void {
  const d = meta.defaults;
  $<HTMLSelectElement>("table").value = meta.default_table;
  $<HTMLTextAreaElement>("table-text").value = "";
  syncTable();
  for (const box of document.querySelectorAll<HTMLInputElement>('#avoid-groups input[type="checkbox"]')) {
    box.checked = false;
  }
  $<HTMLInputElement>("avoid-extra").value = "";
  num("gc-target").value = String(d.gc_target);
  num("gc-low").value = String(d.gc_range[0]);
  num("gc-high").value = String(d.gc_range[1]);
  num("max-homopolymer").value = String(d.max_homopolymer);
  $<HTMLSelectElement>("splice").value = "auto";
  $<HTMLInputElement>("effort-standard").checked = true;
  num("seconds").value = String(d.seconds);
  num("seed").value = "";
  num("gc-window").value = String(d.gc_window);
  num("gc-step").value = String(d.gc_step);
  num("repeat-k").value = String(d.repeat_k);
  num("hairpin-stem").value = String(d.hairpin_stem);
  num("hairpin-min-loop").value = String(d.hairpin_min_loop);
  num("hairpin-max-loop").value = String(d.hairpin_max_loop);
  for (const [name, value] of Object.entries(meta.weights)) num(`weight-${name}`).value = String(value);
  describeEffort(meta);
  describeSplice(meta);
}

/** Show the custom-table box only when "Custom table…" is chosen. */
export function syncTable(): void {
  const custom = $<HTMLSelectElement>("table").value === "custom";
  $("custom-table").hidden = !custom;
}

/**
 * Why the form cannot be submitted as it stands, or null when it can.
 *
 * Only for what the server cannot see: an empty custom table reaches it as no
 * table at all, and it would quietly fall back to the default host.
 */
export function problem(): string | null {
  if ($<HTMLSelectElement>("table").value !== "custom") return null;
  if ($<HTMLTextAreaElement>("table-text").value.trim()) return null;
  return "Paste or open a codon usage table, or pick a built-in host instead of “Custom table…”.";
}

export function describeTable(meta: Meta): void {
  const name = $<HTMLSelectElement>("table").value;
  $("table-hint").textContent =
    name === "custom" ? "" : (meta.tables.find((t) => t.name === name)?.description ?? "");
}

export function describeEffort(meta: Meta): void {
  const level = selectedEffort();
  const e = meta.effort[level];
  if (!e) return;
  const parts = [`${e.max_sweeps} greedy passes`];
  if (e.restarts > 1) parts.push(`${e.restarts} independent restarts`);
  if (e.anneal_steps) parts.push(`${e.anneal_steps.toLocaleString()} annealing steps each`);
  $("effort-hint").textContent = `${parts.join(", ")}.`;
}

export function describeSplice(meta: Meta): void {
  const backend = $<HTMLSelectElement>("splice").value;
  const hint = $("splice-hint");
  if (backend === "maxent" && !meta.maxent_available) {
    hint.textContent = "MaxEntScan is not installed, so this will report an error. Pick Automatic or Consensus screen.";
  } else if (backend === "auto") {
    hint.textContent = meta.maxent_available
      ? "Uses MaxEntScan."
      : "MaxEntScan is not installed, so a simpler consensus screen is used.";
  } else if (backend === "motif") {
    hint.textContent = "A fast consensus screen. It flags likely sites rather than modelling splicing.";
  } else if (backend === "none") {
    hint.textContent = "Splice sites are ignored. Right for sequences that will not be spliced, and faster.";
  } else {
    hint.textContent = "";
  }
}

export function selectedEffort(): string {
  return document.querySelector<HTMLInputElement>('input[name="effort"]:checked')?.value ?? "standard";
}

export function selectedKind(): Kind {
  return (document.querySelector<HTMLInputElement>('input[name="kind"]:checked')?.value ?? "auto") as Kind;
}

/** A number box's value: undefined when empty (use the default), the raw text when not numeric. */
function read(id: string): number | string | undefined {
  const text = num(id).value.trim();
  if (text === "") return undefined;
  const value = Number(text);
  return Number.isFinite(value) ? value : text;
}

/** Collect the form into the request body the server expects. */
export function readPayload(): Record<string, unknown> {
  const avoid = [
    ...document.querySelectorAll<HTMLInputElement>('#avoid-groups input:checked'),
  ].map((box) => box.value);
  avoid.push(
    ...$<HTMLInputElement>("avoid-extra")
      .value.split(/[,\n]/)
      .map((item) => item.trim())
      .filter(Boolean),
  );

  const weights: Record<string, number | string> = {};
  for (const box of document.querySelectorAll<HTMLInputElement>("#weights input")) {
    const value = box.value.trim();
    if (value !== "") weights[box.dataset.weight ?? ""] = Number.isFinite(Number(value)) ? Number(value) : value;
  }

  const table = $<HTMLSelectElement>("table").value;
  const kind = selectedKind();
  return {
    sequence: $<HTMLTextAreaElement>("sequence").value,
    kind,
    add_stop: $<HTMLInputElement>("add-stop").checked,
    table: table === "custom" ? undefined : table,
    table_text: table === "custom" ? $<HTMLTextAreaElement>("table-text").value : undefined,
    avoid,
    gc_target: read("gc-target"),
    gc_range: [read("gc-low"), read("gc-high")],
    max_homopolymer: read("max-homopolymer"),
    splice_backend: $<HTMLSelectElement>("splice").value,
    effort: selectedEffort(),
    seconds: read("seconds"),
    seed: read("seed"),
    gc_window: read("gc-window"),
    gc_step: read("gc-step"),
    repeat_k: read("repeat-k"),
    hairpin_stem: read("hairpin-stem"),
    hairpin_min_loop: read("hairpin-min-loop"),
    hairpin_max_loop: read("hairpin-max-loop"),
    weights,
  };
}
