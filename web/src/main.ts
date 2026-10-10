import { Api, ApiError } from "./api";
import { describe } from "./detect";
import { $ } from "./dom";
import {
  describeEffort,
  describeSplice,
  describeTable,
  populate,
  readPayload,
  reset,
  selectedKind,
  syncTable,
} from "./form";
import { Progress, showError, showResult } from "./results";
import type { Job, Meta } from "./types";

const POLL_MS = 350;

async function main(): Promise<void> {
  const token = document.querySelector<HTMLMetaElement>('meta[name="codonopt-token"]')?.content ?? "";
  const api = new Api(token);
  const output = $("output");

  let meta: Meta;
  try {
    meta = await api.meta();
  } catch (error) {
    showError(output, error instanceof Error ? error.message : String(error));
    return;
  }
  populate(meta);
  describeTable(meta);

  const sequence = $<HTMLTextAreaElement>("sequence");
  const detected = $("detected");
  const run = $<HTMLButtonElement>("run");
  const cancel = $<HTMLButtonElement>("cancel");
  const stopRow = $("add-stop-wrap");
  let running: { id: string; timer: number } | null = null;

  /* ----- sequence box ----- */

  const refresh = (): void => {
    const info = describe(sequence.value, selectedKind(), meta.limits.max_residues);
    detected.textContent = info.message;
    detected.className = `detected ${info.state}`;
    stopRow.hidden = info.kind === "dna";
    run.disabled = running !== null || info.state === "empty";
  };
  sequence.addEventListener("input", refresh);
  for (const radio of document.querySelectorAll('input[name="kind"]')) radio.addEventListener("change", refresh);

  const load = async (file: File | undefined, into: HTMLTextAreaElement): Promise<void> => {
    if (!file) return;
    into.value = await file.text();
    into.dispatchEvent(new Event("input", { bubbles: true }));
  };
  const fileInput = $<HTMLInputElement>("file-input");
  $("btn-open").addEventListener("click", () => fileInput.click());
  fileInput.addEventListener("change", () => {
    void load(fileInput.files?.[0], sequence);
    fileInput.value = "";
  });
  sequence.addEventListener("dragover", (event) => event.preventDefault());
  sequence.addEventListener("drop", (event) => {
    event.preventDefault();
    void load(event.dataTransfer?.files[0], sequence);
  });

  $("btn-example").addEventListener("click", () => {
    sequence.value = `>${meta.example.name}\n${meta.example.protein}\n`;
    $<HTMLInputElement>("kind-auto").checked = true;
    refresh();
    sequence.focus();
  });
  $("btn-clear").addEventListener("click", () => {
    sequence.value = "";
    refresh();
    sequence.focus();
  });

  /* ----- settings ----- */

  const tableSelect = $<HTMLSelectElement>("table");
  tableSelect.addEventListener("change", () => {
    syncTable();
    describeTable(meta);
  });
  const tableText = $<HTMLTextAreaElement>("table-text");
  const tableFile = $<HTMLInputElement>("table-file");
  $("btn-open-table").addEventListener("click", () => tableFile.click());
  tableFile.addEventListener("change", () => {
    void load(tableFile.files?.[0], tableText);
    tableFile.value = "";
  });
  $("splice").addEventListener("change", () => describeSplice(meta));
  for (const radio of document.querySelectorAll('input[name="effort"]')) {
    radio.addEventListener("change", () => describeEffort(meta));
  }
  $("btn-reset").addEventListener("click", () => {
    reset(meta);
    describeTable(meta);
  });

  /* ----- running ----- */

  const setRunning = (active: boolean): void => {
    cancel.hidden = !active;
    cancel.disabled = false;
    run.hidden = active;
    refresh();
  };

  const finish = (job: Job): void => {
    if (running) window.clearInterval(running.timer);
    running = null;
    setRunning(false);
    if (job.status === "error") showError(output, job.error ?? "The optimization failed.");
    else if (job.result) showResult(output, job.result, meta);
  };

  $("form").addEventListener("submit", async (event) => {
    event.preventDefault();
    if (running || run.disabled) return;
    let job: Job;
    try {
      job = await api.start(readPayload());
    } catch (error) {
      showError(output, error instanceof ApiError ? error.message : String(error));
      return;
    }
    const progress = new Progress(output);
    progress.update(job);
    running = { id: job.id, timer: 0 };
    setRunning(true);
    running.timer = window.setInterval(async () => {
      if (!running) return;
      try {
        const latest = await api.poll(running.id);
        if (latest.status === "running") progress.update(latest);
        else finish(latest);
      } catch (error) {
        if (running) window.clearInterval(running.timer);
        running = null;
        setRunning(false);
        showError(output, error instanceof Error ? error.message : String(error));
      }
    }, POLL_MS);
  });

  cancel.addEventListener("click", async () => {
    if (!running) return;
    cancel.disabled = true;
    cancel.textContent = "Stopping…";
    try {
      await api.cancel(running.id);
    } catch {
      /* the next poll reports whatever state the job is in */
    }
  });
  new MutationObserver(() => {
    if (cancel.hidden) cancel.textContent = "Stop";
  }).observe(cancel, { attributes: true, attributeFilter: ["hidden"] });

  refresh();
}

void main();
