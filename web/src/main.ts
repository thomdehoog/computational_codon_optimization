import { Api, ApiError } from "./api";
import { describe } from "./detect";
import { $ } from "./dom";
import {
  describeEffort,
  describeSplice,
  describeTable,
  populate,
  problem,
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

  // The form starts inert: `populate` fills it from the server's defaults, and
  // anything typed before that lands would be overwritten without warning.
  const form = $<HTMLFormElement>("form");
  let meta: Meta;
  try {
    meta = await api.meta();
  } catch (error) {
    showError(output, error instanceof Error ? error.message : String(error));
    return;
  }
  populate(meta);
  describeTable(meta);
  form.removeAttribute("inert");
  form.removeAttribute("aria-busy");

  const sequence = $<HTMLTextAreaElement>("sequence");
  const detected = $("detected");
  const run = $<HTMLButtonElement>("run");
  const cancel = $<HTMLButtonElement>("cancel");
  const stopRow = $("add-stop-wrap");
  /** The one run that owns the screen. Anything else that arrives is stale. */
  interface Run {
    id: string;
    timer: number;
    progress: Progress;
  }
  let running: Run | null = null;
  let starting = false;

  /* ----- sequence box ----- */

  const refresh = (): void => {
    const info = describe(
      sequence.value,
      selectedKind(),
      meta.limits.max_residues,
      $<HTMLInputElement>("add-stop").checked,
    );
    detected.textContent = info.message;
    detected.className = `detected ${info.state}`;
    stopRow.hidden = info.kind === "dna";
    run.disabled = running !== null || starting || info.state === "empty";
  };
  sequence.addEventListener("input", refresh);
  $("add-stop").addEventListener("change", refresh);
  for (const radio of document.querySelectorAll('input[name="kind"]')) radio.addEventListener("change", refresh);

  // Generous next to the 15,000-residue limit, small enough that a misdropped
  // genome cannot hang the tab.
  const MAX_FILE_BYTES = 8_000_000;

  const load = async (file: File | undefined, into: HTMLTextAreaElement): Promise<void> => {
    if (!file) return;
    if (file.size > MAX_FILE_BYTES) {
      showError(output, `${file.name} is ${Math.round(file.size / 1e6)} MB; open a file under 8 MB.`);
      return;
    }
    let text: string;
    try {
      text = await file.text();
    } catch {
      // A folder dropped on the box, or a file that moved since it was picked.
      showError(output, `Could not read ${file.name || "that file"}. Open a text file, not a folder.`);
      return;
    }
    into.value = text;
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
    describeSplice(meta);
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
    cancel.textContent = "Stop";
    run.hidden = active;
    refresh();
  };

  /** Announce a change a sighted user sees in the results region. */
  const say = (message: string): void => {
    $("status").textContent = message;
  };

  /** Let go of a run: stop polling it, and clear the state if it is current. */
  const release = (job: Run): void => {
    window.clearInterval(job.timer);
    if (running !== job) return; // a later run owns the screen now
    running = null;
    setRunning(false);
  };

  const finish = (job: Run, latest: Job): void => {
    const current = running === job;
    release(job);
    if (!current) return; // a stale poll must not paint over a newer run
    if (latest.status === "error") {
      say("The optimization failed.");
      showError(output, latest.error ?? "The optimization failed.");
    } else if (latest.result) {
      say(!latest.result.valid ? "Finished with unresolved mandatory requirements." :
        latest.status === "cancelled" ? "Stopped early; showing the best sequence found." : "Optimization finished.");
      showResult(output, latest.result, meta);
    }
  };

  const poll = async (job: Run): Promise<void> => {
    if (running !== job) return;
    try {
      const latest = await api.poll(job.id);
      if (running !== job) return; // the answer arrived after we moved on
      if (latest.status === "running") job.progress.update(latest);
      else finish(job, latest);
    } catch (error) {
      if (running !== job) return;
      release(job);
      showError(output, error instanceof Error ? error.message : String(error));
    }
  };

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    // `starting` closes the window between this click and the server's reply,
    // where a second click would otherwise start a second job.
    if (running || starting || run.disabled) return;
    const wrong = problem();
    if (wrong) {
      showError(output, wrong);
      return;
    }
    starting = true;
    run.disabled = true;
    let started: Job;
    try {
      started = await api.start(readPayload());
    } catch (error) {
      showError(output, error instanceof ApiError ? error.message : String(error));
      return;
    } finally {
      starting = false;
      refresh();
    }
    say("Optimizing.");
    const progress = new Progress(output);
    progress.update(started);
    const job: Run = { id: started.id, timer: 0, progress };
    running = job;
    setRunning(true);
    job.timer = window.setInterval(() => void poll(job), POLL_MS);
  });

  cancel.addEventListener("click", async () => {
    const job = running;
    if (!job) return;
    cancel.disabled = true;
    cancel.textContent = "Stopping…";
    try {
      await api.cancel(job.id);
    } catch (error) {
      if (running !== job) return;
      // Stop is the only way out of a long run, so never leave it dead.
      cancel.disabled = false;
      cancel.textContent = "Stop";
      job.progress.note(
        `${error instanceof Error ? error.message : String(error)} The run is still going; try Stop again.`,
      );
    }
  });

  refresh();
}

// Without this, any failure while wiring the page would leave it silent and
// dead, with no hint of what went wrong.
void main().catch((error: unknown) => {
  const output = document.getElementById("output");
  const message = error instanceof Error ? error.message : String(error);
  if (output) showError(output, `The interface failed to start: ${message}`);
  else document.title = "codonopt - failed to start";
  console.error(error);
});
