import { expect, test, type Page } from "@playwright/test";

const SHOTS = process.env.CODONOPT_SHOTS ?? "";
const shot = async (page: Page, name: string): Promise<void> => {
  if (SHOTS) await page.screenshot({ path: `${SHOTS}/${name}.png`, fullPage: true });
};

/**
 * Open the page and wait until it is actually usable.
 *
 * The form carries `inert` until the defaults have arrived, which is the only
 * honest readiness signal: until then the buttons exist but do nothing.
 * Returns the list of browser errors seen, which every test asserts is empty.
 */
async function open(page: Page): Promise<string[]> {
  const errors: string[] = [];
  page.on("console", (message) => {
    // A deliberate 400 is logged by the browser as a console error; what we are
    // watching for is thrown exceptions and content-policy refusals.
    if (message.type() === "error" && !/Failed to load resource/.test(message.text())) {
      errors.push(message.text());
    }
  });
  page.on("pageerror", (error) => errors.push(String(error)));
  await page.goto("/");
  await page.waitForFunction(() => !document.getElementById("form")?.hasAttribute("inert"));
  return errors;
}

/** Pick the quick setting so the suite stays fast. */
async function quick(page: Page): Promise<void> {
  await page.locator('label[for="effort-quick"]').click();
}

test("loads with no console errors and a disabled Optimize button", async ({ page }) => {
  const errors = await open(page);
  await expect(page.locator("#empty")).toBeVisible();
  await expect(page.locator("#run")).toBeDisabled();
  await expect(page.locator("#table option")).not.toHaveCount(0);
  await shot(page, "01-empty");
  expect(errors).toEqual([]);
});

test("optimises the example and shows every result section", async ({ page }) => {
  const errors = await open(page);
  await page.locator("#btn-example").click();
  await expect(page.locator("#run")).toBeEnabled();
  await expect(page.locator("#detected")).not.toBeEmpty();
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator(".progress")).toBeVisible();
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 60_000 });
  await expect(page.locator(".stats .stat").first()).toBeVisible();
  await expect(page.locator("table.compare")).toBeVisible();
  await expect(page.locator("svg.chart")).toBeVisible();
  await expect(page.locator(".seq .cell").first()).toBeVisible();
  await expect(page.locator(".callout.error")).toHaveCount(0);

  // The chart reads out a window on hover.
  await page.locator("svg.chart .overlay").hover({ position: { x: 200, y: 40 } });
  await expect(page.locator(".readout")).toContainText("Bases");

  // Download works (a blob, so it has to survive the content security policy).
  const [download] = await Promise.all([
    page.waitForEvent("download"),
    page.getByRole("button", { name: /download fasta/i }).click(),
  ]);
  expect(download.suggestedFilename()).toMatch(/\.(fna|fa|fasta)$/);
  await shot(page, "02-result");
  expect(errors).toEqual([]);
});

test("keeps the protein unchanged (the sequence is a synonymous recoding)", async ({ page }) => {
  const protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV";
  const errors = await open(page);
  await page.locator("#sequence").fill(protein);
  await expect(page.locator("#detected")).toContainText(/protein/i);
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator(".seq .cell").first()).toBeVisible({ timeout: 60_000 });
  const shown = await page.locator(".line.aa").evaluateAll((rows) =>
    rows.map((row) => Array.from(row.querySelectorAll(".cell")).map((cell) => cell.textContent).join("")).join(""),
  );
  // Compared as text, so a failure prints the protein that was actually drawn.
  expect(shown).toBe(`${protein}*`);
  expect(errors).toEqual([]);
});

test("labels an unavoidable motif as an unresolved candidate", async ({ page }) => {
  const errors = await open(page);
  await expect(page.locator("#splice")).toHaveValue("none");
  await page.locator("#sequence").fill("MW");
  await page.locator("#avoid-extra").fill("ATG");
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator(".result-head h2")).toHaveText("Unresolved candidate", { timeout: 60_000 });
  await expect(page.locator(".callout.warn[role=alert]")).toContainText("not a valid design");
  await expect(page.locator("#status")).toContainText("unresolved mandatory requirements");
  await expect(page.locator(".seq .line.dna")).toContainText("ATG");
  expect(errors).toEqual([]);
});

test("explains problems instead of failing silently", async ({ page }) => {
  const errors = await open(page);
  await page.locator("#sequence").fill("MKT*AY1");
  await expect(page.locator("#detected")).toBeVisible();
  // Whatever the badge says, running must either be blocked or give a readable error.
  if (await page.locator("#run").isEnabled()) {
    await page.locator("#run").click();
    await expect(page.locator(".callout.error")).toBeVisible({ timeout: 20_000 });
  } else {
    await expect(page.locator("#detected")).toHaveText(/\S/);
  }
  await shot(page, "03-error");
  expect(errors).toEqual([]);
});

test("can cancel a long run and hands the controls back", async ({ page }) => {
  const errors = await open(page);
  await page.locator("#btn-example").click();
  const example = await page.locator("#sequence").inputValue();
  await page.locator("#sequence").fill(example.repeat(20));
  await page.locator("#seconds").fill("60");
  await page.locator('label[for="effort-thorough"]').click();
  await page.locator("#run").click();
  await expect(page.locator("#cancel")).toBeVisible();
  await page.locator("#cancel").click();
  // What matters is that Stop ends the run and restores the form, whatever the
  // search had managed by then.
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 30_000 });
  await expect(page.locator("#run")).toBeEnabled();
  await expect(page.locator("#cancel")).toBeHidden();
  await expect(page.locator("#run")).toHaveText(/optimize/i);
  expect(errors).toEqual([]);
});

test("two submits in the same moment start only one run", async ({ page }) => {
  const errors = await open(page);
  // Counted at the wire, because the page only ever shows one card either way:
  // the damage of a second job is invisible until its result overwrites yours.
  let started = 0;
  page.on("request", (request) => {
    if (request.method() === "POST" && new URL(request.url()).pathname === "/api/jobs") started++;
  });
  await page.locator("#btn-example").click();
  await quick(page);
  // A double click, or a held Enter key, reaches the handler again before the
  // server has answered the first time.
  await page.evaluate(() => {
    const form = document.getElementById("form") as HTMLFormElement;
    form.requestSubmit();
    form.requestSubmit();
    form.requestSubmit();
  });
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 60_000 });
  expect(started).toBe(1);
  await expect(page.locator(".result-head")).toHaveCount(1);
  await expect(page.locator("#run")).toBeEnabled();
  await expect(page.locator("#cancel")).toBeHidden();
  // A second run must still work, and replace the first result.
  await page.locator("#run").click();
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 60_000 });
  await expect(page.locator(".result-head")).toHaveCount(1);
  expect(started).toBe(2);
  expect(errors).toEqual([]);
});

test("refuses an empty custom codon table instead of using another host", async ({ page }) => {
  const errors = await open(page);
  await page.locator("#btn-example").click();
  await page.locator("#table").selectOption("custom");
  await expect(page.locator("#custom-table")).toBeVisible();
  await page.locator("#run").click();
  await expect(page.locator(".callout.error")).toContainText(/custom table/i);
  await expect(page.locator(".progress")).toHaveCount(0);
  expect(errors).toEqual([]);
});

test("offers every built-in codon usage table, plus a custom one", async ({ page }) => {
  const errors = await open(page);
  const names = await page.locator("#table option").evaluateAll((options) =>
    options.map((option) => (option as HTMLOptionElement).value),
  );
  const expected = await page.evaluate(async () => {
    const token = document.querySelector<HTMLMetaElement>('meta[name="codonopt-token"]')?.content ?? "";
    const meta = await fetch("/api/meta", { headers: { "X-Codonopt-Token": token } }).then((r) => r.json());
    return meta.tables.map((table: { name: string }) => table.name);
  });
  expect(expected.length).toBeGreaterThan(0);
  expect(names).toEqual([...expected, "custom"]);
  expect(errors).toEqual([]);
});

test("warns before sending a sequence that is too long", async ({ page }) => {
  const errors = await open(page);
  const limit = await page.evaluate(async () => {
    const token = document.querySelector<HTMLMetaElement>('meta[name="codonopt-token"]')?.content ?? "";
    const meta = await fetch("/api/meta", { headers: { "X-Codonopt-Token": token } }).then((r) => r.json());
    return meta.limits.max_residues as number;
  });
  // DNA used to slip past the check entirely and show a green badge. The badge
  // only warns: the server stays the authority on what it will accept.
  await page.locator("#sequence").fill("ATG".repeat(limit + 10));
  await expect(page.locator("#detected")).toHaveClass(/error/);
  await expect(page.locator("#detected")).toContainText(/codons/);
  await expect(page.locator("#detected")).toContainText(String(limit.toLocaleString()));

  // And a protein one residue under the limit is fine until "Add stop codon"
  // pushes it over, which the badge used to miss.
  await page.locator("#sequence").fill("M".repeat(limit));
  await page.locator("#add-stop").uncheck();
  await expect(page.locator("#detected")).toHaveClass(/ok/);
  await page.locator("#add-stop").check();
  await expect(page.locator("#detected")).toHaveClass(/error/);
  expect(errors).toEqual([]);
});

for (const scheme of ["light", "dark"] as const) {
  test(`renders in ${scheme} mode`, async ({ browser }) => {
    const context = await browser.newContext({ colorScheme: scheme, viewport: { width: 1280, height: 900 } });
    const page = await context.newPage();
    const errors = await open(page);
    await page.locator("#btn-example").click();
    await quick(page);
    await page.locator("#run").click();
    await expect(page.locator("svg.chart")).toBeVisible({ timeout: 60_000 });
    await shot(page, `04-${scheme}`);
    expect(errors).toEqual([]);
    await context.close();
  });
}

test("works at phone width without horizontal page scroll", async ({ browser }) => {
  const context = await browser.newContext({ viewport: { width: 390, height: 800 } });
  const page = await context.newPage();
  const errors = await open(page);
  await page.locator("#btn-example").click();
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator("svg.chart")).toBeVisible({ timeout: 60_000 });
  const overflow = await page.evaluate(
    () => document.documentElement.scrollWidth - document.documentElement.clientWidth,
  );
  expect(overflow).toBeLessThanOrEqual(1);
  await shot(page, "05-phone");
  expect(errors).toEqual([]);
  await context.close();
});
