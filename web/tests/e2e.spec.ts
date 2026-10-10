import { expect, test, type Page } from "@playwright/test";

const SHOTS = process.env.CODONOPT_SHOTS ?? "";
const shot = async (page: Page, name: string): Promise<void> => {
  if (SHOTS) await page.screenshot({ path: `${SHOTS}/${name}.png`, fullPage: true });
};

/** Pick the quick setting so the suite stays fast. */
async function quick(page: Page): Promise<void> {
  await page.locator('label[for="effort-quick"]').click();
}

test.beforeEach(async ({ page }) => {
  await page.goto("/");
  await expect(page.locator("#run")).toBeEnabled({ timeout: 10_000 }).catch(() => undefined);
});

test("loads with no console errors and a disabled Optimise button", async ({ page }) => {
  const errors: string[] = [];
  page.on("console", (m) => m.type() === "error" && errors.push(m.text()));
  page.on("pageerror", (e) => errors.push(String(e)));
  await page.goto("/");
  await expect(page.locator("#empty")).toBeVisible();
  await expect(page.locator("#run")).toBeDisabled();
  await expect(page.locator("#table option").first()).toBeAttached();
  expect(errors).toEqual([]);
  await shot(page, "01-empty");
});

test("optimises the example and shows every result section", async ({ page }) => {
  await page.locator("#btn-example").click();
  await expect(page.locator("#run")).toBeEnabled();
  await expect(page.locator("#detected")).not.toBeEmpty();
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator(".progress")).toBeVisible();
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 45_000 });
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
});

test("keeps the protein unchanged (the sequence is a synonymous recoding)", async ({ page }) => {
  await page.locator("#sequence").fill("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV");
  await expect(page.locator("#detected")).toContainText(/protein/i);
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator(".seq .cell").first()).toBeVisible({ timeout: 45_000 });
  const shown = await page.locator(".line.aa").evaluateAll((rows) =>
    rows.map((r) => Array.from(r.querySelectorAll(".cell")).map((c) => c.textContent).join("")).join(""),
  );
  expect(shown.startsWith("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV")).toBe(true);
});

test("explains problems instead of failing silently", async ({ page }) => {
  await page.locator("#sequence").fill("MKT*AY1");
  await expect(page.locator("#detected")).toBeVisible();
  // Whatever the badge says, running must either be blocked or give a readable error.
  if (await page.locator("#run").isEnabled()) {
    await page.locator("#run").click();
    await expect(page.locator(".callout.error")).toBeVisible({ timeout: 15_000 });
  } else {
    await expect(page.locator("#detected")).toHaveText(/\S/);
  }
  await shot(page, "03-error");
});

test("can cancel a long run and keeps the best result so far", async ({ page }) => {
  await page.locator("#btn-example").click();
  await page.locator("#seconds").fill("60");
  await page.locator('label[for="effort-thorough"]').click();
  await page.locator("#run").click();
  await expect(page.locator("#cancel")).toBeVisible();
  await page.waitForTimeout(1500);
  await page.locator("#cancel").click();
  await expect(page.locator(".result-head")).toBeVisible({ timeout: 20_000 });
  await expect(page.locator(".callout.warn")).toContainText("Stopped early");
  await expect(page.locator("#run")).toBeEnabled();
  await expect(page.locator("#cancel")).toBeHidden();
});

test("offers several built-in codon-usage tables", async ({ page }) => {
  await page.locator("#btn-example").click();
  await page.locator("#table").selectOption({ index: 0 });
  const options = await page.locator("#table option").allTextContents();
  expect(options.length).toBeGreaterThan(1);
});

for (const [name, scheme] of [["light", "light"], ["dark", "dark"]] as const) {
  test(`renders in ${name} mode`, async ({ browser }) => {
    const context = await browser.newContext({ colorScheme: scheme, viewport: { width: 1280, height: 900 } });
    const page = await context.newPage();
    await page.goto("/");
    await page.locator("#btn-example").click();
    await quick(page);
    await page.locator("#run").click();
    await expect(page.locator("svg.chart")).toBeVisible({ timeout: 45_000 });
    await shot(page, `04-${name}`);
    await context.close();
  });
}

test("works at phone width without horizontal page scroll", async ({ browser }) => {
  const context = await browser.newContext({ viewport: { width: 390, height: 800 } });
  const page = await context.newPage();
  await page.goto("/");
  await page.locator("#btn-example").click();
  await quick(page);
  await page.locator("#run").click();
  await expect(page.locator("svg.chart")).toBeVisible({ timeout: 45_000 });
  const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  expect(overflow).toBeLessThanOrEqual(1);
  await shot(page, "05-phone");
  await context.close();
});
