import { defineConfig } from "@playwright/test";

// A port of this run's own, so a stale `codonopt gui` left running elsewhere
// can never be adopted and tested in place of the current build.
const PORT = Number(process.env.CODONOPT_PORT ?? 8700 + (process.pid % 900));

export default defineConfig({
  testDir: "tests",
  timeout: 90_000,
  fullyParallel: false,
  workers: 1,
  retries: process.env.CI ? 2 : 0,
  // The "github" reporter turns failures into annotations on the run itself,
  // which can be read without digging through the log.
  reporter: process.env.CI ? [["github"], ["list"], ["html", { open: "never" }]] : "list",
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
    launchOptions: process.env.CODONOPT_CHROMIUM ? { executablePath: process.env.CODONOPT_CHROMIUM } : {},
    // So a failure in CI can be read afterwards instead of guessed at.
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
  },
  // The real Python server, serving the committed build: what users get.
  webServer: {
    command: `python -m codonopt gui --no-window --no-browser --port ${PORT}`,
    url: `http://127.0.0.1:${PORT}/`,
    reuseExistingServer: false,
    stdout: "pipe",
    stderr: "pipe",
    timeout: 60_000,
  },
});
