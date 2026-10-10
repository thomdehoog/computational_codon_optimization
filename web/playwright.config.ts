import { defineConfig } from "@playwright/test";

const PORT = 8765;

export default defineConfig({
  testDir: "tests",
  timeout: 60_000,
  fullyParallel: false,
  workers: 1,
  reporter: "list",
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
    launchOptions: process.env.CODONOPT_CHROMIUM ? { executablePath: process.env.CODONOPT_CHROMIUM } : {},
  },
  // The real Python server, serving the committed build: what users get.
  webServer: {
    command: `python -m codonopt gui --no-window --no-browser --port ${PORT}`,
    url: `http://127.0.0.1:${PORT}/`,
    reuseExistingServer: !process.env.CI,
    timeout: 30_000,
  },
});
