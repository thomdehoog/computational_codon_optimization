// Build the interface: TypeScript and CSS in ./src become plain static files
// in ../src/codonopt/web, which are committed. Python users never need Node.
//
//   npm run build          write the files
//   npm run watch          rebuild on change
//   npm run check          fail if the committed files are out of date

import { build, context } from "esbuild";
import { watchFile } from "node:fs";
import { mkdir, readFile, readdir, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
const src = (name) => path.join(here, "src", name);
const outDir = path.resolve(here, "..", "src", "codonopt", "web");

const options = {
  bundle: true,
  minify: true,
  target: "es2020",
  legalComments: "none",
  logLevel: "warning",
};

const scripts = { ...options, entryPoints: [src("main.ts")], format: "iife", outfile: "app.js" };
const styles = { ...options, entryPoints: [src("styles.css")], outfile: "app.css" };

/** Files that ship as they are, rather than being compiled. */
const COPIED = ["index.html", "favicon.svg"];

/** Build everything in memory: a map of output file name to bytes. */
async function artifacts() {
  const files = new Map();
  for (const config of [scripts, styles]) {
    const result = await build({ ...config, write: false });
    if (result.warnings.length) {
      // The build is the gate CI relies on, so a warning must not slip past it.
      throw new Error(`esbuild reported ${result.warnings.length} warning(s); fix them or silence them explicitly.`);
    }
    for (const file of result.outputFiles) files.set(path.basename(file.path), Buffer.from(file.contents));
  }
  for (const name of COPIED) files.set(name, await readFile(src(name)));
  return files;
}

const mode = process.argv[2];

if (mode === "--watch") {
  await mkdir(outDir, { recursive: true });
  const watchers = await Promise.all([
    context({ ...scripts, outfile: path.join(outDir, "app.js") }),
    context({ ...styles, outfile: path.join(outDir, "app.css") }),
  ]);
  await Promise.all(watchers.map((w) => w.watch()));
  // index.html carries the token placeholder and the whole form skeleton, so
  // it is edited as often as the TypeScript and has to be watched too.
  const copy = async () => {
    for (const name of COPIED) await writeFile(path.join(outDir, name), await readFile(src(name)));
  };
  await copy();
  for (const name of COPIED) {
    watchFile(src(name), { interval: 200 }, () => {
      void copy().then(() => console.log(`copied ${name}`));
    });
  }
  console.log(`watching src/ and writing to ${outDir}`);
} else if (mode === "--check") {
  const stale = [];
  const built = await artifacts();
  for (const [name, bytes] of built) {
    const current = await readFile(path.join(outDir, name)).catch(() => null);
    if (!current || !current.equals(bytes)) stale.push(name);
  }
  // A committed file the build no longer produces would otherwise ship forever.
  const present = await readdir(outDir).catch(() => []);
  for (const name of present) {
    if (!built.has(name)) stale.push(`${name} (no longer built)`);
  }
  if (stale.length) {
    console.error(`Out of date: ${stale.join(", ")}\nRun "npm run build" in web/ and commit the result.`);
    process.exit(1);
  }
  console.log("The committed interface matches its sources.");
} else {
  await mkdir(outDir, { recursive: true });
  let total = 0;
  for (const [name, bytes] of await artifacts()) {
    await writeFile(path.join(outDir, name), bytes);
    total += bytes.length;
    console.log(`  ${name.padEnd(14)} ${(bytes.length / 1024).toFixed(1)} kB`);
  }
  console.log(`built ${(total / 1024).toFixed(1)} kB into ${path.relative(process.cwd(), outDir) || "."}`);
}
