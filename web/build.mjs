// Build the interface: TypeScript and CSS in ./src become plain static files
// in ../src/codonopt/web, which are committed. Python users never need Node.
//
//   npm run build          write the files
//   npm run watch          rebuild on change
//   npm run check          fail if the committed files are out of date

import { build, context } from "esbuild";
import { mkdir, readFile, writeFile } from "node:fs/promises";
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

/** Build everything in memory: a map of output file name to bytes. */
async function artifacts() {
  const files = new Map();
  for (const config of [scripts, styles]) {
    const result = await build({ ...config, write: false });
    for (const file of result.outputFiles) files.set(path.basename(file.path), Buffer.from(file.contents));
  }
  for (const name of ["index.html", "favicon.svg"]) files.set(name, await readFile(src(name)));
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
  for (const name of ["index.html", "favicon.svg"]) await writeFile(path.join(outDir, name), await readFile(src(name)));
  console.log(`watching; writing to ${outDir} (index.html and favicon.svg were copied once)`);
} else if (mode === "--check") {
  const stale = [];
  for (const [name, bytes] of await artifacts()) {
    const current = await readFile(path.join(outDir, name)).catch(() => null);
    if (!current || !current.equals(bytes)) stale.push(name);
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
