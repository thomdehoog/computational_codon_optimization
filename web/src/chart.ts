import { fmt, h, int, svg } from "./dom";
import type { GcProfile } from "./types";

const W = 720;
const H = 236;
const M = { left: 46, right: 14, top: 12, bottom: 34 };
const PLOT_W = W - M.left - M.right;
const PLOT_H = H - M.top - M.bottom;

const NICE_STEPS = [50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000, 50000];

/** Pick a tick spacing that gives roughly five to eight ticks. */
function tickStep(length: number): number {
  return NICE_STEPS.find((step) => length / step <= 8) ?? NICE_STEPS[NICE_STEPS.length - 1];
}

/**
 * GC content along the sequence, before and after, against the allowed range.
 *
 * Built as an SVG from attributes only (the page's content security policy
 * forbids inline styles); colours and dashes come from the stylesheet.
 */
export function gcChart(gc: GcProfile, length: number): HTMLElement {
  const total = Math.max(length, 1);
  const x = (position: number): number => M.left + (position / total) * PLOT_W;
  const y = (fraction: number): number => M.top + (1 - fraction) * PLOT_H;
  const centre = (start: number): number => Math.min(start + gc.window / 2, total);

  const line = (values: number[]): string =>
    values
      .map((value, i) => `${i ? "L" : "M"}${x(centre(gc.starts[i])).toFixed(1)} ${y(value).toFixed(1)}`)
      .join("");

  const chart = svg("svg", {
    viewBox: `0 0 ${W} ${H}`,
    role: "img",
    class: "chart",
    "aria-label":
      `GC content in ${gc.window} base windows along the sequence, before and after optimisation, ` +
      `against the allowed range of ${fmt(gc.range[0] * 100, 0)} to ${fmt(gc.range[1] * 100, 0)} percent.`,
  });

  // Allowed range, then grid and labels, so the lines sit on top.
  chart.append(
    svg("rect", {
      class: "band",
      x: M.left,
      y: y(gc.range[1]),
      width: PLOT_W,
      height: Math.max(0, y(gc.range[0]) - y(gc.range[1])),
    }),
  );
  for (const tick of [0, 0.25, 0.5, 0.75, 1]) {
    chart.append(
      svg("line", { class: "grid", x1: M.left, x2: W - M.right, y1: y(tick), y2: y(tick) }),
      svg("text", { class: "tick", x: M.left - 8, y: y(tick) + 4, "text-anchor": "end" }, `${tick * 100}%`),
    );
  }
  const step = tickStep(total);
  for (let position = 0; position <= total; position += step) {
    chart.append(
      svg("line", { class: "grid faint", x1: x(position), x2: x(position), y1: M.top, y2: M.top + PLOT_H }),
      svg("text", { class: "tick", x: x(position), y: H - 12, "text-anchor": "middle" }, int(position)),
    );
  }
  chart.append(
    svg("line", { class: "target", x1: M.left, x2: W - M.right, y1: y(gc.target), y2: y(gc.target) }),
  );

  if (gc.starts.length > 1) {
    chart.append(
      svg("path", { class: "line before", d: line(gc.initial) }),
      svg("path", { class: "line after", d: line(gc.final) }),
    );
  } else if (gc.starts.length === 1) {
    chart.append(
      svg("circle", { class: "dot before", cx: x(centre(gc.starts[0])), cy: y(gc.initial[0]), r: 4 }),
      svg("circle", { class: "dot after", cx: x(centre(gc.starts[0])), cy: y(gc.final[0]), r: 4 }),
    );
  }

  const cursor = svg("line", { class: "cursor", x1: 0, x2: 0, y1: M.top, y2: M.top + PLOT_H });
  cursor.setAttribute("visibility", "hidden");
  const overlay = svg("rect", {
    class: "overlay",
    x: M.left,
    y: M.top,
    width: PLOT_W,
    height: PLOT_H,
  });
  chart.append(cursor, overlay);

  const idle = "Hover the chart to read a window.";
  const readout = h("p", { class: "readout", text: idle });

  overlay.addEventListener("pointermove", (event: PointerEvent) => {
    const box = chart.getBoundingClientRect();
    const position = ((event.clientX - box.left) / box.width) * W;
    const bp = ((position - M.left) / PLOT_W) * total;
    let best = 0;
    for (let i = 1; i < gc.starts.length; i++) {
      if (Math.abs(centre(gc.starts[i]) - bp) < Math.abs(centre(gc.starts[best]) - bp)) best = i;
    }
    if (!gc.starts.length) return;
    const at = x(centre(gc.starts[best]));
    cursor.setAttribute("x1", String(at));
    cursor.setAttribute("x2", String(at));
    cursor.setAttribute("visibility", "visible");
    const start = gc.starts[best];
    readout.textContent =
      `Bases ${int(start + 1)}–${int(Math.min(start + gc.window, total))}: ` +
      `before ${fmt(gc.initial[best] * 100, 0)}%, after ${fmt(gc.final[best] * 100, 0)}%`;
  });
  overlay.addEventListener("pointerleave", () => {
    cursor.setAttribute("visibility", "hidden");
    readout.textContent = idle;
  });

  const legend = h(
    "ul",
    { class: "legend" },
    h("li", null, h("span", { class: "swatch before" }), "Before"),
    h("li", null, h("span", { class: "swatch after" }), "After"),
    h("li", null, h("span", { class: "swatch band" }), "Allowed range"),
    h("li", null, h("span", { class: "swatch target" }), `Target ${fmt(gc.target * 100, 0)}%`),
  );

  return h("figure", { class: "figure" }, chart, legend, readout);
}
