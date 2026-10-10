/** A very small DOM toolkit: enough to build the results without a framework. */

export type Child = Node | string | number | null | undefined | false;

export interface Props {
  class?: string;
  text?: string;
  title?: string;
  id?: string;
  type?: string;
  role?: string;
  hidden?: boolean;
  disabled?: boolean;
  value?: string;
  checked?: boolean;
  attrs?: Record<string, string>;
  on?: Partial<{ [K in keyof HTMLElementEventMap]: (event: HTMLElementEventMap[K]) => void }>;
}

function append(parent: Node, children: Child[]): void {
  for (const child of children) {
    if (child === null || child === undefined || child === false) continue;
    parent.appendChild(typeof child === "object" ? child : document.createTextNode(String(child)));
  }
}

/** Create an HTML element: `h("p", { class: "hint" }, "text", child)`. */
export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  props?: Props | null,
  ...children: Child[]
): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  if (props) {
    const { attrs, on, text, ...rest } = props;
    for (const [key, value] of Object.entries(rest)) {
      if (value === undefined || value === false) continue;
      if (key === "class") el.className = String(value);
      else if (typeof value === "boolean") el.setAttribute(key, "");
      else el.setAttribute(key, String(value));
    }
    for (const [key, value] of Object.entries(attrs ?? {})) el.setAttribute(key, value);
    if (text !== undefined) el.textContent = text;
    if (on) {
      for (const [name, handler] of Object.entries(on)) {
        el.addEventListener(name, handler as EventListener);
      }
    }
  }
  append(el, children);
  return el;
}

const SVG_NS = "http://www.w3.org/2000/svg";

/** Create an SVG element. Attributes only: inline `style` is blocked by the CSP. */
export function svg<K extends keyof SVGElementTagNameMap>(
  tag: K,
  attrs?: Record<string, string | number> | null,
  ...children: Child[]
): SVGElementTagNameMap[K] {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [key, value] of Object.entries(attrs ?? {})) el.setAttribute(key, String(value));
  append(el, children);
  return el;
}

/** Look up a required element by id, with the type the caller expects. */
export function $<T extends HTMLElement = HTMLElement>(id: string): T {
  const el = document.getElementById(id);
  if (!el) throw new Error(`Missing #${id} in the page`);
  return el as T;
}

/** Replace everything inside `parent`. */
export function replace(parent: Element, ...children: Child[]): void {
  parent.replaceChildren();
  append(parent, children);
}

export const fmt = (value: number, digits = 2): string =>
  Number.isFinite(value) ? value.toFixed(digits) : "–";

export const int = (value: number): string => Math.round(value).toLocaleString("en-US");
