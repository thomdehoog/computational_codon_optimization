import type { Job, Meta } from "./types";

/** An error the server explained, safe to show to the user as written. */
export class ApiError extends Error {
  constructor(message: string, readonly status: number) {
    super(message);
    this.name = "ApiError";
  }
}

const TOKEN_HEADER = "X-Codonopt-Token";

/** Thin typed wrapper over the codonopt JSON API. */
export class Api {
  constructor(private readonly token: string) {}

  meta(): Promise<Meta> {
    return this.request<Meta>("GET", "/api/meta");
  }

  start(payload: Record<string, unknown>): Promise<Job> {
    return this.request<Job>("POST", "/api/jobs", payload);
  }

  poll(id: string): Promise<Job> {
    return this.request<Job>("GET", `/api/jobs/${encodeURIComponent(id)}`);
  }

  cancel(id: string): Promise<Job> {
    return this.request<Job>("POST", `/api/jobs/${encodeURIComponent(id)}/cancel`);
  }

  private async request<T>(method: string, path: string, body?: unknown): Promise<T> {
    let response: Response;
    try {
      response = await fetch(path, {
        method,
        headers: {
          [TOKEN_HEADER]: this.token,
          ...(body === undefined ? {} : { "Content-Type": "application/json" }),
        },
        body: body === undefined ? undefined : JSON.stringify(body),
        cache: "no-store",
      });
    } catch {
      throw new ApiError(
        "Cannot reach the codonopt server. Was the window or terminal that started it closed?",
        0,
      );
    }
    let data: unknown = null;
    try {
      data = await response.json();
    } catch {
      /* fall through to the status check */
    }
    if (!response.ok) {
      const message =
        data && typeof data === "object" && "error" in data
          ? String((data as { error: unknown }).error)
          : `The server answered ${response.status}.`;
      throw new ApiError(message, response.status);
    }
    return data as T;
  }
}
