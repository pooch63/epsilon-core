"use client";

import { useState } from "react";
import CoreCanvas from "./core-canvas";
import { parseNodes } from "@/lib/nodes";
import { resolveCore, type ResolvedCore } from "@/lib/resolve";

function amazonUrl(asin: string) {
  return `https://www.amazon.com/dp/${encodeURIComponent(asin)}`;
}

function formatBytes(bytes: number) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(0)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

export default function Viewer() {
  const [nodesName, setNodesName] = useState<string | null>(null);
  const [jsonlName, setJsonlName] = useState<string | null>(null);
  const [status, setStatus] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<ResolvedCore | null>(null);

  async function onSubmit(event: React.FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = event.currentTarget;
    const nodeFile = (form.elements.namedItem("nodes") as HTMLInputElement).files?.[0];
    const jsonlFile = (form.elements.namedItem("jsonl") as HTMLInputElement).files?.[0];
    if (!nodeFile || !jsonlFile) {
      setError("Choose the node file and the raw JSONL.");
      setResult(null);
      return;
    }
    setError(null);
    setResult(null);
    setStatus(`Reading ${formatBytes(nodeFile.size + jsonlFile.size)}…`);
    try {
      const [nodeText, jsonlText] = await Promise.all([nodeFile.text(), jsonlFile.text()]);
      setStatus("Numbering users and ASINs in file order…");
      const resolved = resolveCore(parseNodes(nodeText), jsonlText);
      setResult(resolved);
      setStatus(null);
    } catch (caught) {
      setStatus(null);
      setError(caught instanceof Error ? caught.message : "Could not read those files.");
    }
  }

  return (
    <main className="page">
      <header>
        <h1>Epsilon-core viewer</h1>
        <p>
          Upload <code>nodes.csv</code> from the heuristic and the raw Amazon reviews JSONL.
          Item node <i>n</i> is the <i>n</i>th distinct ASIN in that file, and user node <i>n</i> is
          the <i>n</i>th distinct reviewer. A row appears when a review connects one of those users
          to one of those products.
        </p>
      </header>

      <form onSubmit={onSubmit}>
        <label className="file">
          <span>Node file</span>
          <input
            name="nodes"
            type="file"
            accept=".csv,text/csv,text/plain"
            onChange={(event) => {
              setNodesName(event.target.files?.[0]?.name ?? null);
              setResult(null);
            }}
          />
          <small>{nodesName ?? "nodes.csv"}</small>
        </label>
        <label className="file">
          <span>Raw reviews JSONL</span>
          <input
            name="jsonl"
            type="file"
            accept=".jsonl,.json,application/jsonl,application/json,text/plain"
            onChange={(event) => {
              setJsonlName(event.target.files?.[0]?.name ?? null);
              setResult(null);
            }}
          />
          <small>{jsonlName ?? "Subscription_Boxes.jsonl"}</small>
        </label>
        <button type="submit" className="submit">
          Look up products
        </button>
      </form>

      {status ? <p className="status">{status}</p> : null}
      {error ? <p className="error">{error}</p> : null}

      {result ? <Report result={result} /> : null}
    </main>
  );
}

function Report({ result }: { result: ResolvedCore }) {
  const settings = [
    result.theta ? `θ = ${result.theta}` : null,
    result.epsilon ? `ε = ${result.epsilon}` : null,
    result.metric,
  ].filter(Boolean);

  return (
    <>
      <div className="stats">
        <b>{result.users.length} users</b>
        <b>{result.items.length} products</b>
        <b>{result.edges.length} edges</b>
      </div>
      <p className="meta">
        {settings.length > 0 ? `${settings.join(" · ")}. ` : null}
        Scanned {result.reviews.toLocaleString()} reviews
        {" · "}
        {result.distinctUsers.toLocaleString()} users
        {" · "}
        {result.distinctItems.toLocaleString()} ASINs
        {result.skipped > 0 ? ` · skipped ${result.skipped.toLocaleString()} malformed rows` : ""}.
      </p>
      {result.missingUsers.length > 0 || result.missingItems.length > 0 ? (
        <p className="error">
          Some node numbers never appeared in this JSONL
          {result.missingUsers.length > 0 ? `. Users: ${result.missingUsers.join(", ")}` : ""}
          {result.missingItems.length > 0 ? `. Items: ${result.missingItems.join(", ")}` : ""}.
          The JSONL has to be the same file that was indexed into the graph.
        </p>
      ) : null}

      <section>
        <h2>Edges</h2>
        <CoreCanvas result={result} />
      </section>

      <section>
        <h2>Products in the group</h2>
        <table className="products-table">
          <thead>
            <tr>
              <th>Node</th>
              <th>ASIN</th>
              <th>Core users</th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            {result.items.map((item) => (
              <tr key={item.node}>
                <td>{item.node}</td>
                <td>{item.asin ?? "Not in this file"}</td>
                <td>{item.userCount}</td>
                <td>
                  {item.asin ? (
                    <a href={amazonUrl(item.asin)} target="_blank" rel="noreferrer">
                      View on Amazon
                    </a>
                  ) : null}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </section>

      <section>
        <h2>Products mapped to each user</h2>
        <div className="users">
          {result.users.map((user) => (
            <article className="card" key={user.node}>
              <h3>User {user.node}</h3>
              <p className="userid">{user.userId ?? "Not in this file"}</p>
              {user.edges.length === 0 ? (
                <p className="empty">No review in this JSONL connects this user to a product in the group.</p>
              ) : (
                <table>
                  <thead>
                    <tr>
                      <th>Item</th>
                      <th>ASIN</th>
                      <th>Rating</th>
                      <th>Review title</th>
                      <th></th>
                    </tr>
                  </thead>
                  <tbody>
                    {user.edges.map((edge, index) => (
                      <tr key={`${edge.itemNode}-${index}`}>
                        <td>{edge.itemNode}</td>
                        <td>{edge.asin}</td>
                        <td>{edge.rating ?? "—"}{edge.verified ? " · verified" : ""}</td>
                        <td>{edge.title ?? "—"}</td>
                        <td>
                          <a href={amazonUrl(edge.asin)} target="_blank" rel="noreferrer">
                            View
                          </a>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
            </article>
          ))}
        </div>
      </section>

    </>
  );
}
