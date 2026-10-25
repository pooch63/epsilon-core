"use client";

import { useState } from "react";
import type { ResolvedCore, ReviewEdge } from "@/lib/resolve";

const NODE_R = 13;
const ROW = 56;
const TOP = 64;
const LEFT_X = 110;
const RIGHT_X = 560;
const WIDTH = 960;

type Selection =
  | { kind: "user"; node: number }
  | { kind: "item"; node: number }
  | { kind: "edge"; userNode: number; itemNode: number };

function amazonUrl(asin: string) {
  return `https://www.amazon.com/dp/${encodeURIComponent(asin)}`;
}

function edgeHot(edge: ReviewEdge, selected: Selection | null) {
  if (!selected) return false;
  if (selected.kind === "user") return edge.userNode === selected.node;
  if (selected.kind === "item") return edge.itemNode === selected.node;
  return edge.userNode === selected.userNode && edge.itemNode === selected.itemNode;
}

export default function CoreCanvas({ result }: { result: ResolvedCore }) {
  const [selected, setSelected] = useState<Selection | null>(null);
  const rows = Math.max(result.users.length, result.items.length, 1);
  const height = TOP + (rows - 1) * ROW + 48;

  const userAt = new Map(result.users.map((user, index) => [user.node, { ...user, y: TOP + index * ROW }]));
  const itemAt = new Map(result.items.map((item, index) => [item.node, { ...item, y: TOP + index * ROW }]));

  const drawn = new Map<string, ReviewEdge[]>();
  for (const edge of result.edges) {
    const key = `${edge.userNode}-${edge.itemNode}`;
    const list = drawn.get(key);
    if (list) list.push(edge);
    else drawn.set(key, [edge]);
  }

  const hotEdges = selected ? result.edges.filter((edge) => edgeHot(edge, selected)) : [];

  function faded(side: "user" | "item", node: number) {
    if (!selected) return false;
    if (selected.kind === "user") {
      if (side === "user") return node !== selected.node;
      return !hotEdges.some((edge) => edge.itemNode === node);
    }
    if (selected.kind === "item") {
      if (side === "item") return node !== selected.node;
      return !hotEdges.some((edge) => edge.userNode === node);
    }
    return side === "user" ? node !== selected.userNode : node !== selected.itemNode;
  }

  function toggle(next: Selection) {
    setSelected((current) => {
      if (!current) return next;
      if (current.kind !== next.kind) return next;
      if (current.kind === "edge" && next.kind === "edge") {
        return current.userNode === next.userNode && current.itemNode === next.itemNode ? null : next;
      }
      if (current.kind !== "edge" && next.kind !== "edge" && current.node === next.node) return null;
      return next;
    });
  }

  return (
    <div className="graph-canvas">
      <div className="graph-legend">
        <span><i className="swatch user" /> Users</span>
        <span><i className="swatch item" /> Products</span>
        <span><i className="swatch edge" /> Review edge</span>
      </div>
      <div className="canvas-scroll">
        <svg
          width={WIDTH}
          height={height}
          viewBox={`0 0 ${WIDTH} ${height}`}
          onClick={() => setSelected(null)}
        >
          <defs>
            <pattern id="core-grid" width="48" height="48" patternUnits="userSpaceOnUse">
              <circle cx="1" cy="1" r="1.15" className="grid-dot" />
            </pattern>
          </defs>
          <rect className="canvas-bg" x={0} y={0} width={WIDTH} height={height} fill="url(#core-grid)" />
          <text x={LEFT_X} y={28} textAnchor="middle" className="col-label">USERS</text>
          <text x={RIGHT_X} y={28} textAnchor="middle" className="col-label">PRODUCTS</text>
          {[...drawn.entries()].map(([key, reviews]) => {
            const edge = reviews[0];
            const user = userAt.get(edge.userNode);
            const item = itemAt.get(edge.itemNode);
            if (!user || !item) return null;
            const hot = edgeHot(edge, selected);
            const faded = selected !== null && !hot;
            return (
              <g key={key} onClick={(event) => { event.stopPropagation(); toggle({ kind: "edge", userNode: edge.userNode, itemNode: edge.itemNode }); }}>
                <title>{`${edge.userNode} — ${edge.asin}`}</title>
                <line
                  x1={LEFT_X}
                  y1={user.y}
                  x2={RIGHT_X}
                  y2={item.y}
                  stroke="transparent"
                  strokeWidth={14}
                  style={{ cursor: "pointer" }}
                />
                <line
                  x1={LEFT_X}
                  y1={user.y}
                  x2={RIGHT_X}
                  y2={item.y}
                  stroke={hot ? "var(--solution)" : "var(--edge)"}
                  strokeWidth={hot ? 2.4 : 1.4}
                  opacity={faded ? 0.18 : 1}
                  pointerEvents="none"
                />
              </g>
            );
          })}
          {[...userAt.values()].map((user) => (
            <Node
              key={`u${user.node}`}
              x={LEFT_X}
              y={user.y}
              label={String(user.node)}
              side="user"
              selected={selected?.kind === "user" && selected.node === user.node}
              faded={faded("user", user.node)}
              onClick={() => toggle({ kind: "user", node: user.node })}
            />
          ))}
          {[...itemAt.values()].map((item) => (
            <Node
              key={`v${item.node}`}
              x={RIGHT_X}
              y={item.y}
              label={String(item.node)}
              side="item"
              caption={item.asin ?? "missing"}
              selected={selected?.kind === "item" && selected.node === item.node}
              faded={faded("item", item.node)}
              onClick={() => toggle({ kind: "item", node: item.node })}
            />
          ))}
        </svg>
      </div>
      <SelectionDetail result={result} selected={selected} edges={hotEdges} />
    </div>
  );
}

function Node({
  x,
  y,
  label,
  side,
  caption,
  selected,
  faded,
  onClick,
}: {
  x: number;
  y: number;
  label: string;
  side: "user" | "item";
  caption?: string;
  selected: boolean;
  faded: boolean;
  onClick: () => void;
}) {
  return (
    <g
      transform={`translate(${x},${y})`}
      onClick={(event) => {
        event.stopPropagation();
        onClick();
      }}
      style={{ cursor: "pointer" }}
      opacity={faded ? 0.28 : 1}
    >
      <circle
        r={NODE_R}
        fill={side === "user" ? "var(--u-color)" : "var(--v-color)"}
        stroke={selected ? "#1a1d24" : "none"}
        strokeWidth={selected ? 2.5 : 0}
      />
      <text x={0} y={3.5} textAnchor="middle" className="node-label" fill="#ffffff" fontWeight={700} fontSize={10}>
        {label}
      </text>
      {caption ? (
        <text x={NODE_R + 10} y={4} className="node-label">
          {caption}
        </text>
      ) : null}
    </g>
  );
}

function SelectionDetail({
  result,
  selected,
  edges,
}: {
  result: ResolvedCore;
  selected: Selection | null;
  edges: ReviewEdge[];
}) {
  if (!selected) {
    return <p className="canvas-note">Click a node or an edge. Users are on the left, products on the right.</p>;
  }
  const title = selected.kind === "user"
    ? `User ${selected.node}`
    : selected.kind === "item"
      ? `Product ${selected.node}`
      : `User ${selected.userNode} — product ${selected.itemNode}`;
  const item = selected.kind === "item" ? result.items.find((entry) => entry.node === selected.node) : null;
  return (
    <div className="canvas-detail">
      <p className="canvas-note">
        <b>{title}</b>
        {item?.asin ? (
          <>
            {" · "}
            <a href={amazonUrl(item.asin)} target="_blank" rel="noreferrer">{item.asin}</a>
          </>
        ) : null}
        {" · "}
        {edges.length} {edges.length === 1 ? "edge" : "edges"}
      </p>
      {edges.length === 0 ? <p className="canvas-note">No review connects this node inside the core.</p> : null}
      {edges.map((edge, index) => (
        <p className="canvas-note" key={`${edge.userNode}-${edge.itemNode}-${index}`}>
          {edge.userNode} → {edge.itemNode} {edge.asin}
          {edge.rating != null ? ` · ${edge.rating}` : ""}
          {edge.title ? ` · ${edge.title}` : ""}
          {" · "}
          <a href={amazonUrl(edge.asin)} target="_blank" rel="noreferrer">Amazon</a>
        </p>
      ))}
    </div>
  );
}
