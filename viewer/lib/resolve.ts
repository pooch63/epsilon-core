import type { NodeFile } from "./nodes";

export type ReviewEdge = {
  userNode: number;
  itemNode: number;
  userId: string;
  asin: string;
  rating: number | null;
  title: string | null;
  verified: boolean | null;
};

export type ResolvedCore = {
  theta: string | null;
  epsilon: string | null;
  metric: string | null;
  users: { node: number; userId: string | null; edges: ReviewEdge[] }[];
  items: { node: number; asin: string | null; userCount: number }[];
  edges: ReviewEdge[];
  missingUsers: number[];
  missingItems: number[];
  reviews: number;
  skipped: number;
  distinctUsers: number;
  distinctItems: number;
};

function asId(value: unknown): string | null {
  if (typeof value === "string") {
    const trimmed = value.trim();
    return trimmed.length > 0 ? trimmed : null;
  }
  if (typeof value === "number" && Number.isFinite(value)) return String(value);
  return null;
}

/**
 * Walk a reviews JSONL in file order and number users and ASINs by first
 * appearance. That is the same numbering `index_dataset` writes into the
 * graph, so heuristic node `n` on the item side is the nth distinct ASIN.
 * An edge is a review whose user and ASIN are both in the dumped node sets.
 */
export function resolveCore(nodes: NodeFile, jsonl: string): ResolvedCore {
  const wantUsers = new Set(nodes.users);
  const wantItems = new Set(nodes.items);
  const userIndex = new Map<string, number>();
  const itemIndex = new Map<string, number>();
  const userIdAt = new Map<number, string>();
  const asinAt = new Map<number, string>();
  const edges: ReviewEdge[] = [];
  let reviews = 0;
  let skipped = 0;

  for (const raw of jsonl.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    let row: unknown;
    try {
      row = JSON.parse(line);
    } catch {
      skipped += 1;
      continue;
    }
    if (row === null || typeof row !== "object") {
      skipped += 1;
      continue;
    }
    const record = row as Record<string, unknown>;
    const userId = asId(record.user_id);
    const asin = asId(record.asin);
    if (userId === null || asin === null) {
      skipped += 1;
      continue;
    }
    reviews += 1;

    let userNode = userIndex.get(userId);
    if (userNode === undefined) {
      userNode = userIndex.size + 1;
      userIndex.set(userId, userNode);
      userIdAt.set(userNode, userId);
    }
    let itemNode = itemIndex.get(asin);
    if (itemNode === undefined) {
      itemNode = itemIndex.size + 1;
      itemIndex.set(asin, itemNode);
      asinAt.set(itemNode, asin);
    }

    if (!wantUsers.has(userNode) || !wantItems.has(itemNode)) continue;
    const rating = typeof record.rating === "number" ? record.rating : null;
    const title = typeof record.title === "string" && record.title.trim() ? record.title.trim() : null;
    const verified = typeof record.verified_purchase === "boolean" ? record.verified_purchase : null;
    edges.push({ userNode, itemNode, userId, asin, rating, title, verified });
  }

  const edgesByUser = new Map<number, ReviewEdge[]>();
  for (const edge of edges) {
    const list = edgesByUser.get(edge.userNode);
    if (list) list.push(edge);
    else edgesByUser.set(edge.userNode, [edge]);
  }

  const users = nodes.users.map((node) => ({
    node,
    userId: userIdAt.get(node) ?? null,
    edges: edgesByUser.get(node) ?? [],
  }));

  const usersByItem = new Map<number, Set<number>>();
  for (const edge of edges) {
    let buyers = usersByItem.get(edge.itemNode);
    if (!buyers) {
      buyers = new Set();
      usersByItem.set(edge.itemNode, buyers);
    }
    buyers.add(edge.userNode);
  }

  const items = nodes.items.map((node) => ({
    node,
    asin: asinAt.get(node) ?? null,
    userCount: usersByItem.get(node)?.size ?? 0,
  }));

  return {
    theta: nodes.theta,
    epsilon: nodes.epsilon,
    metric: nodes.metric,
    users,
    items,
    edges,
    missingUsers: nodes.users.filter((node) => !userIdAt.has(node)),
    missingItems: nodes.items.filter((node) => !asinAt.has(node)),
    reviews,
    skipped,
    distinctUsers: userIndex.size,
    distinctItems: itemIndex.size,
  };
}
