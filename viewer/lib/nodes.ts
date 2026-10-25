export type NodeFile = {
  theta: string | null;
  epsilon: string | null;
  metric: string | null;
  users: number[];
  items: number[];
};

const USER_SIDES = new Set(["user", "u", "left"]);
const ITEM_SIDES = new Set(["item", "v", "right", "product"]);

function parseNode(raw: string, line: string): number {
  if (!/^[1-9]\d*$/.test(raw)) {
    throw new Error(`Node numbers are positive integers. Could not read "${line}".`);
  }
  return Number(raw);
}

/**
 * Read the heuristic dump. `nodes.csv` rows are `user,<id>` and `item,<id>`.
 * A one-column `left.csv` (`u`) or `right.csv` (`v`) is accepted too.
 */
export function parseNodes(text: string): NodeFile {
  const users = new Set<number>();
  const items = new Set<number>();
  let theta: string | null = null;
  let epsilon: string | null = null;
  let metric: string | null = null;
  let single: "user" | "item" | null = null;

  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith("#")) continue;
    const parts = line.split(",").map((part) => part.trim());
    const head = parts[0]?.toLowerCase() ?? "";
    const tail = parts[1] ?? "";

    if (parts.length === 1) {
      if (head === "u") {
        single = "user";
        continue;
      }
      if (head === "v") {
        single = "item";
        continue;
      }
      const node = parseNode(head, line);
      if (single === "item") items.add(node);
      else if (single === "user") users.add(node);
      else throw new Error(`Expected side,node rows. Could not read "${line}".`);
      continue;
    }

    if (head === "theta") {
      theta = tail;
      continue;
    }
    if (head === "epsilon") {
      epsilon = tail;
      continue;
    }
    if (head === "metric") {
      metric = tail;
      continue;
    }
    if (head === "side" && tail.toLowerCase() === "node") continue;

    if (USER_SIDES.has(head)) {
      users.add(parseNode(tail, line));
      continue;
    }
    if (ITEM_SIDES.has(head)) {
      items.add(parseNode(tail, line));
      continue;
    }
    throw new Error(`Expected a user or item node. Could not read "${line}".`);
  }

  if (users.size === 0 && items.size === 0) {
    throw new Error("No node numbers found. Upload the nodes.csv written by bin/heuristic.jl.");
  }

  return {
    theta,
    epsilon,
    metric,
    users: [...users].sort((a, b) => a - b),
    items: [...items].sort((a, b) => a - b),
  };
}
