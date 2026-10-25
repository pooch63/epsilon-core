import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { parseNodes } from "../lib/nodes.ts";
import { resolveCore } from "../lib/resolve.ts";

const root = resolve(import.meta.dirname, "../..");
const nodes = parseNodes(
  readFileSync(resolve(root, "results/amazon/boxes/theta-2-epsilon-0.5-overlap/nodes.csv"), "utf8"),
);
const jsonl = readFileSync(resolve(root, "data/amazon/boxes/raw/Subscription_Boxes.jsonl"), "utf8");
const core = resolveCore(nodes, jsonl);

const asins = core.items.map((item) => item.asin);
const expected = ["B097B34JVF", "B097B2N157", "B09M7GDZ46", "B09M7BQDHG"];
if (JSON.stringify(asins) !== JSON.stringify(expected)) {
  throw new Error(`ASINs ${asins.join(",")} !== ${expected.join(",")}`);
}
if (core.edges.length !== 15) throw new Error(`expected 15 edges, got ${core.edges.length}`);
if (core.missingUsers.length || core.missingItems.length) throw new Error("missing nodes");
if (core.distinctUsers !== 15237 || core.distinctItems !== 1128) {
  throw new Error(`counts ${core.distinctUsers} ${core.distinctItems}`);
}
const user630 = core.users.find((user) => user.node === 630);
const items630 = user630?.edges.map((edge) => edge.asin).join(",");
if (items630 !== "B097B34JVF,B097B2N157") throw new Error(`user 630 products ${items630}`);
console.log(
  `ok users=${core.users.length} products=${core.items.length} edges=${core.edges.length} asins=${asins.join(" ")}`,
);
