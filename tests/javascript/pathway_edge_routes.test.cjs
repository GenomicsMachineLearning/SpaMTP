// Run with node tests/javascript/pathway_edge_routes.test.cjs [template] [links.json].
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const template = process.argv[2] || path.join(__dirname, "../../inst/templates/pathway_network.html");
const html = fs.readFileSync(template, "utf8");
const functions = ["assignEdgeRoutes", "edgePath"].map(name => {
  const match = html.match(new RegExp(`    function ${name}\\([^]*?\\n    \\}`));
  assert.ok(match, `Missing ${name} in the production template`);
  return match[0];
}).join("\n");
const context = vm.createContext({});
vm.runInContext(functions, context);
const {assignEdgeRoutes, edgePath} = context;
const node = (id, x, y) => ({id, x, y});
const a = node("A", 10, 20), b = node("B", 210, 120);
const edge = (source, target, reaction) => ({source, target, reaction, head: "arrow"});
const test = (name, body) => { body(); console.log(`PASS ${name}`); };
const midpoint = link => {
  const parts = edgePath(link).match(/-?\d+(?:\.\d+)?(?:e[+-]?\d+)?/gi).map(Number);
  assert.equal(parts.length, 6);
  return [(parts[0] + 2 * parts[2] + parts[4]) / 4,
    (parts[1] + 2 * parts[3] + parts[5]) / 4].join(",");
};

test("all parallel interactions occupy distinct curves", () => {
  const links = ["activation", "inhibition", "phosphorylation"].map(r => edge(a, b, r));
  assignEdgeRoutes(links);
  assert.equal(new Set(links.map(midpoint)).size, links.length);
});
test("opposite directions cannot trace the same geometric curve", () => {
  const links = [edge(a, b, "activation"), edge(b, a, "inhibition"), edge(a, b, "binding")];
  assignEdgeRoutes(links);
  assert.equal(new Set(links.map(midpoint)).size, links.length);
});
test("multiple self-loops remain separate", () => {
  const links = [edge(a, a, "activation"), edge(a, a, "inhibition")];
  assignEdgeRoutes(links);
  assert.equal(new Set(links.map(edgePath)).size, 2);
});
test("reordering the input preserves each interaction's route", () => {
  const links = [edge(a, b, "activation"), edge(b, a, "inhibition"), edge(a, b, "binding")];
  assignEdgeRoutes(links);
  const original = new Map(links.map(e => [e.reaction, edgePath(e)]));
  links.reverse();
  assignEdgeRoutes(links);
  links.forEach(e => assert.equal(edgePath(e), original.get(e.reaction)));
});
test("routing works before and after D3 resolves endpoint IDs", () => {
  const links = [edge("A", "B", "activation"), edge("B", "A", "inhibition")];
  assignEdgeRoutes(links);
  const indices = links.map(e => e.routeIndex);
  links.forEach(e => { e.source = e.source === "A" ? a : b; e.target = e.target === "A" ? a : b; });
  assignEdgeRoutes(links);
  assert.deepEqual(links.map(e => e.routeIndex), indices);
});
test("node identifiers cannot collide through a grouping delimiter", () => {
  const links = [edge("A|B", "C", "activation"), edge("A", "B|C", "inhibition")];
  assignEdgeRoutes(links);
  links.forEach(e => assert.equal(e.routeCount, 1));
});
test("coincident endpoints produce finite paths", () => {
  const links = [edge(a, node("B", a.x, a.y), "activation"), edge(a, a, "inhibition")];
  assignEdgeRoutes(links);
  links.forEach(e => assert.doesNotMatch(edgePath(e), /NaN|Infinity/));
});
test("a single non-loop interaction retains its original arc", () => {
  const link = edge(a, b, "activation");
  assignEdgeRoutes([link]);
  assert.match(edgePath(link), /^M10,20A/);
});

if (process.argv[3]) test("every parallel pair in the supplied graph has distinct geometry", () => {
  const links = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
  const ids = [...new Set(links.flatMap(e => [e.source, e.target]))].sort();
  const nodes = new Map(ids.map((id, i) => [id, node(id, i * 31, (i % 7) * 43)]));
  assignEdgeRoutes(links);
  const groups = new Map();
  links.forEach(e => {
    const key = JSON.stringify([e.source, e.target].sort());
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(e);
    e.source = nodes.get(e.source); e.target = nodes.get(e.target);
  });
  groups.forEach(group => {
    if (group.length < 2) return;
    const paths = group.map(e => e.source.id === e.target.id ? edgePath(e) : midpoint(e));
    assert.equal(new Set(paths).size, group.length, `Overlapping interactions: ${group[0].source.id}`);
  });
  console.log(`Checked ${links.length} edges across ${groups.size} unordered node pairs`);
});
