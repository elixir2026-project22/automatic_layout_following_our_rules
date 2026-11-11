# Chat transcript — SBML pathway layouts

Conversation between the user and Claude about laying out metabolic pathways from SBML files.

**What this file contains:** every user message and every reply from Claude, in order, as text.
**What it does not contain:** the internal tool calls (code runs, file reads, web fetches) and their raw outputs, and the images that were shown inline. All files produced during the conversation are in the same outputs folder; the reply that produced each file names it.

---

## 1. User

How many graph layouts do you know?

## 1. Claude

Quite a lot. If you mean graph drawing (laying out nodes and edges), I know several dozen algorithms, falling into about eight families.

**Force-directed** layouts treat the graph as a physical system: Eades' spring embedder, Fruchterman–Reingold, Kamada–Kawai (based on graph-theoretic distances), ForceAtlas2 (Gephi's default), Yifan Hu's multilevel method, OpenOrd/DrL for very large graphs, the stress majorization used in Graphviz's neato, and fCoSE/CoSE in Cytoscape.

**Hierarchical/layered** layouts suit directed graphs and flowcharts. The Sugiyama framework (layer assignment, crossing minimization, coordinate assignment) is the foundation, as in Graphviz's dot, Dagre, and ELK Layered.

**Tree** layouts include Reingold–Tilford tidy trees, Walker's algorithm, radial trees, balloon layouts, dendrograms, H-trees for compact binary trees, and space-filling variants like treemaps, sunbursts, and icicle plots.

**Circular** layouts place nodes on one or more circles (Graphviz's circo), often with crossing reduction or hierarchical edge bundling. Chord diagrams belong here.

**Spectral and dimensionality-reduction** layouts use eigenvectors of the Laplacian or adjacency matrix, classical MDS, pivot MDS, or embeddings like t-SNE/UMAP applied to node vectors.

**Orthogonal** layouts, used for circuit and UML-style diagrams, draw edges with horizontal and vertical segments only. Tamassia's topology-shape-metrics approach is the classic one.

**Planar** layouts include straight-line embeddings such as de Fraysseix–Pach–Pollack and Schnyder's algorithm, plus Tutte's barycentric embedding for 3-connected planar graphs.

**Other specialized** layouts include arc diagrams, hive plots, BioFabric (nodes as lines), adjacency matrix views, geographic layouts, grid and random layouts as baselines, bipartite and two-column layouts, Sankey diagrams for flows, timeline/storyline layouts for dynamic graphs, and constraint-based layouts like WebCoLa.

Different tools also have their own named implementations of these (NetworkX, igraph, Graphviz's six engines, Cytoscape, D3, ELK), so the count depends on whether you count families or variants. Are you choosing a layout for a particular graph? I can help pick one if you tell me its size and structure.

---

## 2. User

*(uploaded: brenda_citric-acid-cycle.xml)*

Given this model, could you test different layouts on it to visualize the pathway?

## 2. Claude

I tested 12 layouts on your model. The file is an SBML Level 3 model from BRENDA with 34 species and 24 reactions, all marked reversible. I drew it as a bipartite graph, with metabolites as circles and reactions as small squares labelled by enzyme name from the model's EC groups. The 12 tests are in `tca_layout_comparison.png`.

The choice that mattered most wasn't the layout algorithm but how cofactors are handled. NAD+/NADH, CoA, ATP/ADP, H₂O, CO₂ and phosphate each take part in many reactions. If each is drawn as a single shared node, it pulls everything toward it, and every algorithm produced a tangle (panels A–D). Spectral (C) collapsed almost completely, and circular (D) was unreadable. Drawing a separate copy of each cofactor for every reaction (panels E–L) is standard practice for metabolic maps, and it let the cycle appear on its own in most layouts.

With cofactors copied, the general-purpose algorithms ranked roughly like this:
- **Best:** Graphviz `sfdp` (G) and `neato` (F) recovered the ring most clearly. Kamada–Kawai (E) was close behind.
- **Middling:** Radial `twopi` (I) and spectral (K) found the loop but distorted it.
- **Poor:** Layered `dot` (H) is a bad fit because a cycle has no top or bottom, so it has to break the loop and route edges backward. `circo` (J) piled the copied cofactors onto one point.

**Recommendation:** None of the automatic layouts places the intermediates in pathway order. The one I'd use is a hybrid I built for this model, in `tca_cycle_anchored.png`. The nine core intermediates sit on a circle in textbook order, with acetyl-CoA entering at the top. Side branches sit outside the ring: maleate, the citryl-CoA route, and acetoacetate. The glyoxylate shunt (isocitrate lyase and malate synthase) crosses the middle. Alternative enzymes for the same step, such as NAD vs. NADP isocitrate dehydrogenase or ADP- vs. GDP-forming succinyl-CoA synthetase, are spread side by side so they don't overlap.

To make the arrows follow the usual cycle direction, I reversed how four reactions are written: aconitase (isocitrate side), fumarase, and both succinyl-CoA synthetases. This changes only the drawing, since the model already marks every reaction as reversible.

`tca_layouts.py` regenerates everything from the SBML file. It uses networkx, matplotlib and the Graphviz command-line tools.

---

## 3. User

I have a set of rules to apply for the layouts:

* The rendering must be reproducible
* Labels must have fixed positions, according to their type
* No line crossing
* Avoid long lines
* Avoid line bends

## 3. Claude

I don't have a finished layout that meets all five rules yet. Below is what I found, how far the new layout got, and what's left. Your earlier files (`tca_layout_comparison.png`, `tca_cycle_anchored.png`, `tca_layouts.py`) are unchanged.

**None of the earlier layouts pass "no line crossing."** I counted crossings with cofactors drawn as separate copies per reaction: sfdp had 21, neato 43, and my cycle-anchored layout 45.

**Crossings can't be removed by layout alone.** The pathway network is not planar, meaning no arrangement of its nodes draws it without crossings, whether cofactors are shared or copied. Removing any one or two reactions doesn't fix that.

**Drawing acetyl-CoA three times does fix it.** This is the usual clone-marker approach in pathway maps, where one metabolite is shown in several places. A search over the ways to split acetyl-CoA showed:
- No split into two copies works; at least three are needed.
- A sensible split is: one copy for citrate synthase and ATP citrate synthase, one for malate synthase (the glyoxylate shunt), and one for the citryl-CoA and acetate CoA-transferase reactions.
- Duplicating oxaloacetate or citrate instead would also work, but acetyl-CoA is the most natural choice because it is an input coming from outside the cycle.

**How the new layout handles each rule:**
- **Reproducible:** The layout uses no randomness. It still varied between runs at first, because networkx's planarity test depends on Python's per-run string hashing (`PYTHONHASHSEED`). The script now fixes that setting, and the starting layout is identical across runs.
- **Fixed label positions by type:** Metabolite labels sit centred above the node, reaction labels centred below, and cofactor labels to the right. The layout moves nodes to clear overlaps instead of shifting labels.
- **No crossings:** The layout starts from a crossing-free arrangement. The optimiser never accepts a move that creates a crossing, checking both the end position and points along the way.
- **No bends:** All edges are straight segments.
- **Avoid long lines:** The optimiser penalises uneven edge lengths, with an extra penalty for edges more than 1.5 times the target length.

**Progress so far:** The first version had zero crossings and no bends. But edges ranged from 0.25× to 2.7× the target length, and about 50 label positions still overlapped a line. It used a slow optimiser that moved one node at a time and got stuck. I replaced it with one that moves all nodes together and holds back any node whose move would create a crossing. In a 30-step test, its score fell from about 9,400 to 3,400, versus about 1,050 for the old optimiser after a full 12-minute run. It is clearly still improving, but I haven't run it to completion.

**Still to do:**
1. Run the new optimiser to completion and finish with a short node-by-node clean-up.
2. Try each of the four candidate outer boundaries and keep the best.
3. Report the final crossing count, the label overlaps and the edge-length range.
4. Confirm that two separate runs produce byte-identical SVGs.

**Possible conflicts between rules:** No crossings and no long lines may pull against each other. The acetate shuttle (succinyl-CoA ↔ acetate ↔ citrate) has to go around part of the cycle to avoid crossing it, so one or two edges may stay longer than the rest. If that happens, adding one more acetyl-CoA or acetate copy would shorten them. Would you accept that, or should the number of duplicated metabolites stay as small as possible?

---

## 4. User

I have a set of updated rules to apply for the layouts:

* The rendering must be reproducible
* Labels must have fixed positions, according to their type
* No line crossing
* try to avoid long lines
* try to avoid line bends

## 4. Claude

I don't have a finished layout for the updated rules yet, and there are no new files this turn. Your earlier files (`tca_layout_comparison.png`, `tca_cycle_anchored.png`, `tca_layouts.py`) are unchanged.

**What the relaxed rules change.** Allowing some bends and some long lines opens up orthogonal drawing (horizontal and vertical lines only), which is how KEGG-style metabolic maps are drawn. Allowing bends doesn't remove the need to draw acetyl-CoA three times, though: crossings depend on how the network is connected, not on how lines are routed. With those three copies the network can be drawn with no crossings.

**Straight-line attempt.** I finished running the straight-line optimiser.
- It holds the hard rules: zero crossings, no bends, and identical output from run to run.
- It doesn't produce a good picture. Edges range from 0.2× to 3.4× the target length, and about 50–75 labels still sit on top of a line.
- The cause is that labels hang outside their nodes, so the optimiser has to push lines away from them one small move at a time, and it gets stuck.

**Orthogonal attempt.** I switched to OGDF, a C++ graph-drawing library that I installed and called from Python. Its planarization layout guarantees zero crossings for this network and reduces bends as far as it can.
- Each metabolite, reaction and cofactor is drawn as a box with its name inside. Label position is then fixed by type: bold in a red rounded box for cycle intermediates, a blue box for other metabolites, a dashed outline for the acetyl-CoA copies, italic for enzyme names, small grey for cofactors.
- Because a line can never pass through a box and boxes never overlap, labels can't collide with lines or each other.
- OGDF's random seed is fixed, so the output should be reproducible.
- The script tries six ways of arranging the cycle's regions (OGDF's embedders) and checks each result independently for crossings, lines through boxes, overlapping boxes, bend count, longest edge and area.

A small OGDF test drew correctly with zero crossings. The full model crashed on its first layout call, most likely because Python and OGDF both try to free the same embedder object. That needs a one-line change, but I hit the tool limit before I could make it.

**Next steps:**
1. Fix the crash and run all six embedder options.
2. Pick the best by the order of your rules: zero crossings and fixed labels as hard checks, then fewest long lines, then fewest bends.
3. Confirm two runs produce byte-identical SVGs.
4. Deliver the winning drawing, a comparison sheet and the script.

One choice is yours: should the labels sit inside the nodes (boxes, which guarantee no label collisions) or beside them (dots with labels outside, which look lighter but need more space)?

---

## 5. User

Continue

## 5. Claude

I have a layout that passes all your rules, but it isn't saved as a file yet. I hit the tool limit before exporting it, so your files from before (`tca_layout_comparison.png`, `tca_cycle_anchored.png`, `tca_layouts.py`) are still the only ones available.

**How the drawing is built.** I used OGDF, a C++ graph-drawing library, called from Python. Its layout guarantees no crossings whenever the network can be drawn without them, and it routes lines horizontally and vertically while keeping bends to a minimum. The crash from last turn is fixed: Python and OGDF were both trying to free the same object. I also fixed a drawing bug that left some arrowheads floating short of their box.

**Labels.** Every name sits inside its node's box, so labels can't touch lines or each other. The look is fixed by node type:
- Cycle intermediates: bold, red rounded box.
- Other metabolites: blue rounded box.
- The three acetyl-CoA copies: dashed blue box.
- Reactions: enzyme name in italics, with the cofactor change written on a small grey second line (e.g. "NAD+ → NADH + CO2").

**Cofactors.** Writing cofactors inside the reaction boxes, instead of drawing them as separate nodes, made the biggest difference. It cuts the network from 104 lines to 59 and the bends from 35 to 9–14.

**Search for the best arrangement.** I tried five OGDF arrangement strategies, two routing modes and 120 fixed node orderings, about 500 layouts in all. Every one had zero crossings, zero lines through boxes and zero overlapping boxes. The orderings are generated from fixed numbers, so any chosen layout can be reproduced exactly. The best so far (Simple strategy, ordering 31) has:

| Check | Result |
|---|---|
| Crossings | 0 |
| Lines through a box / overlapping boxes | 0 / 0 |
| Straight lines (no bend) | 48 of 59 |
| Bends in total | 14 |
| Lines longer than 250 pt | 3 (the citrate → aconitase step and the two routes around the cycle to oxaloacetate) |

In this layout the cycle reads clockwise from oxaloacetate. The glyoxylate shunt runs along the top, and the citryl-CoA and acetate shuttles sit in the middle.

**Not done yet:**
1. The final ranking with corrected line lengths. One of four search batches is complete; the others were stopped when the tool time ran out.
2. Proof that two separate runs produce byte-identical SVG files.
3. Saving the winning PNG and SVG, a version with cofactors drawn as separate nodes for comparison, and the cleaned-up script.

Ask me to continue and I'll finish those and deliver the files.

---

## 6. User

Continue

## 6. Claude

The finished layout, `tca_pathway.svg` / `.png`, meets all five rules. The script `tca_orthogonal_layout.py` regenerates it from your SBML file.

| Rule | Result |
|---|---|
| Reproducible | Three separate runs produced byte-identical SVG and PNG files. |
| Fixed label positions by type | Every name sits centred in its node's box, styled by type. No box overlaps another and no line passes through one, so labels can't collide with anything. |
| No line crossing | 0 crossings, rechecked geometrically rather than trusting the layout library. |
| Try to avoid long lines | 3 of 59 lines are longer than 250 pt. |
| Try to avoid bends | 52 of 59 lines are straight; 9 bends in total. |

**How to read it:**
- **Red boxes** are the cycle intermediates. The loop runs isocitrate → 2-oxoglutarate → succinyl-CoA → succinate → fumarate → malate → oxaloacetate → citrate → cis-aconitate.
- **Blue boxes** are the other metabolites. Acetyl-CoA appears three times, in dashed boxes; without those copies, the network can't be drawn with no crossings.
- **White boxes** are reactions, with the cofactor change on a small grey second line. Four reversible reactions are drawn in the usual cycle direction.

**How the winner was chosen:** I generated about 500 layouts by varying the layout settings and the order nodes are fed in (fixed numbers, so reproducible). All of them met the three hard rules. I ranked them by long lines, then bends, then total and longest line length.

**Cofactors as separate nodes:** `tca_pathway_cofactor_nodes.svg` / `.png` draws each cofactor as its own small node. It also has zero crossings but costs 104 lines, 33 bends and 8 long lines, so I'd only use it if cofactors have to be separate nodes.

**Remaining compromises:**
- The three long lines are citrate → aconitase → cis-aconitate, the isocitrate lyase connection, and the route from oxaloacetate to citrate synthase. Shortening them further would need another duplicated metabolite.
- Lines are horizontal and vertical only, so the cycle reads as a rectangle rather than a circle.
- The script needs `pip install ogdf-wheel cppyy networkx matplotlib numpy`. `--search N` reruns the search, for example on a different model.

---

## 7. User

Improve the rendering with the cycle drawn as a cycle

## 7. Claude

The cycle is now drawn as a circle, in `tca_cycle.svg` and `tca_cycle.png`, and it still meets all your rules. `tca_cycle_layout.py` regenerates it from the SBML file and needs only networkx, numpy and matplotlib, not OGDF.

| Rule | Result |
|---|---|
| Reproducible | Positions are fixed in the script and nothing is random. Two runs produced byte-identical SVG and PNG files. |
| Fixed label positions by type | Names are centred inside their boxes, with the same styles as before. |
| No line crossing | 0 crossings, 0 lines through a box, 0 overlapping boxes. The script refuses to write output if any check fails. |
| Try to avoid long lines | The longest line is 298 pt, versus 717 pt in the previous layout. The median is 89 pt. |
| Try to avoid bends | 0 bends. Lines are straight, or smooth arcs following the circle. |

**How it's arranged:**
- The nine intermediates sit on the circle in pathway order, with acetyl-CoA entering at the top. The main enzyme for each step sits on the circle between its two metabolites, and those connections are drawn as thick red arcs.
- Alternative enzymes for the same step sit just inside or just outside the circle beside the main one. Examples are IDH (NADP) next to IDH (NAD), and MDH (NADP) and MQO next to MDH (NAD).
- The glyoxylate shunt crosses the middle: isocitrate lyase → glyoxylate → malate synthase → malate.
- The citryl-CoA bypass sits outside at the top. Maleate is on the left, and the acetoacetate enzyme is inside near the bottom.

**One change to the model's drawing:** acetate now appears twice and acetyl-CoA four times, in dashed boxes. Before, the acetate CoA-transferase had to connect back to the citryl-CoA bypass with lines running halfway around the cycle. Giving it its own copies removes those long lines entirely.

**Remaining long line:** the isocitrate lyase → succinate chord (298 pt) is the longest line. It comes from the shunt crossing the circle, so it stays unless you'd rather route the shunt outside, which would make it longer.

**Trade-off:** the positions were tuned by hand for this model, with the checker confirming each change. Unlike the previous OGDF script, this one won't automatically lay out a different pathway. If you edit `POS` in the script, the checker will report any rule you break.

---

## 8. User

I need the script to be able to read any SBML file

## 8. Claude

The general SBML script works end to end and meets all five rules on your file, but the drawing it produces is far too large to use, so I haven't delivered it yet. Your earlier files are unchanged.

**What the script does now (`sbml_cycle_layout.py`):**
- **Reads any SBML file** (levels 1–3), using the names in the file. When a reaction's name is only an identifier, as in your file (BR29304…), it takes a readable name from the SBML groups, e.g. "isocitrate dehydrogenase (NAD+)". Same-named species in different compartments get the compartment added to their label.
- **Side compounds** (H2O, ATP, NAD(P)H, CoA, CO2, quinones, ferredoxins and similar) are recognised by name and written inside the reaction box. `--side` and `--keep` adjust the list.
- **Main cycle** is found automatically. Each step must turn a substrate into a product, with partners in multi-substrate reactions matched by name (succinyl-CoA → succinate, acetate → acetyl-CoA). Among the remaining cycles it picks the one whose steps have the most alternative enzymes. On your file this gives exactly the nine TCA intermediates. `--cycle "A,B,C"` overrides the choice.
- **No crossings by construction:** metabolites off the main cycle are duplicated only where needed. On your file that is acetyl-CoA in three places, the same result as my manual search earlier.
- **Side routes** are split inside or outside the circle so none interleave, and each side is drawn with a method that is mathematically guaranteed not to cross. Long-span routes go inside as straight chords; short local ones go outside and follow the circle. Small branches are attached in free gaps, and every later adjustment rejects any move that would create a crossing.
- **Fixed labels:** a final step only scales the drawing up until no box overlaps another and no line passes through a box. Scaling can't create crossings, and the script refuses to write output if any check fails.
- **Reproducible:** nothing is random, and Python's hash seed is fixed at startup.
- Networks without a cycle get a radial tree layout, and separate sub-networks are laid out side by side.

**The problem:** the clean-up gets the TCA drawing down to one remaining fault, a line through the acetoacetyl-CoA box in the crowded corner between succinyl-CoA and succinate. Four alternative enzymes connect those two metabolites, plus the acetoacetate branch. That one fault makes the final scaling enlarge the whole drawing about 7.7×, so the result is several metres across with tiny boxes.

**Next steps:**
1. Fix the last fault locally, by spacing the alternative enzymes for one step along the circle as well as inward and outward, instead of enlarging everything.
2. Shorten the curves outside the circle, which still swing too wide.
3. Test on a few other SBML files: a genome-scale BiGG model, one with several compartments, and one without a cycle.
4. Confirm byte-identical output across runs and deliver the script with sample drawings.

Even when finished, an automatic layout won't match the hand-tuned `tca_cycle.png` for your TCA model; the general script trades some polish for working on any file.

---

## 9. User

Continue

## 9. Claude

The general script now works on your TCA file and on a second, cycle-free model, but the larger E. coli core model was still running when I hit the tool limit. I haven't delivered it yet; your earlier files are unchanged.

**Your TCA file:** passes every check (0 crossings, 0 lines through boxes, 0 overlapping boxes, no bends). The drawing is a normal size again: placing small branches where they cause the fewest problems removed the last fault that had forced the 7.7× enlargement.
- **Main cycle:** found automatically as the nine TCA intermediates. Isocitrate lyase and malate synthase run across the middle, and parallel enzymes sit on both sides of the circle.
- **Duplicates:** acetyl-CoA is drawn four times in total.
- **Reaction names:** taken from the file's EC groups, choosing the group whose name mentions the reaction's own compounds. Isocitrate dehydrogenase (NADP+) is now named correctly instead of "homoisocitrate dehydrogenase". One name is still wrong: one aconitase step is labelled "2-methylcitrate dehydratase". I've changed the matching since, but haven't rerun to check it.

**COBRA test model (18 reactions, glycolysis):** passes every check. Glycolysis isn't a cycle, so it gets a radial tree layout, with exchange and transport reactions as separate pieces.
- An earlier version drew a false "cycle" that went through PEP and glucose 6-phosphate together. Stopping "phosphate" in compound names from counting as a match fixed it.
- Long names such as "Nicotinamide adenine dinucleotide - reduced" are now recognised as side compounds.

**E. coli core (95 reactions, two compartments):**
- It exposed three problems, all now fixed in the code:
  - The biomass pseudo-reaction connects to nearly everything and pulled the cycle through it. Reactions with seven or more compounds are now given their own copies of each compound.
  - A cycle could reuse the same reaction for two steps; now each step must use a different one.
  - Some branches could slip to the wrong side of the circle; the method now prevents that.
- The layout reaches the final clean-up in about 22 seconds, but the clean-up was too slow, so I reduced it for networks over 80 nodes. That run was still going when I ran out of tool calls. The final enlargement step still guarantees the rules either way, at the cost of a larger drawing for big models.

**Next steps:**
1. Finish the E. coli core run and check its size and readability.
2. Confirm the aconitase label fix.
3. Confirm two runs give identical files.
4. Remove the timing printouts and deliver the script with the three test drawings.

**Limitations to expect:**
- Duplicated metabolites increase quickly on big models: 57 copies in an early E. coli run, before the biomass fix.
- Genome-scale models (thousands of reactions) will be slow and hard to read as a single drawing.
- The automatic cycle choice is a heuristic; `--cycle` lets you set it yourself.

---

## 10. User

I don't find the script

## 10. Claude

The script is now in your outputs as `sbml_cycle_layout.py`, with its drawing of your model as `tca_auto.svg` and `tca_auto.png`.

**Usage:**
```
python3 sbml_cycle_layout.py model.xml [--out name] [--cycle "A,B,C,..."] [--side "x,y"] [--keep "x,y"]
```
- It needs only networkx (3.1 or later), numpy and matplotlib.
- It reads SBML levels 1–3, plain or `.gz`.
- `--cycle` sets the main cycle yourself.
- `--side` and `--keep` change which compounds are written inside reaction boxes.
- Progress messages go to stderr; the rule check result goes to stdout. If any rule fails, no output is written.

**Your TCA file:**
- Passes every check: 0 crossings, 0 lines through boxes, 0 overlapping boxes, no bends.
- Two runs gave byte-identical SVG and PNG files.
- It takes about 2½ minutes.
- The main cycle was found automatically, and both mislabelled reactions are now correct.

**COBRA glycolysis test model (18 reactions):** passes every check. It has no cycle, so it gets a tree layout.

**Not yet confirmed: E. coli core (95 reactions).** The layout stage finishes in about 20 seconds, but the clean-up after it hadn't finished when my test run was cut off. Expect long run times and large drawings on models of that size or bigger.

**Compared with the hand-tuned `tca_cycle.png`:** the automatic drawing is busier and some lines are longer, for example the acetate and citryl-CoA routes run across the inside of the circle. `tca_cycle_layout.py` remains the better choice for this particular model; `sbml_cycle_layout.py` is the one that works on any file.

---

## 11. User

Some models have cycles, some don't. I think cycles have to be drawn separately, then integrated back with the non-cycle parts

## 11. Claude

I agree, and it fixes the main weakness of the current script. Right now it picks one main cycle per connected network and opens every other cycle by duplicating a metabolite, so a model with two real cycles (say TCA and the urea cycle) gets only one circle.

**How I'd structure it.** The standard way to split a network into cyclic and non-cyclic parts is its *biconnected components*: maximal pieces that stay connected if any single node is removed.
- Any piece containing a cycle is a **cyclic module**.
- Everything else is single connections, which together form the **tree parts** (linear pathways and branches).
- Modules and tree parts meet only at single shared metabolites, the *cut points*. The map of which piece meets which at which cut point is always a tree (the *block-cut tree*).

The pipeline would then be:
1. **Expose the modules.** Duplicate hub metabolites first (pyruvate, acetyl-CoA, glutamate and the like, chosen by connection count). Otherwise genome-scale models collapse into one giant module.
2. **Draw each cyclic module on its own.** Its main cycle goes on a circle, with any extra routes inside it inside or outside the circle, using the current method.
3. **Draw the tree parts on their own.** Linear pathways become straight chains rather than the current radial fans, which suits glycolysis-style pathways much better.
4. **Assemble along the block-cut tree.** Each piece connects to the rest at only one metabolite, so it can be rotated into a free angular gap at that metabolite, outside its parent circle. That can't create crossings, so the no-crossing guarantee holds at the assembly level as well as within each piece.
5. **Clean up and check** as now: spread nodes, repair faults, scale only if necessary, then verify all rules.

**Further benefits:**
- Fewer duplicates.
- Every genuine cycle drawn as a circle.
- Faster: each module is small, so the slow clean-up runs on dozens of nodes instead of hundreds. That should also solve the E. coli core run time.

**The case to decide:** cycles that share two or more metabolites form a single module. In your model, the glyoxylate shunt shares isocitrate, succinate and malate with the TCA cycle. They can be drawn as:
- **One circle with the shunt across it**, as now. It shows the true connectivity with no extra copies.
- **Two separate circles**, with the shared metabolites drawn twice. Each cycle is cleaner on its own, at the cost of duplicates.

Which do you prefer for overlapping cycles? With your answer I'll restructure the script along these lines and test it on your TCA file, the glycolysis model and E. coli core.

---

## 12. User

Is there anything from https://github.com/forxhunter/MetaCarto that could improve your code?

## 12. Claude

Yes. MetaCarto tackles the same problem with a different philosophy, and several of its ideas address exactly the weak points of my script. I read its README and design document, not the source code itself.

**What it does:** it reads BiGG/SBML models and aims to draw them as a curator would, with straight pathway backbones, cycles as rings and cofactors as side branches. The layout is constructive and deterministic: no annealing and no random seed, so the same model always produces the same map. That matches your reproducibility rule.

**Ideas worth adopting:**

1. **Choosing the main substrate → product pair by chemistry, not names.** My script pairs compounds by name similarity, which needed patches (the "phosphate" problem in glycolysis). MetaCarto scores each substrate/product pair by shared heavy atoms from the molecular formulas, weighting carbon most, and breaks ties by lower overall connectivity.
   - Most SBML models carry formulas, and your BRENDA file has InChI strings that contain them, so this would generalise much better.
   - Two refinements matter. A cofactor is never picked as a main compound while a non-cofactor alternative exists; otherwise CoA → acetyl-CoA, sharing 21 carbons, would beat pyruvate → acetyl-CoA.
   - When chemistry is nearly tied, as with citrate synthase, candidates are re-scored by whether they close a cycle, join two components or extend an existing chain, which correctly picks oxaloacetate → citrate.

2. **A computed cofactor list instead of a fixed one.** A compound counts as a cofactor if it is among the most connected 10%, is on a curated list, or is never chosen as a main compound. The name list stays only as a starting point. This would replace my hard-coded side-compound list, which already missed several variants.

3. **Rings detected as separate motifs.** This is the idea you proposed last turn. MetaCarto looks for directed cycles of 3–12 steps with no shortcuts and at least two outside connections, contracts each into a single unit with known entry and exit points, and expands it back onto a circle at the end.
   - Each ring is rotated so its entry point faces the incoming pathway and its exit faces the outgoing one. My script places the cycle with no regard to where the rest of the network attaches.
   - It also settles the open question about overlapping cycles. It keeps only rings that share no nodes, preferring longer ones. In your model the glyoxylate shunt would therefore be a chord of the TCA ring, not a second circle.

4. **Layered placement for everything that isn't a ring.** My radial tree for linear pathways is the weakest part of my output. MetaCarto lays the contracted network out in layers, with flow running downward, then uses the Brandes–Köpf method to place nodes. The design document stresses that this step is what straightens long chains such as glycolysis. networkx doesn't include this, so it would mean writing roughly 300 lines.

5. **Direction from flux.** Models store reversible reactions in whatever direction the curator wrote them, so MetaCarto runs parsimonious FBA and draws each reaction in the direction it actually carries flux.
   - That would have fixed the four backwards TCA reactions I flipped by hand.
   - It needs cobrapy and a model with an objective, which BRENDA pathway files don't have. In that case the fallback in their design applies: choose directions that minimise backward edges.

6. **Clustering by annotation before structure.** Reactions are grouped into maps by the model's subsystem annotation first, then a KEGG lookup, and only failing both by network structure.
   - My script already reads SBML groups for reaction labels; the same data could split large models into readable per-pathway drawings.
   - This is the realistic answer to E. coli core and anything bigger. MetaCarto itself states that whole-model maps aren't print figures and that the per-cluster maps are the readable output.

7. **More quality checks.** It scores every map on orthogonality, crossings per edge, longest straight run, node separation, three kinds of label collision, local density, occupancy and aspect ratio. Adding the density, separation and aspect-ratio checks to my rule checker would catch cramped or badly shaped drawings, not just rule violations.

**What not to copy, because it conflicts with your rules:**
- **Crossings:** its crossing target is under 0.05 per edge, not zero. My planarity step should stay.
- **Labels:** where a dense cluster leaves no room, it drops a cofactor label rather than forcing it. That breaks your fixed-label rule; my labels-inside-boxes approach avoids it.
- **Bends:** every edge becomes an orthogonal line with up to two bends. Your rules prefer no bends, so straight lines and ring arcs should stay.
- **Hidden compounds:** H⁺ and H₂O are removed from every map. My script keeps them in the reaction box text.

**Licence:** the code is released under CC BY 4.0, which requires attribution, and the author notes it is a content licence with no patent grant and is not OSI-approved. If we borrow code rather than just ideas, your script must credit the repository. Reimplementing the ideas from the design document avoids most of that, but citing it is still good practice.

**Proposed next version:**
1. Formula-based main pairs with cofactor tiers, plus the computed cofactor list.
2. Ring motifs with no shared nodes, each drawn separately and rotated toward its connections, as you suggested.
3. Layered placement of the contracted network for the parts that aren't rings.
4. Assembly with the crossing, label and box checks that enforce your rules.
5. Splitting by subsystem or groups for large models.

Flux-based direction would be an optional extra when cobrapy is installed. Shall I go ahead in that order?

---

## 13. User

Yes, please go

## 13. Claude

The new script is `sbml_pathway_layout.py`, delivered with your TCA drawing (`tca_pathway_v2.svg` / `.png`) and two test outputs. It's much faster than the previous version: 8 seconds for your TCA file, about 30 for all of E. coli core.

| Model | Result |
|---|---|
| Your TCA file (24 reactions) | One 9-member TCA ring, side pathways placed radially next to the metabolites they attach to, no lines cut, 0 bends |
| Glycolysis test (18 reactions) | No cycle, so a layered drawing with straight vertical backbones |
| E. coli core (95 reactions) | 7 drawings; the pentose phosphate pathway comes out as a 7-member ring with glycolysis rising straight from it (`ecoli_core_pentose_phosphate.png`) |

All drawings pass the checks (0 crossings, 0 lines through boxes, 0 box overlaps), and two runs gave byte-identical files.

**What changed:**
1. **Main compounds chosen by chemistry.** Formulas are read from the file, either directly or from InChI strings; all 34 were found in yours and all 72 in E. coli.
   - The CoA part is ignored when comparing acyl-CoAs, so succinyl-CoA pairs with succinate.
   - Cofactors can never be the main compound, and near-ties go to the pair that continues the pathway.
   - This chose oxaloacetate → citrate for citrate synthase and isocitrate → succinate for isocitrate lyase without any hand-tuning.
2. **Cofactors** are found from the built-in list plus very highly connected compounds. Other secondary compounds become small dashed side boxes, such as glyoxylate and the acetyl-CoA copies.
3. **Rings** are found separately and never share metabolites; each is rotated toward the pathways it connects to.
4. **Layered placement** handles everything that isn't a ring: straight chains, with bends only where a line spans several layers (10 in the glycolysis test).
5. **Large models** are split by subsystem, or by network communities when the file has no subsystems (the case for this E. coli file).

The script credits MetaCarto in its header as the source of these ideas; no code was copied from it.

**Limitations:**
- Directions come from the file plus a step that minimises backward-pointing edges. Flux-based direction would need cobrapy, which isn't installed here.
- Ring steps with three or more alternative enzymes put the extras further out, which can make long lines (succinate–CoA ligase (GDP-forming) in your TCA drawing).
- The glycolysis test splits into several side-by-side columns rather than one tall spine.

**Older files:** your outputs folder still contains every earlier version. `tca_cycle_layout.py` (hand-tuned, TCA only) still gives the most polished TCA picture. `sbml_pathway_layout.py` replaces `sbml_cycle_layout.py` as the general script.

---

## 14. User

Could the layout be exported by the script to be reused later?

## 14. Claude

Yes. The script now exports every layout in two formats alongside the SVG/PNG drawings. For your TCA model those are `tca_pathway_v2.layout.json` and `tca_pathway_v2.layout.sbml`, delivered with the updated `sbml_pathway_layout.py`.

**`NAME.layout.json`: for reuse with this script.** It stores, for every drawing:
- each box's position, size, type, label and cofactor line, with the SBML id it refers to;
- every connection, including bend points and which lines are ring arcs;
- each ring's centre and radius, and the rule-check results.

To redraw from it without recomputing the layout:
```
python3 sbml_pathway_layout.py --from-layout NAME.layout.json [--out NEWNAME]
```
A redraw of the TCA layout is byte-identical to the original. You can also edit positions in the JSON by hand. The redraw re-checks the rules and warns you if your edits broke one, but still draws the result so you can see the problem.

**`NAME.layout.sbml`: for other SBML tools.** This is a copy of your model with the drawing embedded using the SBML Layout package, the standard way to store diagrams in SBML files.
- Each metabolite box becomes a species glyph, with a separate glyph for each duplicate copy.
- Each reaction box becomes a reaction glyph, with its connections marked as substrate or product (side compounds as side-substrate or side-product).
- Ring arcs are stored as cubic Bézier curves, and labels as text glyphs.
- Large models get one layout per drawing (7 for E. coli core).
- I validated the output with libSBML: 0 errors on your TCA file and on E. coli core. The glycolysis test file has one error, but it is already in the original file.

**Limitations:**
- SBML Layout export needs a Level 3 model; for Level 2 files only the JSON is written.
- The SBML file stores geometry only. Colours and line styles would need the separate SBML Render package, which I haven't added.
- Reusing a layout after the model changes isn't supported yet. That would mean keeping saved positions for unchanged parts and placing only the new reactions, which is possible but would need more work.

---

## 15. User

*(uploaded: test2.svg, brenda_urea-cycle.xml)*

I use `sbml_pathway_layout.py` but it splits reactions for some cases
Example in test2.svg from the SBML file attached

## 15. Claude

Fixed. The updated `sbml_pathway_layout.py` now draws your urea cycle model as one connected network (`urea_cycle.svg` / `.png`).

**What caused the split:** each reaction keeps one main substrate → product pair, and every other compound used to become a private dashed copy. That disconnected compounds that link different parts of the network:
- **Carbamoyl phosphate** is made by carbamoyl-phosphate synthase and carbamate kinase and used by ornithine carbamoyltransferase. Its synthesis ended up drawn on its own.
- **Urea** is released by arginase and used by urease and urea carboxylase. Its branch floated free.

**The fix:** a secondary compound now gets a line to its shared box whenever that line joins parts of the network that would otherwise be separate. Grey lines mark these connections; dark lines are main steps. If the parts are already connected, the compound stays a private copy. That keeps lines from crossing the ring, which matters in your TCA model, where acetyl-CoA, acetate and glyoxylate would otherwise cut across it.

**Two bugs found while testing, also fixed:**
- Reactions sitting on a ring could be moved off it during clean-up, leaving them visually disconnected.
- A shared box whose every connection had been replaced by a copy was left behind as an empty duplicate.

**In your urea drawing:**
- Carbamoyl phosphate and urea are single shared boxes, with all their reactions attached.
- Citrullinase now sits inside the ring, beside the ornithine–citrulline step, away from the carbamoyl phosphate connection.
- Only L-aspartate and fumarate remain as dashed copies; each appears in just one reaction.
- No connections had to be cut to avoid crossings.

**Regression tests:** TCA, E. coli core and glycolysis all still pass every check, the urea drawing is byte-identical across runs, and the SBML Layout exports validate with 0 errors. The TCA drawing changed slightly: the crossing repair now makes 2 copies instead of 0, and the citrate-CoA ligase branch has longer lines than in the previous version.

---

## 16. User

I would also like to see some compounds as Urea is in the graph. E.g. NH3 is important for the urea cycle. So an option would be nice to see or not some specific compounds

## 16. Claude

The updated `sbml_pathway_layout.py` has two new options for this, and `urea_cycle_show_NH3.svg` / `.png` shows the urea cycle drawn with NH3 shown.

```
python3 sbml_pathway_layout.py brenda_urea-cycle.xml --show "NH3"
python3 sbml_pathway_layout.py model.xml --show "NH3,CO2" --hide "L-aspartate"
```
- **`--show`** draws the listed compounds as shared boxes connected to their reactions, like urea. They're taken out of the reaction-box text and are never used as the main step of a reaction, so they can't distort the ring.
- **`--hide`** always writes the listed compounds inside the reaction box. It replaces `--side`, which still works.
- **Names** are matched loosely: case, spaces and hyphens are ignored. A name not found in the model gets a warning, so a typo isn't silently ignored.
- **Saved layouts:** the chosen options are recorded in the `.layout.json` file.

**With `--show NH3`:**
- NH3 is one shared box connected to urease, allophanate hydrolase and carbamate kinase.
- Carbamoyl phosphate, urea and allophanate stay connected exactly as without the option.
- NH3's other two uses (carbamoyl-phosphate synthase and citrullinase) are dashed copies. A line to the shared box would have crossed another line, and your no-crossing rule takes priority. The more compounds you show, the more such copies appear: 3 with `--show "NH3,CO2"`, for example.

**Fixes made along the way:**
- Shown compounds were first linked before deciding which other compounds connect separate parts. NH3 then made carbamoyl phosphate and urea look redundant, and they were split into copies. Shared compounds are now decided first, as if `--show` weren't given.
- Each shared compound's box is now placed by testing a grid of candidate positions and keeping the one where its lines cross nothing and pass through no box. This keeps drawings a normal size; one test had grown the drawing to many times its normal size.
- If several shown compounds would land on the same spot, they're spread apart.
- A line running almost exactly along a ring arc could become a measured crossing after rescaling, through floating-point rounding. The final repairs now repeat until all three checks pass together.

**Regression tests:** all test models still pass every check, the output is byte-identical across runs, and the SBML Layout files validate with 0 errors. The default urea drawing, without `--show`, is unchanged.

---

## 17. User

The cycles look good, but I would need the non-cycle reactions to be displayed more vertically without too many angles

## 17. Claude

The pathways branching off the cycles now run vertically, with straight vertical or horizontal connections to the ring. The updated script is `sbml_pathway_layout.py`, with new drawings `urea_cycle.png` and `tca_pathway_v2.png` (SVGs alongside).

**What changed:**
- **Top-to-bottom side pathways.** Each pathway attached to a ring gets its own top-to-bottom layout instead of being fanned out radially.
- **Placement by anchor position.**
  - Anchor at the ring's top or bottom: the pathway goes above or below, lined up so the connecting line is vertical.
  - Anchor at the left or right: the pathway goes beside the ring, level with the anchor, so the line is horizontal.
- **Ring rotation** now favours putting anchors at the top and bottom.
- **Shared compound boxes** (urea, NH3 with `--show`) are placed preferring horizontal or vertical lines.
- **Crowded rings** are retried with a larger radius instead of enlarging the whole drawing, which had made the TCA drawing balloon on the first try.

**Urea cycle:** the urea branch (urease, urea carboxylase, allophanate) rises straight from arginase, and carbamoyl phosphate joins ornithine carbamoyltransferase with a horizontal line. Aspartate and fumarate also connect straight, and no lines had to be cut.

**TCA:** the citryl-CoA bypass sits vertically above citrate, and malate synthase and maleate isomerase connect horizontally beside the ring.

The step from succinyl-CoA to succinate still has two long diagonal lines. Four enzymes catalyse that one step, and the extra ones are placed inside and outside the ring, the furthest well below it. Those are alternative enzymes within the cycle, not side pathways, so this change didn't affect them. If you'd like, I can stack them in a short vertical column under the step instead.

**Regression tests:** all test models still pass every check (E. coli core, the glycolysis model, and the `--show NH3` run included). Output is byte-identical across runs, and the SBML Layout files validate with 0 errors.

---

## 18. User

I need an option to see the legend of the drawings at the bottom right of the graphs.

## 18. Claude

The updated `sbml_pathway_layout.py` has a `--legend` option that puts a legend at the bottom right of each drawing. The examples are `urea_cycle_legend.png` and `tca_pathway_legend.png`, with SVGs alongside.

```
python3 sbml_pathway_layout.py model.xml --legend
python3 sbml_pathway_layout.py --from-layout NAME.layout.json --legend    # also works on a saved layout
```

**What it shows:**
- **Only the elements used in that drawing,** so a drawing without a cycle won't list "cycle intermediate" or "cycle step".
- **Boxes:** cycle intermediate, metabolite, and dashed copy, meaning a compound drawn elsewhere too or used only once.
- **Reaction box:** enzyme name, with the cofactors consumed → released on its second line.
- **Lines:** cycle step (red), main substrate → product (dark), secondary compound (grey).

Each sample is drawn with exactly the same style as the drawing, so the two can't drift apart if styles change later.

**Placement:** the legend goes inside the drawing's bottom-right corner when nothing is there, as in the TCA drawing. Otherwise the canvas is extended downward and the legend placed at the bottom right, so it never covers a box or line. That happened in the urea drawing, where the L-aspartate connection occupied the corner.

**Checks:**
- Without `--legend`, the output is byte-identical to before.
- With it, two runs and a redraw from the saved layout gave identical files, and all 7 E. coli core drawings work.
- The legend is added after the rule checks, so it can't affect the layout.
- The setting is recorded in the `.layout.json` file. It's left out of the `.layout.sbml` file, since a legend isn't part of the SBML Layout standard.

---

## 19. User

Is the layout.json output compliant with Escher https://github.com/opencobra/escher ?

## 19. Claude

No. `NAME.layout.json` is the script's own format, not an Escher map. I checked it against Escher's official schema (`jsonschema/1-0-0` in the repository) with a JSON Schema validator, and it fails right away, because the overall structure is different.

**How they differ:**

| | `layout.json` (this script) | Escher map |
|---|---|---|
| Top level | one object with `drawings`, `settings`, `model` | an array `[header, body]`, with fixed `schema` and `homepage` values in the header |
| Compounds | boxes with position, size, type and label | `nodes` with `node_type: "metabolite"`, `bigg_id`, `name`, label position, and `node_is_primary` |
| Reactions | boxes connected by links | a list of metabolites with coefficients plus `segments` between nodes; the reaction itself is a small `midmarker` point, not a box |
| Curves | straight lines, ring arcs, bend points | each segment is one cubic Bézier (`b1`/`b2`); bends need extra `multimarker` nodes |
| Y axis | up | down |
| Several drawings | all in one file | one map per file |

**Can it be made compatible?** Yes, as an additional export like the SBML Layout one; I'd add a `--escher` option that writes one Escher map per drawing:
- Metabolite boxes become metabolite nodes, and dashed copies become extra nodes with `node_is_primary: false`.
- Reaction boxes become midmarkers, with the label placed at the box position.
- Ring arcs become Bézier segments, and bend points become multimarkers.
- The y axis is flipped.

The export would be checked against the schema you linked, and I can also check that Escher's own Python package (`pip install escher`) loads the maps.

**Trade-offs to be aware of:**
- **Escher draws in its own style:** reactions as lines with a small dot and a floating label, not as boxes. The positions carry over, but box shapes, colours, the cofactor line inside boxes and the legend don't.
- **Your fixed-label and no-overlap rules** are checked against the box sizes, which Escher doesn't use. In Escher, labels could therefore sit slightly differently relative to lines than in the SVG.
- **Data overlays** (fluxes, expression data) work when node `bigg_id`s match the model's ids. I'd use the SBML species and reaction ids with their `M_`/`R_` prefixes removed, which matches BiGG-style models such as E. coli core. For BRENDA files the ids are BRENDA's own, so overlays would need data keyed the same way.

Would you like me to add the `--escher` export?

---

## 20. User

Yes, add the `--escher` export

## 20. Claude

The updated `sbml_pathway_layout.py` now has an `--escher` option that writes each drawing as an Escher map. Examples: `urea_cycle.escher.json` and `tca_pathway.escher.json`.

```
python3 sbml_pathway_layout.py model.xml --escher
python3 sbml_pathway_layout.py model.xml --from-layout NAME.layout.json --escher   # from a saved layout
```
- It writes one `NAME.escher.json` per drawing, so E. coli core gives 7 maps.
- To open one, go to escher.github.io, then Map → Load map JSON.
- Redrawing from a saved layout also needs the SBML file, because Escher maps carry each reaction's full stoichiometry.

**How the drawing maps onto Escher:**
- **Metabolite boxes** become metabolite nodes. Dashed copies become extra nodes marked as secondary (`node_is_primary: false`), which is Escher's own notion of a duplicated compound.
- **Reactions** become midmarkers, with all participants and coefficients listed, cofactors included.
- **Ring arcs** become Bézier segments, at most 90° each, and bend points become multimarkers.
- **Ids** have the SBML `M_`/`R_` prefixes removed, so E. coli core maps use BiGG ids and Escher can overlay flux or expression data on them.

**Validation:**
- All 9 maps (urea cycle, TCA, 7 × E. coli core) pass the official schema from the repository you linked.
- They also pass Escher's own `validate_map`. This additionally checks that every segment connects existing nodes and that every connected metabolite has non-zero stoichiometry in its reaction.
- I plotted a map independently from its nodes and Bézier segments, and it matches the drawing.
- Runs are reproducible: only `map_id` changes, because it is taken from the output name.
- Without `--escher`, the other outputs are byte-identical to before.

**Not tested:**
- **Escher's viewer itself:** I couldn't run it here, so the first real load should be in Escher.
- **Fonts and labels:** Escher uses its own styles and fonts, which are larger than mine, so labels can run into each other even though the positions are taken from the drawing. Reaction labels are placed inside the space their box occupied, which is free of other boxes and lines, but text bigger than the original box can still overflow it.
- **What doesn't carry over:** box shapes, colours, the cofactor line and the legend; Escher shows the cofactors through the reaction's stoichiometry instead.

---

## 21. User

How the reaction direction is inferred?

## 21. Claude

The script takes each reaction's direction from the SBML file exactly as written, and nothing currently changes it. Three steps involve direction, but only the first one sets the arrows:

**1. Arrows follow the written direction.**
- For each reaction, the main pair is always a reactant → a product (from `listOfReactants` to `listOfProducts`). The arrow runs metabolite → reaction box → metabolite in that order.
- Secondary compounds and dashed copies point into the reaction if they are reactants and out of it if they are products.
- The `reversible` attribute is read, but it doesn't change any arrow. A reversible reaction is drawn in its written direction with a single arrowhead.

**2. The ring is oriented by majority vote.**
- Each ring is ordered so that most of its steps, as written, run with the ring.
- A step written the other way still keeps its own arrow, so it appears to point against the flow of the ring.
- Mirroring the ring during rotation only changes clockwise vs. counter-clockwise on the page, never an arrow.

**3. The layered part only borrows direction for placement.**
- To lay out top to bottom, the network must have no loops. A standard heuristic (Eades–Lin–Smyth) temporarily reverses the fewest connections needed.
- That only decides what goes above what. When the lines are drawn, each reversed connection is flipped back, so arrows still show the written direction.
- In an earlier summary I described this as a step "that minimises backward-pointing edges". That was misleading: it minimises connections running *upwards* on the page, not arrows that point the wrong way.

**In the exports:**
- **Escher:** the `reversibility` field is copied from the file. Escher then uses that and the coefficient signs to draw its own arrowheads, so reversible reactions may show arrowheads at both ends there, unlike in my SVG/PNG.
- **SBML Layout:** substrate and product roles follow the written direction.

**Limitations:**
- **For BRENDA files, written direction is often arbitrary.** In your TCA file every reaction is marked reversible, and some are written against the textbook direction. In the old TCA-only scripts I fixed four of them with a hard-coded list; the general script does no such correction.
- **No flux-based direction.** I mentioned using pFBA flux directions as an option, but it isn't implemented: it needs cobrapy and a model with an objective, which BRENDA pathway files don't have.

**Options I could add:**
- **`--orient-reversible`:** flip reversible reactions to agree with their ring's direction, or with the overall top-to-bottom flow outside rings. This is a general version of the hand-made TCA fixes; irreversible reactions would never be flipped.
- **`--show-reversible`:** draw reversible reactions with arrowheads at both ends, matching Escher's convention.
- **`--flux-direction`:** when cobrapy is installed and the model has an objective, run pFBA and draw each reaction in the direction it carries flux. Reactions with zero flux would keep their written direction.

Would any of these help?

---

## 22. User

Could you save this chat entirely in a file?

## 22. Claude

*(This file.)*
