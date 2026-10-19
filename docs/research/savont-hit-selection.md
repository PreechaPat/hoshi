# How Savont picks one reference per ASV

Investigated 2026-10-19 for savont **v0.8.0** (the version installed in the `bio`
env), in the context of
[savont#18](https://github.com/bluenote-1577/savont/issues/18). Source links are
pinned to the `v0.8.0` tag commit
[`4e242ed`](https://github.com/bluenote-1577/savont/tree/4e242ed9e20f65583759e2eedd8a18dd58d02295);
`src/classify.rs`, `src/merge.rs` and `src/taxonomy.rs` are unchanged on `main`
(`5a5a1ce`) at the time of writing.

Short answer: selection happens in two stages.

1. **Alignment filter:** only hits with the same mismatch count (`nm`) as
   minimap2's first hit are kept. Every row of an ASV in `asv_mappings.tsv`
   therefore has the same `number_mismatches`; they are all alignment ties by
   construction.
2. **EM decides among them:** `savont classify` picks the hit whose `tax_id` has
   the highest **EM-estimated abundance across the whole sample**.

In this note **"tie" means an EM-abundance tie**: two or more of those
candidates whose `tax_id`s have exactly equal EM abundance. On an EM tie,
classify picks the **last** tied hit. `asv_mappings.tsv` lists hits highest-EM
first, with EM ties kept in minimap2 order, so:

- no EM tie → the first row **is** savont's pick;
- EM tie → the first row is **not** savont's pick (the last tied row is).

`savont export` does not run EM. It reads `asv_mappings.tsv` and keeps the
first row per ASV, so it agrees with classify except on EM ties.

## Pipeline

All links below point to savont `v0.8.0` source at commit `4e242ed`.

### 1. Candidate hits: equal mismatch count to minimap2's first hit

[`collect_best_mappings`, classify.rs L138–L229](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L138-L229)

- Each ASV is aligned with minimap2 (`map_ont`).
- `min_nm` is taken from `mappings.first()`
  ([L186–L190](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L186-L190)), i.e. the first record minimap2
  returns, **not** a true minimum over all records despite the comment
  "Find the minimum NM value".
- Every alignment with `nm == min_nm` is kept
  ([L195](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L195)); others are dropped. These kept hits are exactly
  the rows that later appear in `asv_mappings.tsv`.
- Identity is `100 * (1 - nm / aligned_query_length)`
  ([L197–L198](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L197-L198)). Hits share `nm` but can differ in
  alignment length, so identity can differ between candidates; identity plays no
  part in choosing between them.
- The lookup key for EMU databases is the `tax_id` before the first `:` of the
  reference name ([taxonomy.rs L724–L727](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/taxonomy.rs#L724-L727), wired in
  [databases.rs L26](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/databases.rs#L26)). Several references of one species
  therefore share one key, which is why the same `tax_id` can appear on several
  rows of one ASV.

### 2. EM over all ASVs, per `tax_id`

[`run_em_algorithm`, classify.rs L24–L134](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L24-L134), called at
[L539–L545](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L539-L545)

- One abundance value per unique key (`tax_id`)
  ([index built at L506–L516](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L506-L516)), initialised uniformly
  ([L31](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L31)).
- Each iteration, every ASV splits its read depth over its candidate hits in
  proportion to the current abundance of each hit's `tax_id`
  ([L61–L73](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L61-L73)); totals are divided by total reads
  ([L79](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L79)). Duplicate rows of the same `tax_id` within an ASV
  each count separately in this split.
- Stops when the largest change is below `0.1 / total_reads` or after 1000
  iterations ([L116](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L116), threshold at
  [L539](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L539)).
- Abundances below that threshold are then set to exactly `0.0`
  ([L126–L132](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L126-L132)).

This is the "prefer taxonomy seen across all ASVs" behaviour the author
describes in the issue: a `tax_id` supported unambiguously by other ASVs gains
abundance and wins the shared ASVs.

Exact ties come from symmetry. If two `tax_id`s are always hit by the same set
of ASVs (e.g. strains with identical 16S, or one species under two names), they
receive identical updates every iteration and stay bit-for-bit equal. Several
taxa that were all zeroed in the last step also tie at `0.0`.

### 3. Writing `asv_mappings.tsv` (all candidates)

[classify.rs L574–L612](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L574-L612), written at
[L710–L711](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L710-L711) by
[`write_asv_mappings`, taxonomy.rs L903](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/taxonomy.rs#L903)

- The ASV's hits are sorted by EM abundance, descending
  ([`sort_by`, L576–L580](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L576-L580)). `slice::sort_by` is stable
  ([Rust docs](https://doc.rust-lang.org/std/primitive.slice.html#method.sort_by)),
  so tied hits keep minimap2's order.
- Every sorted hit becomes one row
  ([L603–L611](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L603-L611)). The EM abundance itself is **not**
  written; it only appears in the debug log as
  `Hit #n: ... EM_abundance=...` ([L583–L591](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L583-L591)).

### 4. The final pick (`species_abundance.tsv`, `genus_abundance.tsv`)

[classify.rs L615–L648](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L615-L648)

```rust
// Take the mapping with highest EM-estimated abundance (or first if tied)
let best_mapping = asv_mappings
    .iter()
    .max_by(|a, b| tax_abundances[a.index].partial_cmp(&tax_abundances[b.index]).unwrap())
    .unwrap();
```

- `asv_mappings` here is the **unsorted** candidate list (minimap2 order).
- `Iterator::max_by` returns the **last** element when several are equally
  maximum ([Rust docs](https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.max_by)),
  so the code comment "or first if tied" is wrong.
- Taxonomy levels are then filled or set to `UNCLASSIFIED` from the chosen hit's
  identity: species at `>= 99.0`, genus at `>= 94.5` (defaults,
  [cli.rs L205–L210](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/cli.rs#L205-L210); applied in
  [`from_taxonomy_entry`, taxonomy.rs L592–L640](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/taxonomy.rs#L592-L640)).
- `species_abundance.tsv` sums ASV depth per lineage string, not per `tax_id`
  ([key built at taxonomy.rs L791](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/taxonomy.rs#L791)), so the chosen `tax_id`
  is not recoverable from it.

Consequence: the stable sort keeps EM-tied hits in the same relative order in
the file as in the candidate list, so the **last EM-tied row** in
`asv_mappings.tsv` is the hit classify picks. With a single top EM value, the
first row is the pick.

### 5. `savont export` takes the first row

[`export`, merge.rs L475–L489](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/merge.rs#L475-L489)

Export does no alignment or EM; it only re-reads classify's output files.
`hash_to_lineage.entry(hash).or_insert_with(|| lineage)` keeps the first
`asv_mappings.tsv` row per ASV and ignores the rest. Because that file is
already sorted by EM abundance, the first row is a top-EM hit, and it matches
classify's pick unless there is an EM tie.

## Evidence

Each sample was rerun with `savont classify -l debug -d ~/works/db/savont/emu-2026`;
in every case the rerun's `species_abundance.tsv` and `asv_mappings.tsv` were
identical to the original, so the debug EM values describe the original runs.
`savont export --no-fuzzy` was run on the original folder.

"Tied top EM" counts ASVs where two or more rows share the maximum EM
abundance; most of these are the same `tax_id` on several references, which
cannot change the result.

| Sample | ASVs with tied top EM | ...of which tied hits differ in `tax_id` | Classify chose last tied row | `species_abundance.tsv` rebuilt from first row / last tied row | Export == first row |
|---|---|---|---|---|---|
| `test_data/.../savont-out-soil` | not counted | 2 | 2/2 | ✗ / ✓ | ✓ |
| si_260109 HN141 | 59 | 5 | 59/59 | ✗ / ✓ | ✓ |
| si_260109 HG113 | 15 | 2 | 15/15 | ✗ / ✓ | ✓ |
| si_260109 HN122 | 38 | 4 | 38/38 | ✓ / ✓ ¹ | ✓ |
| si_260109 HN138 | 22 | 4 | 22/22 | ✗ / ✓ | ✓ |
| si_260109 HN146 | 33 | 5 | 33/33 | ✗ / ✓ | ✓ |

¹ HN122's differing ties were all species-`UNCLASSIFIED` with the same
lineage, so first vs last makes no difference to the output.

Examples where the first row is wrong (first row → classify's pick):

- soil `final_consensus_8`: Bradyrhizobium sp. 200 → Bradyrhizobium sp. 170
  (both EM 0.006358)
- soil `final_consensus_71`: Bacillus sp. FSL R5-0397 → Bacillus sp. FSL K6-0046
  (five hits at EM 0.000874)
- HN138 `final_asv_24`: Simiaoa sunii → Waltera intestinalis
- HG113 `final_asv_21`: Bacteroides sp. CACC 737 → Bacteroides humanifaecis

## Implications for hoshi

- The "last tied row" rule needs the EM abundances to know which rows are tied.
  They are not in `asv_mappings.tsv` and not in info-level logs (those only show
  `Starting EM algorithm with N taxa...` / `EM converged after N iterations`).
- Ways to get savont's actual pick:
  1. Run classify with `-l debug` and read
     `Classified ASV ASV_n: tax_id_acc=<tax_id>` ([L650](https://github.com/bluenote-1577/savont/blob/4e242ed9e20f65583759e2eedd8a18dd58d02295/src/classify.rs#L650)).
  2. Re-run the EM from `asv_mappings.tsv` + `feature-table.tsv` (they contain
     every input the EM uses) and apply max-with-last-on-tie. Exact ties depend
     on reproducing savont's float arithmetic.
- On a tie the pick is an artefact of minimap2 ordering, not evidence for one
  species, so reporting the tied set may be more honest than either rule.
