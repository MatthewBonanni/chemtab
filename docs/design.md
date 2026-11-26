# chemtab and Mallard's table format: joint design

Status: agreed 2026-11-25. Every decision in section 4 took the recommended option. Not yet implemented; decision 9's spike comes first.

Sources read: Mallard `origin/main` (`docs/design/table_format.md`, `docs/design/flamelet.md`,
`docs/guide/flamelet.md`, `src/chemistry/flamelet_table.{h,cpp}`, `src/physics/flamelet.{h,cpp}`,
`tools/docs/flamelet.py`) and FPVgen at `cbc9b6b` (2025-06-18), cloned read-only. Citations are
`file:line`; `FPVgen:` means `src/fpvgen/flamelet_table_generator.py` unless another file is named.

Two things are **decided** by the user and are recorded here as such (section 2.0): the format is
revised in place and stays `version = 1`, since no Mallard release has shipped it; and table names
follow Mallard's all-caps field convention. "Current format" below means what `origin/main` reads
today; "revised format" means this proposal.

The one-paragraph conclusion: the current format's *contract* (tags, units, species block, provenance, checked
physics) is right and stays. Its *layout* couples three things that should be separate: the file
(what a generator and a human want), the lookup record (what a GPU wants) and the set of fields the
solver happens to read today. Mallard's reader already copies and sometimes converts the table before
any lookup, so the file does not need to be the lookup layout. Proposed: one HDF5 dataset per
field, each naming its own axes; the reader packs the handful of fields the solver reads into a
device record of its own design. That makes the generator general (fields on a subset of axes,
mixed precision, any number of species at no solver cost) and makes the GPU layout a Mallard-only
decision that can change without touching the format.

---

## 1. The current format against what a general generator needs

### What it gets right (keep)

| Aspect | Where | Why it holds up |
|---|---|---|
| Model tag + thermodynamics tag, refused by name | table_format.md §3, §5.2, §6.2; `flamelet_table.cpp:63-87` | Orthogonal: FGM reuses `saghafian2015`; a `low-mach` tag reuses `fpv`. |
| Exact units strings, no conversion | §3.1; `flamelet_table.cpp:231-233` | A generator cannot be silently off by a density factor. FPVgen's `SRC_PROG` is 1/s and `HeatRelease` W/kg (`FPVgen:1468`, `:1478`); the format's are per volume. The check catches exactly that port error. |
| Tensor grid with a normalized progress axis | §5.1; flamelet.md §4.3 | `Lambda` is what turns the non-rectangular `(Z, C)` domain into a tensor grid with no search on the device. This is the right answer to "tensor vs non-tensor" and it generalizes (section 2.3). |
| Required provenance with the generator input verbatim and the mechanism SHA-256 | §3.2 | A table regenerates from itself. FPVgen has nothing like it: an auto-built file name (`FPVgen:1311-1317`) is its only record. |
| Readers ignore what they do not know | §3, §8 | Lets chemtab add diagnostics without a format change. |
| float64 axes, lookups summed in double, fixed corner order | `flamelet_table.h:86-127` | The bitwise contract; unaffected by anything below. |
| `omega_C` tabulated at `rho0`, scaled by the solver | §5.1; `flamelet.h:135-140` | Keeps the compressible scaling a solver decision. |

### What it gets wrong or leaves out

1. **The file layout is the lookup layout, and the species block rides in the hot record.** One
   `/values` array, node-major, every field in every node (§2). The solver reads 11 fields
   (`flamelet.h:89-91` plus `C`); a table with all 53 GRI species has a 64+-field node, ~520 B in
   float64, and the 11 hot values sit wherever the writer put them. The "2^d short runs" argument
   (flamelet.md §4.2) only holds for small field counts. FPVgen's default is *all* species and all
   production rates (`sample/input.toml:50-51`), so the first real tables will be the bad case. The
   design note itself lists node-major vs field-major as "to be measured" (flamelet.md §4.2).
2. **The layout is already not what runs.** float32 tables are widened to float64 while being copied
   to the device (`flamelet.cpp:25-57`), so float32 today halves the file and nothing else, contrary
   to table_format.md §2 ("halve the memory and bandwidth of lookups"). On GPUs every table is
   deep-copied (`flamelet_table.h:251-254`); `mmap` only stages that copy. The reader is one small
   step from packing its own record.
3. **`C` is stored as a 3-D field that must be linear in `Lambda`,** then checked to 1e-9
   (`flamelet_table.cpp:245-262`), and the map reads its two ends from 8 nodes that are far apart in
   memory (`flamelet_table.h:158-159`: offsets `0` and `(n_Lambda-1)*stride`, across 4 `(Z, Sz)`
   corners). The information is two 2-D fields, the least and largest `C` over `(Z, Sz)`. The format
   cannot say that, because every field must live on every axis. The same limitation blocks any other
   reduced-dimension data (per-flamelet diagnostics, stream states along a second mixture fraction).
4. **The map from transported variables to axes is hard-wired per model tag** (`FPV::locate`,
   `flamelet_table.h:154-164`). Heat loss (`h` normalized between adiabatic and coldest), two
   streams (`Z2/(1-Z1)`) and FGM are the same two operations, normalized variance and
   normalize-between-bounds, in a different order. The format has no vocabulary for it, so each becomes a
   new tag with bespoke reader code.
5. **Axis lookup is bisection** (`flamelet_table.h:64-76`), about 8 dependent loads per axis for 200
   nodes, and the point is located three times per cell per RHS (`mixture.h:68`,
   `mixture_viscous_flux.h:98`, `solver_flamelet.cpp:41`). Not measured, likely small next to the
   corner reads, but avoidable. FPVgen's coordinates are analytic (linear, power law, two-linear,
   linear-then-stretched: `coordinate.py:7-14`) and carry their type and parameters to the file
   (`coordinate.py:73-85`, `:154-163`, `:262-272`) because CharLES^X inverts them in O(1). The format
   drops that information.
6. **Saghafian's model is only half there.** The format tabulates `gamma0`, `a_gamma` but not the
   temperature exponents of viscosity and `lambda/cp` nor the activation temperature of the source,
   all of which FPVgen computes (`AMU`, `ALOC`, `TA`: `FPVgen:1452-1465`, `:1524-1529`). Mallard's
   `kappa/cp` therefore "does not depend on `T`" (flamelet.md §5.5, Boundaries) and `omega_C`
   scales with density only. Away from the flamelet temperature (shocks, compression) that is an
   accuracy loss the format currently makes unfixable.
7. **The transport model of the flamelets is optional provenance** (§3.2, "anything else"), but the
   solver assumes unity Lewis number (flamelet.md decision 8). FPVgen's sample default is
   `mixture-averaged` (`sample/input.toml:10`). A mismatched table loads silently.
8. **One pressure, as a root scalar** (§3). A pressure axis is deferred to "stage 6" with no format
   provision; it would need the root attribute redefined.
9. **Field-name namespaces are thin.** Only `Y_` is reserved. Species sources (`omega_<k>`, needed
   for transported NO, flamelet.md §1) and the production/consumption split have no convention.
10. **Alignment is "should, 8 bytes"** (§2). Irrelevant once the reader packs; insufficient
    (cache-line, page) if it does not.
11. **No structural link from fields to axes for tools.** `f["values"][..., fields.index("T0")]`
    (§9) works but no generic tool (xarray, ParaView, `h5dump -d`) can label it.

Things the current format is *not* wrong about, despite appearances: tensor grids (non-tensor tables buy nothing the
normalized axis does not); one table per file (multi-table files complicate provenance and have no
user); HDF5 itself.

---

## 2. The format, revised before its first release

### 2.0 Decided

- **`version` stays 1.** No release has shipped the table format, so it is redesigned in place:
  one format, one reader path, no migration and no compatibility with what `origin/main` reads
  today. The six fixtures in `test/data/flamelet/` are rewritten.
- **Field names are all-caps, as Mallard's output fields are** (`RHO`, `P`, `T`, `HRR`, `CP`,
  `Y_<species>`; #1468). Names are uppercased faithfully and keep their meaning: the table's
  reference-state quantities keep their `0`, so they cannot be taken for the solver's own `RHO` and
  `T`. Field names match `[A-Z][A-Z0-9_]*`; the species part of `Y_<species>` and
  `OMEGA_<species>` keeps the mechanism's spelling, as in Mallard's output (`"Y_" + species_names[k]`).

| Current name | Revised name | Note |
|---|---|---|
| `T0` | `T0` | Unchanged. Mallard's *output* catalog also has a `T0`, the stagnation temperature (`derived_catalog.h`), which is why the table's is written out as `T_TABLE`; table fields and output fields are separate namespaces, so nothing collides, but the docs must keep saying which is which |
| `e0` | `E0` | |
| `R` | `R` | |
| `gamma0` | `GAMMA0` | |
| `a_gamma` | `A_GAMMA` | |
| `rho0` | `RHO0` | Not `RHO` |
| `mu` | `MU` | |
| `kappa` | `KAPPA` | |
| `cp` | `CP` | Same quantity as Mallard's output `CP` |
| `omega_C` | `OMEGA_C` | |
| `HRR` | `HRR` | Unchanged; already Mallard's output name |
| `C` | `C_MIN`, `C_MAX` | Replaced by its bounds (section 2.1, decision 2) |
| `Y_<species>` | `Y_<species>` | Unchanged |
| new | `A_MU`, `A_KAPPA_CP`, `T_A` | Optional (decision 6) |
| new | `OMEGA_<species>`, suffixes `_PROD`, `_CONS` | Reserved namespace |

  In Mallard's reader this changes the contract tables (`flamelet_table.cpp:63-93`), the name
  check `is_name` (`:147`), the lookups by name in `flamelet.cpp` (`:69-95`: `"C"`, the thermo,
  transport and source names) and the tests and fixtures. No input key and no output name changes:
  `T_TABLE`, `HRR`, `Y_*`, `C_NORM` are already what Mallard writes.
- **Proposed with it (decision 3): axis names follow suit**, `Z`, `SZ`, `LAMBDA`. Structural names
  (groups and attributes: `model`, `species`, `pressure`, `progress_variable`) stay lower-case, as
  Mallard's TOML keys are; the one mixed-case attribute, `Z_st`, becomes
  `stoichiometric_mixture_fraction`.

### 2.1 Layout

```
/                          attrs: format = "mallard-table", version = 1,
                                  model, thermodynamics, species[], axes[],
                                  pressure (unless a "P" axis exists)
/axes/<NAME>               float64, 1-D, strictly increasing, >= 2 nodes; an HDF5 dimension scale
                           attrs: units
                                  quantity   what the solver enters this axis with ("Z", "ZVAR", "C", "H", "P")
                                  transform  "identity" | "normalized-variance" | "bounds" | "log"
                                  of         for normalized-variance: the axis whose variance this is ("Z")
                                  grid       optional, descriptive: kind + parameters (section 2.4)
/fields/<NAME>             float64 or float32; shape = node counts of its axes, C order
                           attrs: axes[]     names of its axes, a subset of the table's, in table order
                                  units
                           dimension scales attached; contiguous, unfiltered, little-endian
/model_data                as today, plus required attr transport = "unity-lewis" | "mixture-averaged" | "multicomponent"
                           and group pdf with one string attr per axis ("beta", "delta", "none")
/model_data/flamelets      optional diagnostics, one row per flamelet: chi_st, T_max, C_st, branch
/provenance                as today, plus required dependency versions (cantera, chemtab source revision)
```

- **The table's axes** are the datasets of `/axes`, ordered by the root attribute `axes[]`. A field
  lives on any ordered subset. Full-rank fields are the tabulated state; reduced fields carry
  bounds and other per-column data.
- **Bounds replace the linearity contract.** For an axis with `transform = "bounds"` and
  `quantity = Q`, the fields `Q_MIN` and `Q_MAX` exist on the axes that precede it, with
  `Q_MAX >= Q_MIN`. For `fpv`: `C_MIN(Z, SZ)`, `C_MAX(Z, SZ)`; the 3-D `C` field is gone, and
  `C = C_MIN + LAMBDA (C_MAX - C_MIN)` is exact by construction instead of checked to 1e-9. Same
  arithmetic as today (`at_last_node` is a bilinear read of exactly these values), so results are
  unchanged.
- **Field-name namespaces:** `Y_<species>` (mass fraction), `OMEGA_<species>` (species source at
  `RHO0`, `kg m-3 s-1`), suffixes `_PROD`/`_CONS` for split sources, `<Q>_MIN`/`<Q>_MAX` for
  bounds. Every other name is free.
- **Dimension scales** cost a few lines in the writer and make the file open as a labeled dataset in
  generic tools (expected to work with xarray through h5netcdf; to be verified in chemtab's tests,
  not promised in the spec). Mallard ignores them and reads `axes[]`.
- **No filters.** Mallard's HDF5 may be built without zlib, so compression would make a table
  build-dependent. Unused species cost disk only, since the reader no longer loads them.
- **Versioning rule** as written today (table_format.md §8): additions that readers can ignore
  need no new version.

Sizes, 200 x 20 x 200 nodes: a float64 field is 6.4 MB. Ten hot fields are 64 MB whatever the
species count; the current layout with 53 species is ~410 MB loaded on every rank and GPU.

### 2.2 Fields and tags

`fpv` (axes `Z`, `SZ`, `LAMBDA` with transforms `identity`, `normalized-variance of Z`, `bounds`
on `C`): required `C_MIN`, `C_MAX`, `OMEGA_C`, `RHO0`; transport `MU`, `KAPPA`, `CP`; optional
`HRR`. `SZ`'s last node may be below 1 (lookups clamp); `Z` and `LAMBDA` still span exactly [0, 1].

`saghafian2015`: required `T0`, `E0`, `R`, `GAMMA0`, `A_GAMMA`. New optional fields, written by
chemtab always and read by Mallard when it adopts them (decision 6): `A_MU`, `A_KAPPA_CP`
(exponents, units `1`: `mu = MU (T/T0)^A_MU`) and `T_A` (K: source factor
`exp(-T_A (1/T - 1/T0))`).

Reserved, specified when read: models `fgm`, `fpi`, `slfm`, `fpv-h`; thermodynamics `low-mach`.
With transforms in the file, `fpv-h` is "`fpv` plus an axis `H_NORM` with `transform = bounds` on
`H`" and needs no new mechanism, only its contract.

### 2.3 Why declared transforms

Every model on the roadmap maps transported quantities to a unit-box coordinate with the same
three operations: pass through, normalize a variance by `q (1 - q)`, normalize between tabulated
bounds that depend on the earlier axes. Declaring them per axis gives the generator one way to
describe FPV, FGM, heat loss, pressure and second streams, and gives the reader something to check
against its hand-written map. Mallard does **not** need a generic transform interpreter on the
device now: `FPV::locate` stays, and the reader verifies that the file's declared chain is the one
that code implements. A generic chain evaluator can come with the second model, if at all.

### 2.4 Axis lookup: a reader matter, not a format one

Analytic inverse maps in the file (the CharLES^X route) are the wrong trade for Mallard: a power-law
inverse costs a `pow` per lookup, and the cell must still be confirmed against the stored nodes to
keep the result bitwise equal to bisection. The reader can instead build, from the nodes alone, a
uniform bucket index per axis (bucket -> lowest overlapping cell, then at most one or two compares);
for FPVgen's grids that is 2,000-3,000 one-byte entries per axis. An axis that is uniform to rounding
(chemtab's `Lambda`) needs no index, only `floor` and the same fix-up compare. Both return exactly
the cell bisection returns, so `t` and the sum are bitwise unchanged, and an exhaustive test over
nodes and their neighbors in floating point proves it.

So the format carries `grid` (kind and parameters) as *descriptive* metadata only: for regeneration, plots
and the later CharLES^X exporter, which needs coordinate types. Whether to build the index at all
should wait for a profile; it is independent of the format.

### 2.5 What changes in Mallard's reader

| File | Change |
|---|---|
| `src/chemistry/flamelet_table.h` | `Table` holds named fields, each with its axes and dtype. New `pack(names, Value)` returns the node-major `TableView<D, Value>` of just those fields (the existing struct, `:47-128`, unchanged). `FPV` gains a `TableView<2>` of `[C_MIN, C_MAX]`; `locate` (`:154-164`) reads it instead of `at_last_node`. The test constructor (`:203`) takes fields with axes. |
| `src/chemistry/flamelet_table.cpp` | `read()` walks `/fields`, reading only contract fields plus requested outputs. The `mmap` path (`:395-440`) and the linearity check (`:245-262`) are deleted; `C_MAX >= C_MIN` and the transform-chain check are added. Contract tables (`:63-93`) and `is_name` (`:147`) take the all-caps names. |
| `src/physics/flamelet.cpp` | `make_gas` packs the solver's ten fields under their new names; the float32 widening copy (`:25-57`) is deleted, since packing converts. Output fields (`outputs`, `:220`) read host arrays per field. Warn (`W-TABLE-TRANSPORT`) unless `transport = "unity-lewis"`. |
| Kernels (`mixture.h`, `mixture_viscous_flux.h`, `solver_flamelet.cpp`) | None: they call `Gas::locate/thermo/transport/source`. |
| Tests and fixtures | `test/data/flamelet/*.h5` (6 files) rewritten by their h5py scripts; refusal cases in `test/flamelet_table_test.cpp` updated (drop "C not linear", add bounds and transform cases). |
| Docs | `docs/design/table_format.md` rewritten in place (still "version 1"); flamelet.md decision 4, §3, §4.2; guide. |

Core impact: confined to the two flamelet files and their tests; no solver, state or kernel change;
runs without `gas = "flamelet"` untouched. The packed record's layout (field order, padding to a
cache line, float32 or float64, e.g. thermo + transport as exactly one 64-byte line of doubles with
the two source fields beside it) becomes private to Mallard and tunable under GPU profiling without
touching chemtab. Device code is touched (the bounds view in `FPV`), so the GPU CI gate applies.
The CPU "read into memory for huge pages" finding (table_format.md §2) carries over: packing
allocates anonymous memory.

The conservative alternative, **keep the single `/values` array** (add reduced-rank side datasets
for bounds, required `transport`, optional Saghafian fields, writer rule "hot fields first, record
padded to 64 B"), fixes items 3, 6, 7 but leaves the writer responsible for a GPU layout it cannot
measure and keeps species in the loaded array. It is less work in the reader (about a third) and
worse everywhere else. With no released format to stay compatible with, nothing argues for it but
that saving.

---

## 3. chemtab

Designed from what an optimal generator is, not from FPVgen. FPVgen was read as a reference; the
verdict on it is in 3.2.

### 3.1 Design from first principles

**What a table generator is.** Four separable problems: (1) produce a family of 1-D
thermochemical states over a composition coordinate, labeled by a parameter; (2) turn each state
into the fields a solver reads; (3) integrate presumed PDFs; (4) choose a tensor grid, remap onto
it and write. Only (1) is model-specific and expensive. Everything in chemtab hangs off that split.

**Flamelets in mixture-fraction space, not physical space.** The FPV model *is* the steady
flamelet equations, `rho chi(Z)/2 d2Y_k/dZ2 + omega_k = 0`, with unity Lewis number and total
enthalpy linear in `Z`. Solving them directly, on a fixed `Z` grid with `chi_st` as the parameter
and pseudo-arclength continuation through the turning points, is the right primary source:

- The solution is already on `Z`: no interpolation from `x`, no monotone-`Z(x)` assumption, no
  domain width, strain-to-`chi_st` fit or control points. Those four are most of FPVgen's solver
  code (`FPVgen:224-601`) and all of its fragility.
- Arclength continuation traces stable, unstable and (for hydrocarbons with low-temperature
  chemistry) further branches as one curve, with the step set by the curve, not by a temperature
  decrement tuned per fuel (`sample/input.toml:36-41`).
- Unity Lewis number holds exactly, which is what Mallard's flamelet gas assumes
  (flamelet.md decision 8), instead of approximately through a transport-model switch.
- Heat loss is a boundary condition on `h(Z)`, a pressure sweep is a parameter, and the unsteady
  flamelets needed later for NO and autoignition are the same operator with a time term. The
  roadmap's axes come from one solver.
- It is small: with `h` linear in `Z` only the species are unknown; the Jacobian is block
  tridiagonal and Cantera supplies rates, thermodynamics and rate derivatives. A few hundred lines
  around a banded solve.

Cost and risk, stated plainly: chemtab then owns a 1-D solver that must be verified, and its
flamelets differ slightly from a physical counterflow flame's because `chi(Z)` is prescribed
(the usual `erfc` form) rather than computed. This is unproven here; a spike on the H2/N2 sample
and on a methane mechanism should come before committing (decision 9). Cantera's
`CounterflowDiffusionFlame` stays in chemtab as an independent reference in the tests, and as the
source for differential-diffusion studies, where the `Z`-space equations are no longer simple.

**Sources by model.** `fpv`, `slfm`, `fpv-h`: the `Z`-space solver. `fgm`, `fpi`: Cantera's free and
burner-stabilized premixed flames (nothing to gain by owning those). Mixing line, equilibrium and
homogeneous reactors: 0-D Cantera. All behind one interface that returns states on a coordinate
with a label, so stages 2 to 4 never know which ran.

**Error-controlled tabulation.** The lookup is multilinear, so its error against the fine manifold
is computable before writing. chemtab places the nodes of each axis to equidistribute that error
and reports the largest interpolation error per field. The user gives a tolerance (default) or
node counts (expert). This replaces hand-tuned stretching (FPVgen: a third of the `Z` nodes below
`Z_st`, a power law of 2.7 in variance, `FPVgen:1321-1323`) and makes "never lose accuracy" a
number in the report rather than a hope.

**Exact where exactness is available.** Beta-PDF means by incomplete-beta weights
(flamelet.md §2), thermodynamic coefficients in closed form from the NASA polynomials
(`GAMMA0 = cp/cv`, `A_GAMMA = -R (dcv/dT) / cv^2`), bounds instead of a linearity check.

**A library first.** `chemtab.build(config) -> Table` is the product; the CLI is a thin layer, and
a `Table` is plain numpy arrays with named axes, so notebooks and other writers need nothing else.

### 3.2 FPVgen: what to keep and what to throw away

Kept, as ideas, because they are right:

- A resumable on-disk store of solved flamelets, so tabulation never re-solves (`FPVgen:1008-1151`).
- Clustering `Z` nodes around `Z_st` (`coordinate.py:215-260`), as the starting guess of the
  error-controlled grid, and grids that record their kind and parameters (`coordinate.py:73-85`).
- The complete Saghafian field set, including the transport exponents and the activation
  temperature (`FPVgen:1452-1465`, `:1524-1529`), which Mallard's format had dropped.
- TOML input, and the diagnostic plots as a product: S-curve, profiles over `Z`, table slices
  (`FPVgen:1683-1987`).
- The mixing solution closing the table at `C = 0` (`FPVgen:873-1006`), computed in 0-D instead.

Thrown away:

- The physical-space S-curve machinery as the table source: two-point control with adaptive
  `delta_T` and error back-off (`:651-871`), control-point placement (`:532-601`), flame-width
  adaptation (`:409-531`), the strain/`chi_st` fit (`:224-348`, `:1655`). It is FPVgen's best
  engineering and it solves a problem the `Z`-space formulation does not have.
- The whole assembly (`:1289-1653`): the variance axis tiled from the laminar column (`:1598`, with
  2 variance nodes by default, `sample/input.toml:47`); `Lambda` from flamelets sorted by peak `C`,
  normalized per `Z` with 1e-4 tolerances, a uniform fallback, a `force_monotonicity` patch and
  extrapolation (`:1531-1596`); Saghafian coefficients by finite differences of a +-5000/rho J/kg
  perturbation in a Python loop over grid points (`:1391-1445`); `interp1d` over `Z(x)`
  (`:1357-1358`); sources per unit mass (`:1468`, `:1478`, `:1491`).
- The CharLES^X and FlameMaster writers (`:1153-1287`, `:1611-1653`). A CharLES^X exporter, when
  wanted, is written against chemtab's `Table`.
- The structure: one 1,987-line class holding solve, storage, assembly, writers and plots
  (`FPVgen:40`), no tests (`tests/` holds only `__init__.py`), console scripts in a top-level
  package named `scripts` (`pyproject.toml:44-47`).

Net: under this design chemtab shares **no code** with FPVgen, so the README says "inspired by
FPVgen", with its citation. If decision 9 instead keeps the Cantera continuation as the table
source and ports it, the README says "derived from FPVgen" and carries its MIT notice.

### 3.3 Architecture

```
chemtab/
  api.py            build(config) -> Table; the stages below as functions
  cli.py            init | run | check | info | plot
  config.py         TOML -> frozen dataclasses; every key read is recorded, unused keys are errors
  mechanism.py      mechanism + SHA-256, streams, Bilger Z, Z_st, progress variable
  sources/
    base.py         FlameletSource: solve() -> states on a coordinate, each with a label
    zspace.py       steady flamelet equations in Z, pseudo-arclength continuation in chi_st
    premixed.py     later: Cantera free and burner-stabilized flames (FGM/FPI)
    zero_d.py       mixing line, equilibrium, homogeneous reactors
    counterflow.py  Cantera counterflow flames: test reference; source for differential diffusion
  library.py        resumable HDF5 store of solved states, keyed by a hash of what determines them
  manifold.py       states -> fields on the fine coordinate: thermo coefficients, transport, sources
  pdf.py            exact beta weights, delta; later beta in C, top-hat
  remap.py          label -> normalized axis with monotonicity repair and report
  grid.py           error-controlled node placement; fixed kinds (uniform, power, piecewise, log)
  table.py          Table: axes, fields with their axes, tags, model data; invariants
  io/mallard.py     writer, reader, and check(): the spec's "what Mallard checks", independently
  validate.py       interpolation error, off-table flamelets, PDF vs quadrature; the report
  plot.py           S-curve, profiles, table slices (matplotlib is an optional extra)
```

Pipeline, each stage a pure function of its inputs with its output cached under a content hash, so
changing the tolerance or the species list never re-solves flamelets:

```
mechanism + streams -> source -> library of states phi(Z; label) -> manifold fields
  -> PDF integration phi(Z, SZ, label) -> grid + remap to the table's axes -> check -> write
```

A model is a small class naming its source, label, axes with transforms, PDF per axis and
contract fields. FGM swaps the source; heat loss adds an enthalpy parameter to the same source
and a `bounds` axis; an equilibrium table uses `zero_d` alone.

Determinism: solves are serial and repeatable; parallelism is over independent continuations
(pressures, enthalpy levels) in processes. Two builds from the same library give bitwise equal
fields (`chemtab check --same a.h5 b.h5` ignores only `created`).

### 3.4 Config and CLI

TOML in Mallard's style: lower-case sections, SI units, strict keys, defaults documented in one
generated page. The whole beginner input:

```toml
[mechanism]
file = "h2o2.yaml"                  # any Cantera mechanism; bundled names resolve
transport = "unity-lewis"           # default; what Mallard's flamelet gas assumes

[streams.fuel]
X = { H2 = 0.5, N2 = 0.5 }
T = 300.0

[streams.oxidizer]
X = { O2 = 0.21, N2 = 0.79 }
T = 300.0

[model]
type = "fpv"
pressure = 101325.0
progress_variable = { H2O = 1.0 }   # default chosen from the fuel's elements if omitted

[table]
output = "h2_fpv.h5"
tolerance = 1e-3                    # largest relative interpolation error; or nodes = { Z = 201, SZ = 21, LAMBDA = 101 }
species = ["H2O", "OH"]             # or "all"; costs disk only
```

Expert keys live in `[flamelets]` (source, continuation steps, `Z` resolution, branch limits, restart)
and `[table.grid.<axis>]` (kind, clustering), all defaulted. A first table in three commands:

```bash
pip install chemtab
chemtab init h2-air            # writes chemtab.toml from a preset
chemtab run chemtab.toml       # solve, build, check, write h2_fpv.h5 and report.html
```

then `table = "h2_fpv.h5"` in Mallard. `run` prints progress along the S-curve (branch, `chi_st`,
`T_max`), ends with the checks it passed and the off-table error, and on failure says which stage
and which key to change. `chemtab info table.h5` prints tags, axes, fields and provenance;
`chemtab check` validates any file against the spec; `chemtab plot` draws the S-curve and slices.

### 3.5 Testing and validation

Each test must fail under a plausible bug (the admission rule Mallard uses).

- **Flamelet solver:** residual of every stored solution below tolerance on a refined `Z` grid;
  the analytic Jacobian against finite differences; one-step chemistry with a known S-curve
  (turning points in closed form); the stable branch and the extinction `chi_st` against Cantera's
  counterflow flames within a stated band (they differ by the shape of `chi(Z)`, so this is a
  consistency check, not an identity).
- **Units of the pipeline, no chemistry:** a synthetic flamelet family (Burke-Schumann
  profiles with an analytic `C_st`) drives manifold, PDF, remap and writer in seconds. Beta weights
  against high-precision quadrature to 1e-9, including `Sz -> 0`, `Sz = 1` and means at 0 and 1;
  the remap returns `C` exactly at the bounds and reports non-monotone columns; closed-form
  `a_gamma` against the finite difference; units of every written field.
- **Contract:** write then read back bitwise; `chemtab check` refuses each violation Mallard's
  reader refuses (the list in flamelet.md §9), and reads Mallard's reference fixtures at a pinned
  tag; Mallard reads a chemtab-written table in its own CI. Two independent implementations of
  the checks is the point.
- **Regression on the documentation's sample** (H2/N2 1:1 against air, 300 K, 1 atm, `h2o2.yaml`,
  unity Lewis, `C = Y_H2O`). The documented numbers are Cantera counterflow results: extinction at
  `chi_st` = 2777 1/s and `T_max` = 1246 K, stable branch from 0.32 1/s, unstable branch back to
  0.52 1/s at 986 K and `C_st` = 0.059 (flamelet.md §2). The `cantera-counterflow` source must
  reproduce them; the `Z`-space source gets its own pinned values, and the difference between the
  two is reported in the docs rather than hidden.
- **Validation targets** (flamelet.md §10): an off-table flamelet against a direct solve of the
  same source, `T` within 0.5% of the peak; the table's interpolation error against the fine
  manifold below the requested tolerance; PDF integration to 1e-9. Reported in every `run`.
- **Deliverable to Mallard** (flamelet.md §5.5): the committed H2/N2 table with `HRR`, `Y_OH`,
  `Y_H2O` for `examples/counterflow_diffusion`, a few MB, with its input in `/provenance/input`.

### 3.6 Packaging

New repository (not a GitHub fork), MIT. The README credits FPVgen and cites its DOI: as
"inspired by" under the recommended design, which shares no code with it, or as "derived from"
(keeping its MIT notice) if decision 9 keeps FPVgen's Cantera continuation as the table source. PyPI
`chemtab`, import `chemtab`, one console script `chemtab`. Python >= 3.10; dependencies `cantera >=
3.1`, `numpy`, `scipy`, `h5py`; extras `plot` (matplotlib) and `dev`. `uv` for
development, `ruff`, `mypy`, `pytest`. Versioned independently of Mallard: the table's `version`
is the coupling.

---

## 4. Joint decisions

Already decided (section 2.0), not open: the format is revised in place as `version = 1` with no
compatibility path; field names are all-caps, uppercased faithfully (`RHO0`, `T0`, `OMEGA_C`, ...).

Open, ordered so that each depends only on those before it.

1. **File layout.**
   (a) Current: one node-major `/values`, the file is the lookup layout.
   (b) Current plus additions: side datasets for bounds, writer rule "hot fields first, padded".
   (c) One dataset per field with its own axes; Mallard packs its own device record.
   Tradeoffs: (a) zero reader work, species bloat in the hot path, writer owns a GPU layout it
   cannot measure. (b) small reader change, still couples the two. (c) rewrites `Table::read` and
   the fixtures, deletes the `mmap` and widening paths, and makes the GPU layout Mallard-private.
   **Recommend (c).**

2. **Progress variable: stored `C` or bounds, and declared axis transforms.**
   (a) 3-D `C` with the linearity check (current). (b) `C_MIN`, `C_MAX` on `(Z, SZ)` plus per-axis
   `quantity`/`transform` attributes, checked against Mallard's hand-written map. (c) as (b) with a
   generic transform interpreter on the device.
   Tradeoffs: (b) is exact by construction, one field smaller, the map reads a 64 KB table instead
   of 8 scattered nodes, and heat loss, FGM and pressure need no new mechanism; (c) is premature.
   **Recommend (b).** Requires 1(b) or 1(c).

3. **Axis names: do they follow the all-caps rule?**
   (a) `Z`, `SZ`, `LAMBDA` (and later `H_NORM`, `P`). (b) Keep FPVgen's `Z`, `Sz`, `Lambda`.
   Tradeoffs: (a) is one rule for every name a user types or greps, and reduced fields named after
   an axis quantity (`C_MIN`) read consistently; prose and equations can still write `Sz` and
   `Lambda` as symbols. (b) keeps the literature's look but leaves two conventions in one file.
   `LAMBDA` does not collide with Mallard's output `LAMBDA2`, and axes are not output fields
   (the normalized progress variable is output as `C_NORM`).
   **Recommend (a),** with structural names (groups, attributes) lower-case like TOML keys.

4. **O(1) axis lookup.**
   (a) Keep bisection. (b) Normative analytic mappings in the file. (c) Reader-built bucket index,
   bitwise equal to bisection; file carries grid kind only as description.
   Tradeoffs: (b) puts `pow` on the device, constrains chemtab's grids and risks cell disagreement;
   (c) needs no format support and any axis works.
   **Recommend (c), implemented only if a profile shows `locate` matters;** format carries
   descriptive `grid` metadata either way.

5. **Precision.**
   (a) File dtype decides the device dtype. (b) A Mallard input key decides, independent of the
   file. (c) float64 only.
   Tradeoffs: today float32 files are widened on the device, so the option has no runtime effect.
   (b) lets a user narrow a float64 table silently, against "never lose accuracy".
   **Recommend (a):** chemtab writes float64 by default with a per-table `precision` key; Mallard's
   packed record takes the dtype of the fields it packs and logs it; a true float32 device path is
   built when lookups show in a profile.

6. **Completing Saghafian's model: `A_MU`, `A_KAPPA_CP`, `T_A`.**
   (a) Leave out. (b) Optional fields now, always written by chemtab, adopted by Mallard in a
   later, separately validated change. (c) Required now and used at once.
   Tradeoffs: without them transport and the source ignore the difference between local and
   flamelet temperature; (c) couples the format change to a physics change that needs its own
   validation.
   **Recommend (b).** A physics call for the user: whether Mallard's `(rho/RHO0)^2` source scaling
   also takes the Arrhenius factor.

7. **More dimensions: pressure, heat loss, second stream.**
   (a) Extra axes of one tensor table. (b) Several tables per file, the solver interpolating
   between them. (c) Several files.
   Tradeoffs: (a) is uniform with everything else and memory grows by the axis length (fine for 3
   to 5 pressures); (b), (c) add a second interpolation path and split provenance.
   **Recommend (a),** one table per file; root `pressure` is "required unless a `P` axis exists";
   Mallard instantiates `TableView<4>` when it reads such a model.

8. **Home of the specification and conformance.**
   (a) Mallard owns `table_format.md` and reference fixtures; chemtab implements an independent
   `check` and tests against Mallard's fixtures at a pinned commit. (b) chemtab owns the spec and
   Mallard links to it. (c) a third repository.
   Tradeoffs: the reader is what refuses a table, so it is the authority; (b) lets generator
   releases move the contract.
   **Recommend (a),** format name stays `mallard-table`. Optional small addition on Mallard's side:
   a `mallard table check <file>` command so the real reader is available as a validator.

9. **chemtab's flamelet source, and with it its relation to FPVgen.**
   (a) Cantera counterflow flames with two-point control as the table source, ported from FPVgen
   (README: "derived from FPVgen"). (b) The same, written new from the 60-line continuation of
   `tools/docs/flamelet.py:106-140` ("inspired by"). (c) chemtab's own `Z`-space flamelet solver
   with arclength continuation as the source; Cantera counterflow only as a test reference and for
   differential diffusion ("inspired by"; no shared code).
   Tradeoffs: (a) and (b) need no solver verification and match Mallard's counterflow validation
   case by construction, but inherit physical-space fragility, per-fuel tuning, the `x -> Z`
   interpolation and only approximate unity-Lewis consistency; heat loss and unsteady flamelets
   need further machinery. (c) is the model's own equations, robust, fast and extensible, at the
   price of owning and verifying a solver, and of tables that differ slightly from counterflow
   flames through the prescribed `chi(Z)`.
   **Recommend (c), gated on a short spike** (H2/N2 sample and one methane mechanism: full S-curve,
   extinction `chi_st` against Cantera's). If the spike fails, (b). Everything downstream of the
   source (manifold, PDF, grid, writer) is new in every option; FPVgen's assembly is not kept.

10. **Transport default and its place in the contract.**
    (a) chemtab defaults to unity Lewis; `transport` is a required attribute; Mallard warns on
    anything else. (b) Mallard refuses non-unity tables. (c) Leave it as optional provenance.
    Tradeoffs: differential-diffusion tables are a legitimate study (flamelet.md decision 8), so
    refusing is too strong; silence is the current hazard.
    **Recommend (a).**

11. **The documentation's sample table.**
    (a) Keep `tools/docs/flamelet.py` as a private 467-line generator for figures. (b) Once chemtab
    reproduces its numbers, commit the chemtab-made table and reduce the script to plotting it;
    the same file serves the counterflow validation.
    Tradeoffs: (b) removes a second generator from Mallard and makes the docs show what users will
    actually get; it adds no build dependency since the table is committed.
    **Recommend (b),** after chemtab's regression test passes.
