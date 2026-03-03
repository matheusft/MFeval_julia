# MFeval.jl

A high-performance Julia implementation of the **Pacejka Magic Formula** tyre
model (MF 5.2, 6.1 and 6.2), reimplemented from the MATLAB
[mfeval](https://uk.mathworks.com/matlabcentral/fileexchange/63618-mfeval)
toolbox by Marco Furlan.

> **Reference equations**  
> Pacejka, H.B. — *Tyre and Vehicle Dynamics*, 3rd ed., Elsevier, 2012  
> Besselink et al. — *An improved Magic Formula/Swift tyre model that can handle
> inflation pressure changes*, Vehicle System Dynamics, 48:1, 2010.
> DOI: [10.1080/00423111003748088](https://doi.org/10.1080/00423111003748088)

---

## Project structure

```
MFeval.jl/
│
├── Project.toml               ✅  Dependencies: StaticArrays, Parameters
├── README.md
│
├── src/
│   ├── MFeval.jl              ✅  Package root — exports, module imports
│   │
│   ├── types/
│   │   ├── TireParams.jl      ✅  TireParams{Version} immutable struct + constructor
│   │   ├── MFInputs.jl        ✅  MFInputs struct (Fz, kappa, alpha, gamma, phit, Vx, P, omega)
│   │   ├── MFOutputs.jl       ✅  MFOutputs struct (30 fields, mirrors mfeval columns)
│   │   └── MFModes.jl         ✅  MFModes struct (useLimitsCheck, useAlphaStar, useTurnSlip)
│   │
│   ├── io/
│   │   └── read_tir.jl        ✅  TIR file parser → TireParams
│   │
│   ├── solver/
│   │   ├── parse_inputs.jl    ⬜  Limit clamping, low-speed smoothing, alpha_star, pressure
│   │   ├── basic_vars.jl      ⬜  Shared pre-calculations (dfz, dpi, star/prime vars, slip velocities)
│   │   ├── Fx0.jl             ⬜  Pure longitudinal force
│   │   ├── Fy0.jl             ⬜  Pure lateral force
│   │   ├── combined_slip.jl   ⬜  Gxα and Gyκ weighting functions
│   │   ├── Mx.jl              ⬜  Overturning moment
│   │   ├── My.jl              ⬜  Rolling resistance moment
│   │   ├── Mz.jl              ⬜  Self-aligning moment (trail + residual torque)
│   │   ├── turn_slip.jl       ⬜  ζ factors for turn slip (φ, ψ_dot)
│   │   ├── geometry.jl        ⬜  Re, Rl, ρ, contact patch (a, b), Cz
│   │   └── relaxation.jl      ⬜  Cx, Cy, σx, σy, instantaneous Kya
│   │
│   └── api/
│       ├── mfeval.jl          ⬜  Main scalar entry point — calls full pipeline
│       └── mfeval_batch.jl    ⬜  Batch (matrix) entry point with Threads.@threads
│
└── test/
    ├── runtests.jl            ✅  Test runner — include all test files here
    ├── test_phase1.jl         ✅  Types, I/O, struct correctness
    ├── fixtures/
    │   ├── MagicFormula52_Parameters.tir   ✅
    │   ├── MagicFormula61_Parameters.tir   ✅
    │   └── MagicFormula62_Parameters.tir   ✅
    ├── test_mf52.jl           ⬜  MF 5.2 specific solver cases
    ├── test_mf61.jl           ⬜  MF 6.1 specific solver cases
    ├── test_mf62.jl           ⬜  MF 6.2 + iterative Rl solver
    ├── test_turn_slip.jl      ⬜  Turn slip edge cases
    ├── test_low_speed.jl      ⬜  Low-speed smoothing and singularity protection
    ├── test_vs_matlab.jl      ⬜  Numerical regression vs. MATLAB mfeval (tol < 1e-6)
    └── benchmarks.jl          ⬜  BenchmarkTools suite
```

✅ complete · ⬜ planned

---

## Running the tests

### Prerequisites

Install Julia 1.9+ from [julialang.org](https://julialang.org/downloads/), then
instantiate the project once:

```bash
cd MFeval.jl/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run the full test suite

```bash
# From the MFeval.jl/ root directory
julia --project=. test/runtests.jl
```

Expected output (Phase 1):

```
============================================================
Phase 1 — Types, I/O and structs
============================================================
Test Summary:    | Pass  Total  Time
MFModes decoding |   15     15  0.0s
Test Summary:         | Pass  Total  Time
MFInputs constructors |   12     12  0.0s
Test Summary:        | Pass  Total  Time
MFOutputs round-trip |    1      1  0.0s
Test Summary:                   | Pass  Total  Time
TireParams default construction |   21     21  0.1s
Test Summary:  | Pass  Total  Time
read_tir MF6.1 |   18     18  3.0s
Test Summary:  | Pass  Total  Time
read_tir MF5.2 |    2      2  0.0s
Test Summary:  | Pass  Total  Time
read_tir MF6.2 |    3      3  0.0s
Test Summary:           | Pass  Total  Time
read_tir error handling |    1      1  0.0s
Test Summary:                | Pass  Total  Time
TireParams type concreteness |  269    269  0.0s
```

### Run a single test file

```bash
julia --project=. test/test_phase1.jl
```

### Run from the Julia REPL

```julia
julia> using Pkg
julia> Pkg.activate(".")       # activate the project
julia> Pkg.instantiate()       # first time only
julia> Pkg.test()              # runs test/runtests.jl
```

### Verify type stability (no `Any` fields)

Type stability is a hard requirement — a single `Any` field disables
specialisation for that code path. Check it with:

```julia
julia> using MFeval
julia> p = read_tir("test/fixtures/MagicFormula61_Parameters.tir")
julia> @code_warntype MFModes(111)           # should show all Bool
julia> @code_warntype read_tir("test/fixtures/MagicFormula61_Parameters.tir")
```

No yellow/red output means the compiler has fully resolved all types.

---

## Development phases

### Phase 1 — Foundation ✅

All data structures and the TIR file parser.

- **`TireParams{V}`** — a fully concrete, immutable struct with 270 `Float64`
  fields. The type parameter `V ∈ {:MF52, :MF61, :MF62}` encodes the Magic
  Formula version at compile time, resolved from `FITTYP` by `read_tir`. This
  means every solver function that dispatches on `TireParams{V}` generates a
  fully specialised, branch-free compiled method — there is no `if FITTYP == 61`
  in the hot path.
- **`TireMetadata`** — string fields (TYRESIDE, unit labels, etc.) separated
  from the arithmetic struct so they never affect memory layout or JIT
  specialisation.
- **`MFInputs`** — named scalar struct for one row of the input matrix, with
  6-, 7- and 8-argument convenience constructors.
- **`MFModes`** — decodes the 3-digit `useMode` integer once at construction,
  giving three plain `Bool` fields to the solver.
- **`MFOutputs`** — 30 `Float64` fields in MATLAB column order; `to_vector` and
  `from_matrix_row` provide lossless interop with the matrix-based API.
- **`read_tir`** — robust two-pass parser: first pass builds a
  `Dict{String,Any}`, second pass calls the `TireParams` keyword constructor.
  Missing keys fall back to MATLAB-equivalent defaults. Handles inline `$`
  comments, the obsolete `[SHAPE]` section, and the `MASS`/`MASS1` rename.

### Phase 2 — Scalar solver core ⬜

Each calculation module is a standalone `@inline @fastmath` function taking
only `Float64` scalar arguments. No arrays, no heap allocations, no runtime
version branches (dispatched via `TireParams{V}`).

| File | Equations |
|---|---|
| `parse_inputs.jl` | Input clamping, low-speed cosine/linear reductions, alpha_star, pressure delta |
| `basic_vars.jl` | dfz, dpi, Vcx, LMUX_star/prime, slip velocities |
| `Fx0.jl` | Pure longitudinal force (4.E1–4.E19) |
| `Fy0.jl` | Pure lateral force (4.E20–4.E30) |
| `combined_slip.jl` | Gxα, Gyκ weighting functions (4.E54–4.E65) |
| `Mx.jl` | Overturning moment (Eqn 49 from Besselink paper) |
| `My.jl` | Rolling resistance moment (4.E67–4.E68) |
| `Mz.jl` | Self-aligning moment — trail + residual torque (4.E31–4.E50) |
| `turn_slip.jl` | ζ1–ζ8 reduction factors (4.E81–4.E107), only when `useTurnSlip=true` |
| `geometry.jl` | Re, Rl, ρ, contact patch a/b, Cz; iterative secant solver for MF6.2 Rl |
| `relaxation.jl` | Cx, Cy, σx, σy, instantaneous Kya |

### Phase 3 — Public API ⬜

- `mfeval(params, input::MFInputs, mode::MFModes) → MFOutputs` — scalar entry
  point, calls the full pipeline, returns a named struct.
- `mfeval(params, inputs::Matrix, mode) → Matrix{Float64}` — batch entry point,
  pre-allocates an `N×30` output matrix, processes rows with `Threads.@threads`.

### Phase 4 — Validation ⬜

- Generate reference outputs from MATLAB mfeval across a full parameter sweep
  for all three TIR fixtures.
- Regression tests: all 30 output columns within tolerance `< 1e-6`.
- Edge cases: `Fz = 0`, `Vx = 0` (low-speed), locked wheel (`κ = −1`), max
  camber, inputs outside valid range.

### Phase 5 — Performance ⬜

- Benchmark scalar path — target `< 1 µs` per evaluation.
- Confirm `@code_warntype` shows zero `Any` / type instability across the full
  hot path.
- Verify linear thread scaling for the batch entry point.
- Profile and annotate `@simd` where beneficial (e.g. inner loops in contact
  patch geometry).

---

## Key design decisions

| Decision | Rationale |
|---|---|
| `TireParams{V}` type parameter | MF version resolved at compile time — zero runtime branching in solver |
| All 270 coefficients as `Float64` fields | Flat struct in memory, no pointer chasing, fully concrete for JIT |
| `TireMetadata` separate from `TireParams` | String fields never pollute the arithmetic struct's memory layout |
| `@inline @fastmath` on all solver kernels | Fused multiply-add, scalar reassociation; safe for tyre model arithmetic |
| Immutable structs for all I/O | Stack-allocated, zero GC pressure |
| `MFModes` decoded once at construction | Decoding cost paid once; solver sees plain `Bool` fields |
| Secant solver for MF6.2 Rl | Simple scalar loop — no closures, no allocations |
| `Threads.@threads` on batch loop | Embarrassingly parallel: each row is fully independent |
| Pre-allocated output matrix | Caller provides the buffer; zero allocation per batch call |