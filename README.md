MATLAB implementations and reproducibility scripts accompanying the manuscript:

> **Order-preserving safeguards and sharp convergence bounds for Matérn RBF time stepping**  
> Roberto Cavoretto and Alessandra De Rossi

The package investigates Taylor-expanded Matérn radial basis function (MA-RBF) corrections to fixed-step methods for ordinary differential equation initial-value problems. It includes classical and multiquadric RBF (MQ-RBF) comparisons, the original MA-RBF recurrences, and order-preserving safeguarded Adams methods for scalar problems and systems.

The scripts cover **Tables 1-4, Figure 1(a-c), and the additional numerical experiments described in the manuscript's text**. Each numbered table or figure has an independent driver with all parameters already specified. Here, *adaptive* refers to the step-dependent algebraic correction, not to adaptive time stepping.

## Quick start

Extract the complete package, preserving its subfolders. In MATLAB, set **Current Folder** to the directory containing `setup_package.m`, then enter:

```matlab
run_core
```

This runs the consistency checks and generates all four numbered tables and the complete three-panel figure. Only base MATLAB is needed for this workflow, subject to the version requirements below.

To generate just one item, run its script directly, for example:

```matlab
script_table1
```

Each driver calls `setup_package` automatically. There is no need to run a setup script first, add subfolders manually, edit parameters, or generate another output beforehand. The package directory must be writable.

For all numbered outputs, all textual experiments, and the additional 80-digit checks, use:

```matlab
run_all
```

**`run_all` requires Symbolic Math Toolbox.** It checks for the required symbolic functions before starting the experiments and does not silently substitute double precision.

## Software requirements

| Component | Requirement |
| --- | --- |
| MATLAB | R2017b or later is the source-based compatibility baseline. The stability figure uses `polyshape`, introduced in R2017b, and the drivers contain local functions in scripts. This is not a claim that every release has been tested. |
| Base-MATLAB workflow | Sufficient for `run_core`, `verify_package`, every numbered table/figure script, and the double-precision textual experiments listed below. |
| Symbolic Math Toolbox | Required by `script_text_order_preserving.m`, `script_text_sharp_asymptotics.m`, `script_verify_80digits.m`, and therefore `run_all.m`. These scripts request 80 significant decimal digits using `vpa` and `digits`. |
| LaTeX, optional | Needed only to compile the generated table code. The receiving document must load `booktabs` and `adjustbox`. LaTeX is not required to run the MATLAB experiments or export figures. |

The supplied drivers do not call third-party MATLAB packages or require external numerical datasets. MATLAB graphics support is required for figure export. GNU Octave compatibility is not established.

The release and precision requirements are based on the source and the MathWorks documentation for [polyshape](https://www.mathworks.com/help/matlab/ref/polyshape.html), [local functions in scripts](https://www.mathworks.com/help/matlab/matlab_prog/local-functions-in-scripts.html), and [variable-precision arithmetic](https://www.mathworks.com/help/symbolic/sym.vpa.html).

## Folder structure

```text
package-root/
|-- README.md
|-- setup_package.m
|-- run_core.m
|-- run_all.m
|-- verify_package.m
|-- script_table1.m
|-- script_table2.m
|-- script_table3.m
|-- script_table4.m
|-- script_stability_regions.m
|-- script_text_*.m
|-- script_verify_80digits.m
|-- table_latex.m
|-- Classical methods for IVP/
|   `-- Euler, AB2, AB3, and AB2-AM1 reference methods
|-- RBF methods for IVP/
|   `-- Original RBF formulas, MA_adams, and implicit-step solver
|-- Stability regions/
|   `-- Formal-target stability data and plotting routines
|-- Reproducibility/
|   `-- Test problems, error measures, checks, and output helpers
`-- output/                       [created automatically]
    |-- tables/                   LaTeX tables and text snippets
    |-- figures/                  PDF, EPS, and PNG exports
    |-- data/                     Computed CSV/MAT data and metadata
    `-- checks/                   Timing and verification reports
```

The three method/stability directories retain the original package organization. Shared experiment and reporting utilities are collected in `Reproducibility/`. Keep the directory names, including spaces, unchanged.

## Reproduce the numbered tables and figure

All paths in the following table are relative to the package root. Run the script in the middle column; no parameters need to be supplied.

| Manuscript item | MATLAB script | Generated output |
| --- | --- | --- |
| Table 1(a,b): Example 1, Euler and standalone AM1 comparisons | `script_table1.m` | `output/tables/table1.tex`; data stems `table1a` and `table1b` |
| Table 2(a,b): Example 2, unsafeguarded AB2 and AB2-AM1 comparisons | `script_table2.m` | `output/tables/table2.tex`; data stems `table2a` and `table2b` |
| Table 3(a,b): Example 3, active original safeguards | `script_table3.m` | `output/tables/table3.tex`; data stems `table3a` and `table3b` |
| Table 4: Example 3, nondegenerate control interval | `script_table4.m` | `output/tables/table4.tex`; data stem `table4` |
| Figure 1(a-c): formal-target stability diagnostics | `script_stability_regions.m` | `output/figures/figure1.*`, `figure1a.*`, `figure1b.*`, `figure1c.*`; `output/data/figure1_stability_data.mat` |

Each numerical data stem denotes both `output/data/<stem>.csv` and `output/data/<stem>.mat`. For figures, `*` denotes `pdf`, `eps`, and `png`: the script exports vector PDF/EPS files and 600-dpi PNG files. It creates and closes figures automatically, so the absence of an open figure window does not indicate failure.

### Settings used by the numbered drivers

| Output | Problem and interval | Meshes and principal settings |
| --- | --- | --- |
| Table 1 | Example 1: `u' = -u^2`, `u(0) = 1`, on `[0,1]` | `N = 10,20,...,1280`; classical, MQ-RBF, and MA-RBF Euler and standalone AM1; endpoint errors. |
| Table 2 | Example 2: `u' = (2t^2-u)/(t^2u-t)`, `u(1) = 2`, on `[1,3]` | `N = 10,20,...,320`; classical, MQ-RBF, and MA-RBF AB2 and single-correction AB2-AM1; no safeguard; endpoint errors. |
| Table 3 | Example 3: `u' = t(1+t)u^2`, `u(-1/2) = 12/11`, on `[-1/2,1/2]` | `N = 128,256,...,8192`; original MA-RBF AB2/PC safeguards with `delta_f = 0.1,0.2`; maximum nodal errors and activation counts; unsafeguarded failure checks. |
| Table 4 | Same equation as Example 3, restarted at `u(1/4) = 192/185`, on `[1/4,1/2]` | `N = 16,32,64,128,256`, `h = 1/(4N)`; both fixed thresholds and unsafeguarded comparisons; maximum nodal errors and zero activations. |

Tables 1-4 concern the **original recurrences**, not the order-preserving reformulations. Those reformulations are tested by the additional scripts below.

**Figure 1 is a formal-target diagnostic.** Its panels compare Euler, AB2, and AM1 formulas after substitution of the formal correction targets. They are not stability regions of the actively safeguarded, history-dependent maps. In panel (c), the MQ-RBF and MA-RBF AM1 boundaries coincide. The plotting grids and amplification factors are defined in `Stability regions/formal_stability_data.m`.

## Additional experiments reported in the text

These scripts address numerical results discussed in Section 5 but not assigned separate numbered tables or figures. Each can be run independently.

| Manuscript experiment | MATLAB script | Arithmetic | Output stems |
| --- | --- | --- | --- |
| Example 1: same-order AB2 reference versus MA-RBF Euler | `script_text_same_order.m` | Double | `text_same_order` |
| Example 1: two RK4 startup steps versus exact starts | `script_text_startup.m` | Double | `text_startup` |
| Example 3: unsafeguarded even/odd mesh robustness | `script_text_robustness.m` | Double | `text_robustness` |
| Example 3: fixed-threshold order preservation and AB3 reference | `script_text_order_preserving.m` | 80 digits | `text_order_preserving` |
| Scalar forcing: sharp logarithmic coefficients, endpoint cancellation, and direct-recurrence identities | `script_text_sharp_asymptotics.m` | 80 digits | `text_sharp_asymptotics`, `text_sharp_identity` |
| Scalar forcing: shrinking thresholds `delta_f = h^alpha`, `alpha = 1/4,1/2,3/4` | `script_text_shrinking_thresholds.m` | Double; compensated residual sums | `text_shrinking_thresholds` |
| Vector forcing: logarithmic loss and tangent-projection cancellation | `script_text_vector_forcing.m` | Double; compensated residual sums | `text_vector_forcing` |
| Coupled system: order preservation and rotation/translation equivariance | `script_text_coupled_system.m` | Double | `text_coupled_system` |

Every listed stem produces CSV/MAT files in `output/data/`. The same-order, startup, order-preserving, sharp-asymptotics, and coupled-system drivers additionally write their matching `text_*.tex` files to `output/tables/`. `text_sharp_identity` has data files only; the robustness, shrinking-threshold, and vector-forcing drivers do not export LaTeX tables or figures.

`script_verify_80digits.m` separately checks Tables 3 and 4, the coupled-system reference values, and the vector recurrence/residual identities. Its CSV/MAT stems are `verification_80digits_table3`, `verification_80digits_table4`, `verification_80digits_coupled`, and `verification_80digits_vector_identity`.

For the shrinking-threshold experiment, the manuscript specifies the exponents but not a mesh sequence. The supplied driver fixes `N = 2^7,...,2^20`; this is an explicitly documented package choice rather than an additional setting claimed from the paper.

## Numerical conventions

**Time grids and starts.** `N` is the number of uniform subintervals and `h = (b-a)/N`. The numbered tables use exact additional starting values: one for RBF Euler and classical AB2/PC, and two for the RBF Adams methods and all standalone AM1 comparisons. The classical Euler method starts from the initial value alone. The separate startup experiment uses two classical fourth-order Runge-Kutta steps.

**Error measures.** Tables 1 and 2 use the endpoint error `E_N = |u(b)-v_N|`. Tables 3 and 4 use the maximum error over all mesh nodes. System experiments use the Euclidean norm. Rates are calculated from unrounded errors as `log(E_previous/E_current)/log(N_current/N_previous)`, reducing to the manuscript's base-2 formula on dyadic meshes. The asymptotic diagnostics use the signed error **exact solution minus computed solution**.

**Safeguards.** `RBF methods for IVP/MA_adams.m` implements the `none`, `original`, and `order-preserving` choices. The system extension uses a radial reciprocal, not componentwise division. Activations are counted on the computed histories at the manuscript indices `n = 1,...,N-2`; no clipping, restart, denominator floor, or step rejection is added. An unsafeguarded trajectory that meets an undefined reciprocal is marked incomplete and receives no full-interval error estimate. Corresponding `NaN` entries in robustness data are intentional.

**Implicit and predictor-corrector steps.** Table 1(b) uses Newton iteration with the analytic derivative, a forward Euler guess, and an absolute increment tolerance of `1e-14`. The driver checks that its displayed solves take at most four iterations. The quadratic implicit equation in the order-preserving Example 3 test uses its cancellation-free local root. Coupled-system AM steps solve a 2-by-2 linear system directly. Each PC update uses exactly one corrector evaluation and stores the corrected value in its history.

**Stable asymptotic evaluation.** The scalar and vector forcing diagnostics accumulate analytically derived local residuals rather than subtract nearly equal states on the finest meshes. The scalar sharp-asymptotics driver uses harmonic-number identities in 80-digit arithmetic and checks them against direct recurrences at `N = 512`. The vector identities are cross-checked by `script_verify_80digits.m`.

## Generated files and LaTeX use

All output directories are created automatically under the package root, regardless of the current working directory once a driver has been located. Rerunning a script overwrites its own output files. Previously generated files are never inputs, and the complete `output/` directory can be deleted and regenerated.

CSV files contain column headings and numerical summaries written with `%.17g`. The matching MAT files contain `headers`, `values`, and `metadata`, including the MATLAB version, computer identifier, generation timestamp, and experiment settings. Figure 1 instead saves its stability arrays and figure metadata in `figure1_stability_data.mat`.

The 80-digit computations perform the relevant arithmetic before converting final error summaries to double precision for reporting. Their CSV/MAT summaries are **not archives of all 80 digits** or full high-precision trajectories.

The numbered LaTeX tables use scientific notation with two digits after the decimal point in the mantissa and four decimal places for rates. A missing first-row rate is printed as `--`. To include a generated table, load these packages in the document preamble:

```latex
\usepackage{booktabs}
\usepackage{adjustbox}
```

Then insert, for example:

```latex
\input{output/tables/table1.tex}
```

Each numbered file already contains a complete `table` environment, caption, and label; do not wrap it inside another `table`. Adjust the input path only when the receiving LaTeX document is stored elsewhere. The package generates table source but does not invoke a LaTeX compiler.

## Verification and computational cost

| Command | Scope | Reports in `output/checks/` |
| --- | --- | --- |
| `run_core` | Algebraic consistency checks, all numbered outputs, and basic output-file checks | `core_timings.csv` |
| `verify_package` | The same numbered outputs plus all six double-precision textual drivers | `verification_core_timings.csv`, `verification_core.txt` |
| `run_all` | All numbered outputs, all eight textual drivers, and the additional 80-digit checks | `full_timings.csv`, `verification_full.txt` |

The verification utilities check such properties as algebraic identities, agreement between applicable implementations, expected failure indices, selected numerical reference values, activation counts, and the presence/basic format of generated outputs. `Reproducibility/published_values.m` supplies regression targets only: neither the integrators nor the exporters copy these targets to construct numerical results.

Double-precision comparisons allow explicit roundoff tolerances. Small last-digit differences across MATLAB releases and platforms are therefore possible; a successful check is not a promise of bit-for-bit equality with every printed value. Assertion failures stop the run rather than silently report success. Some previously completed outputs may remain after a failure, so inspect the error message before treating a partial run as complete.

**Execution status.** The manuscript/source correspondence was inspected when preparing this README; the supplied MATLAB scripts were not executed in a native MATLAB session during that documentation preparation. No tested MATLAB release or measured runtime is certified here. Run `verify_package` or `run_all` in the target environment to obtain native checks and timing records.

The 80-digit recurrence checks are expected to cost substantially more than the double-precision table runs. The shrinking-threshold test reaches `N = 1,048,576` and allocates mesh-sized arrays. Actual time and memory use depend on the machine and MATLAB release. After a batch completes successfully, its driver saves the elapsed time for each script. No unmeasured runtime estimates are supplied.

### Coupled-system precision note

The supplied manuscript prints AM and PC maximum errors of `4.5402e-10` and `4.5691e-10` at `N = 1024`. The separate `coupled_80` targets in `published_values.m` are `4.5403e-10` and `4.5690e-10`. The associated displayed rates are unchanged. The double-precision driver and the 80-digit verification intentionally use their respective targets; exported values always come from the calculation.

For the equivariance check, the manuscript reports a discrepancy below `3e-14`. The double-precision script warns when this paper-specific value is exceeded but asserts the more permissive platform tolerance `1e-12`.

An existing comment/metadata field in `script_verify_80digits.m` refers to `NUMERICAL_NOTES.md`. That separate file is absent from the supplied code archive and is not a runtime dependency; the relevant final-digit differences are documented in this section.

## Citation and reuse

When using these experiments or methods, cite the accompanying manuscript by **Roberto Cavoretto and Alessandra De Rossi**, *Order-preserving safeguards and sharp convergence bounds for Matérn RBF time stepping*. The title and author information above follow the supplied manuscript; no publication status, DOI, or release identifier is assumed.

The manuscript identifies the existing MA-RBF repository as the source of the original formula implementations. The supplied code archive contains no license file, and this README does not assign a software license or additional reuse permissions.
