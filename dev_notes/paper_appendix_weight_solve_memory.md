# Appendix draft: Memory-efficient orbit-weight solving for large-magnitude NNLS problems

*Draft appendix text describing the 2026-08 fused-X / streamed-reads work.
All numbers are measured on the production omega Centauri library and
instrumented end-to-end runs (dev_notes/weight_solve_rss_profile.md,
rss_fused_sprint_log.md). Written to be adapted into LaTeX; a ready-made
table skeleton is included at the end.*

---

## A.x Memory-efficient orbit-weight solving

### A.x.1 Motivation

The dynamical models in this Paper use orbit superposition: a library of
orbits, integrated in a trial potential, is assigned non-negative weights
such that the weighted sum of the orbits' observable properties reproduces
the data. Finding those weights is a non-negative least-squares (NNLS)
problem,

$$
\min_{\mathbf{w}\,\geq\,\mathbf{0}} \;\; \tfrac{1}{2}
\left\| \mathbf{A}\mathbf{w} - \mathbf{b} \right\|^2 ,
$$

where each column of $\mathbf{A}$ holds one orbit's contribution to every
constraint and $\mathbf{b}$ stacks the observed constraints, both normalised
by their measurement errors. For the omega Centauri models presented here,
four kinematic data sets are fitted simultaneously — two line-of-sight
velocity-distribution sets and two proper-motion data sets — together with
the intrinsic and projected mass constraints. With $N_{\rm orb}=45\,000$
orbit bundles (after tube-orbit mirroring) and $371\,212$ constraint rows,
$\mathbf{A}$ alone occupies 124 GiB in double precision. Because the models
of this Paper require several hundred such solves across the explored
parameter space, the weight solver's memory footprint, not its floating-
point throughput, set the pace of the entire campaign: as we describe below,
the unmodified pipeline peaked at 626 GiB per solve on our 1416 GB compute
node, allowing exactly one solve at a time.

This appendix documents three changes to the solving pipeline that reduce
the per-solve peak by a factor of 3.3 (to 190 GiB) and raise the number of
concurrent solves on that node from one to six. None of them changes the
optimisation problem being solved: every construction that claims to
rearrange arithmetic rather than redefine it is gated by bit-level
regression tests against the previous implementation.

### A.x.2 Anatomy of the baseline footprint

We instrumented an end-to-end weight solve with phase markers and a 1 Hz
resident-memory sampler to attribute the footprint precisely. Five
components contribute:

| component | size | resident during |
| --- | --- | --- |
| kinematic histograms (all sets, retained) | ~129 GiB | library read → end of solve |
| histogram read transients (raw tube+box families) | +35 GiB | library read only (peak 164 GiB) |
| constraint matrix $\mathbf{A}$ | 124 GiB | assembly → end of solve |
| augmented matrix $\mathbf{X}$ (Sect. A.x.3) | 124 GiB | X-build → end of solve |
| coordinate-descent working arrays | +125–250 GiB | each ALM multiplier iteration |
| **total** | **626 GiB peak** (472 GiB median in the solve loop) | |

Two structural facts stand out. First, the two large matrices coexist only
because $\mathbf{X}$ is *derived from* $\mathbf{A}$: the classic pipeline
assembles $\mathbf{A}$ block by block, then copies it into the augmented
array and column-scales it there. Second, the retained histograms are needed
only transiently — their information enters $\mathbf{A}$ once, at assembly —
yet they are held for the entire solve because the read step precedes
assembly and has no reason to release them. These are the two inefficiencies
exploited below. The third change, single precision, simply halves whatever
remains.

### A.x.3 The augmented formulation and fused assembly

Following Appendix~\ref{app:alm}, the total-mass row of $\mathbf{A}$ — whose
error normalisation $\epsilon_0=10^{-8}$ raises the condition number to
${\sim}5\times10^{22}$ and defeats the inner solver's convergence test — is
replaced by an augmented-Lagrangian penalty row. Writing
$\mathbf{A}_{\rm rest}$ for $\mathbf{A}$ without row 0, the inner problems
are bound-constrained least-squares solves in the scaled variables

$$
\mathbf{X} =
\begin{pmatrix} \sqrt{\mu}\,\mathbf{1}' \\ \mathbf{A}_{\rm rest} \end{pmatrix},
\qquad
\tilde{\mathbf{X}} = \mathbf{X}\,\mathrm{diag}(c^{-1}),
\qquad
c_j = \|\mathbf{X}_{\cdot j}\|_2,
$$

$$
\min_{\boldsymbol{\beta}\geq\mathbf{0}} \tfrac{1}{2}
\left\|\tilde{\mathbf{X}}\boldsymbol{\beta} - \mathbf{y}(\lambda)\right\|^2,
\qquad
\mathbf{w} = \mathrm{diag}(c)\,\boldsymbol{\beta},
\qquad
y_0(\lambda) = \sqrt{\mu}\,(1+\lambda/\mu),
\qquad
\lambda \leftarrow \lambda - \mu\,(\mathbf{1}^{\rm T}\mathbf{w}-1).
$$

Unit-$L_2$ column scaling is an exact change of variables: positive diagonal
rescaling preserves the non-negativity constraint, and the objective values
are preserved since $\tilde{\mathbf{X}}\boldsymbol{\beta}
=\mathbf{X}\,\mathbf{w}$ evaluated at $\boldsymbol{\beta}=c\odot\mathbf{w}$.
Because the multiplier update changes only the scalar $y_0$, the matrix is
built once and warm-started subproblems converge in one to two inner
iterations.

The memory work lies in how $\mathbf{X}$ comes into being. The previous
pipeline materialised $\mathbf{A}$ first and copied; we now allocate
$\mathbf{X}$ once and stream every constraint block directly into it,
applying each row's error division as the block is written. The penalty row
replaces the total-mass row, which survives only as the small vector
$\mathbf{a}_0=\mathbf{1}/\epsilon_0$. Column norms are computed afterwards
from the finished array by the same operation the two-step path used —
deliberately not accumulated incrementally during assembly, which would
reorder floating-point additions — so $\mathbf{X}$, $c$ and $\mathbf{y}$ are
*bitwise identical* to the classic construction. We verify this by SHA-256
digests of all 182 contiguous row slabs of the production matrix between
implementations.

Removing $\mathbf{A}$ then requires replacing its remaining uses downstream,
which follow from two exact identities. First, the plain residual needs no
evaluation for rows $i\geq1$, since the target vector satisfies
$y_i=b_i$: the solver's internal residual already holds
$r_i=b_i-(\mathbf{A}\mathbf{w})_i$. Row 0 contributes one scalar dot
product against $\mathbf{a}_0$. The full $\chi^2$ vector is therefore

$$
\chi^2_i \;=\;
\begin{cases}
\big(\mathbf{a}_0^{\rm T}\mathbf{w}-b_0\big)^2, & i=0,\\[4pt]
r_i^2, & i\geq1 .
\end{cases}
$$

Second, the KKT certificate requires the gradient
$\mathbf{g}=\mathbf{A}^{\rm T}\mathbf{r}$ and the column norms of
$\mathbf{A}$. Both are recovered without $\mathbf{A}$:

$$
g_j \;=\; a_{0j}\,r_0 \;+\; c_j\,
\big(\tilde{\mathbf{X}}_{\rm rest}^{\rm T}\mathbf{r}_{\rm rest}\big)_j,
\qquad
\|\mathbf{A}_{\cdot j}\|^2 \;=\; a_{0j}^2 \;+\; c_j^2 \;-\; \mu ,
$$

the first by linearity and the scaling relation
$\mathbf{A}_{{\rm rest},j}=c_j\,\tilde{\mathbf{X}}_{{\rm rest},j}$, the
second because $c_j^2=\mu+\|\mathbf{A}_{{\rm rest},j}\|^2$ includes the
penalty row. The scaled violation
$\max_j g_j / (\|\mathbf{A}_{\cdot j}\|\,\|\mathbf{r}\|)$ retains its
Cauchy-Schwarz bound in $[0,1]$. Together these remove three full passes
over $\mathbf{A}$ per solve and, more importantly, allow $\mathbf{A}$ never
to exist: **peak memory falls from 626 to 502 GiB**, exactly the 124 GiB of
the eliminated matrix.

### A.x.4 Streamed library reads

Assembly does not need the kinematic histograms of all four data sets to be
resident simultaneously, yet the baseline holds ${\sim}129$ GiB of them
through assembly and the hours-long solve. We restructured the read so that
one kinematic set at a time is decompressed, combined and mirrored,
transformed to observables, written into its block of $\mathbf{X}$, and
freed before the next is read. On glibc systems the freed multi-GiB pages
are returned to the operating system immediately — verified directly with
an allocation probe up to 128 GiB — so the streaming frees are reliable.

One subtlety proved instructive. The bulk binary parsers walk each archive's
sequential Fortran record structure, and their record count is derived from
the aperture list they are handed: handing them a *subset* of a file's
apertures silently desynchronises the walk and corrupts every orbit after
the first (this failure occurred, was diagnosed from its symptom — a
$\chi^2$ three orders of magnitude too large with no exception raised — and
is now pinned by a regression test). The parsers therefore accept the
complete per-file aperture enumeration plus a boolean mask of which
apertures to scatter; the walk stays synchronised while only the requested
set allocates memory. Streamed assembly is bitwise identical to the
non-streamed construction on the real library.

The trade-off is repeated decompression: streamed assembly takes 1888 s
against 249 s when the library is read once ($\sim$27 min per model on our
network storage), negligible beside the multi-hour solve but avoidable when
many models share one library, so the mode is optional.
**Peak memory falls from 502 to 380 GiB.**

### A.x.5 Single precision

With $\mathbf{A}$ gone, every large array on the solved path can be held in
single precision, halving the remaining footprint: **380 → 190 GiB**
($\mathbf{X}\approx62$ GiB plus per-iteration working arrays, which shrink
correspondingly). The inner solver accumulates reductions in double
precision; the converged kinematic chi-square agrees with the double-
precision path to $1.6\times10^{-5}$ relative on omega Centauri, consistent
with earlier single-precision validations on NGC 6278.

### A.x.6 Validation methodology and reproducibility floor

Our acceptance criteria distinguish three levels of agreement, and it is
worth stating them because they generalise to any performance refactoring of
numerical infrastructure:

1. *Bitwise identity* is required of constructions claiming to rearrange,
   not redefine, arithmetic: verified between fused and two-step assembly,
   between streamed and non-streamed reads, and between repeat runs of the
   same binary — using synthetic libraries with interleaved multi-set files
   and SHA-256 slab digests of the full production matrices.
2. *Round-off-level agreement* ($\leq10^{-11}$ relative in $\chi^2$) holds
   between old and new code paths within a single code version.
3. Comparison against results recorded under *earlier versions* reaches
   only ${\sim}10^{-6}$ relative in $\chi^2$. The threaded coordinate-
   descent solver does not fix the order of its floating-point reductions,
   so successive runs diverge in the last significant digit and the ALM
   iterate selection redistributes weight at the ${\sim}10^{-6}$ level
   among essentially unchanged orbit sets: 6537 of ${\sim}6540$ positive
   weights are shared between runs, with maximum absolute difference
   $2.8\times10^{-6}$ and identical support to within $0.2\%$. This
   reproducibility floor is a property of any threaded implementation of
   this algorithm, not of the memory work described here.

### A.x.7 Results

Table~\ref{tab:rss} summarises the measured peaks and what each change
buys. The combined effect is a factor-3.3 reduction in per-solve footprint
and a rise from one to six concurrent solves on our 192-core, 1416 GB node;
double-precision modes remain available for verification runs.

*[LaTeX table skeleton]*

```latex
\begin{table}
  \centering
  \caption{Peak resident memory per weight solve, measured end-to-end on
    the omega Centauri production library ($371\,212\times45\,000$
    float64), and the resulting number of concurrent solves on a
    1416\,GB node at 85\% safe utilisation.}
  \label{tab:rss}
  \begin{tabular}{llcc}
    \hline
    Step & Configuration & Peak memory & Concurrent solves \\
    \hline
    --    & Two-step construction, float64      & 626 GiB & 1 \\
    (1)   & Fused assembly                      & 502 GiB & 2 \\
    (2)   & \quad + streamed library reads      & 380 GiB & 3 \\
    (3)   & Single precision                    & 251 GiB & 4 \\
          & \quad + streamed library reads      & 190 GiB & 6 \\
    \hline
  \end{tabular}
\end{table}
```

Reproducibility: all changes, regression tests, profiling scripts and raw
measurement timelines are part of our public DYNAMITE fork (branch
`orblib-performance`); the profiling methodology, memory-anatomy
decomposition and the streaming-desynchronisation post-mortem are documented
alongside the code.
