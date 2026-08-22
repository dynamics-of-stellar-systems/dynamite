# Appendix draft: Memory-efficient orbit-weight solving for large-magnitude NNLS problems

*Draft appendix text describing the 2026-08 fused-X / streamed-reads work.
Numbers are measured on the production omega Centauri library
(see dev_notes/weight_solve_rss_profile.md). Written to be adapted into
LaTeX; a ready-made table skeleton is included at the end.*

---

## A.x Memory-efficient construction and solution of the weight-solving problem

### A.x.1 The memory wall

The weight-solving step of our models casts the search for non-negative
orbit weights as a non-negative least-squares (NNLS) problem whose constraint
matrix $\mathbf{A}$ has one row per observable constraint and one column per
orbit bundle. For the omega Centauri models presented in this Paper — four
simultaneous kinematic data sets comprising two line-of-sight velocity-
distribution sets and two proper-motion sets, with $N_{\rm orb} = 45\,000$
orbit bundles after tube-orbit mirroring — the matrix has
$371\,212 \times 45\,000$ elements and occupies 124 GiB in double precision.

Solving this problem with the augmented-Lagrangian (ALM) scheme of
Appendix~\ref{app:alm} requires a second, equally large working array: the
augmented matrix

$$
\mathbf{X} \;=\;
\begin{pmatrix}
\sqrt{\mu}\,\mathbf{1}' \\
\mathbf{A}_{\rm rest}
\end{pmatrix},
$$

which replaces the total-mass row of $\mathbf{A}$ (whose error normalisation,
$10^{-8}$, raises the condition number of the naive problem to
$\sim 5\times 10^{22}$) by a finite penalty row, rescaled to unit Euclidean
norm column by column. Until recently both matrices were built explicitly
and held resident for the duration of the solve, together with all orbit
libraries' kinematic histograms ($\sim$165 GiB). Profiled end to end, a
single weight solve occupied 472 GiB of main memory in steady state and
peaked at 626 GiB (Fig.~\ref{fig:rss}); on our 1416 GB compute node this
restricted execution to one weight solve at a time, and because the full
parameter-space exploration of this Paper requires several hundred solves,
the campaign was computationally infeasible in this configuration.

We addressed this with three changes to the weight solver: fused assembly of
the augmented matrix, streamed reading of the orbit library, and single-
precision evaluation. All three are individually gated by bit-level
regression tests against the previous implementation, and none changes the
mathematical statement of the optimisation problem.

### A.x.2 Fused assembly

In the previous implementation, $\mathbf{A}$ was assembled block by block
(one block per kinematic set, plus intrinsic and projected mass rows),
error-normalised, and only then copied into the augmented buffer
$\mathbf{X}$ and column-scaled. We now allocate $\mathbf{X}$ once and write
every block directly into it, applying each row's error division as the
block is written; the penalty row replaces the total-mass row, which
survives only as a small vector $\mathbf{a}_0 = \mathbf{1}/\epsilon_0$. The
column norms are computed afterwards from the finished array with the same
operation the two-step path used, so $\mathbf{X}$, its column norms and the
target vector are *bitwise identical* to the previous construction; we have
verified this with SHA-256 digests of the full production matrix (182
contiguous row slabs) between the two implementations.

Removing $\mathbf{A}$ from the solve requires replacing its two remaining
uses. Both follow from the exact relations

$$
\mathbf{A}_{{\rm rest},j} = c_j\,\tilde{\mathbf{X}}_{{\rm rest},j},
\qquad
\|\mathbf{A}_{\cdot j}\|^2 = a_{0j}^2 + c_j^2 - \mu,
$$

where $c_j$ is column $j$'s norm, $\tilde{\mathbf{X}}$ the column-scaled
body of $\mathbf{X}$, and $\mathbf{r}$ the plain residual
$\mathbf{A}\mathbf{w}-\mathbf{b}$. The residual itself is available without
evaluation, since the solver's internal residual satisfies
$r_i = b_i - (\mathbf{A}\mathbf{w})_i$ exactly for all rows $i\geq1$, and the
total-mass row contributes a scalar dot product with $\mathbf{a}_0$. The
final $\chi^2$ vector and the Karush-Kuhn-Tucker (KKT) optimality
certificate are therefore computed from these identities instead of from
matvecs against $\mathbf{A}$, which additionally removes three full passes
over the matrix (the certificate previously recomputed the residual, the
gradient $\mathbf{A}^{\rm T}\mathbf{r}$, and the column norms).

### A.x.3 Streamed library reads

Assembly no longer requires the orbit libraries' histograms to be resident
simultaneously. The weight solver reads one kinematic set at a time —
decompressing the archive, combining and mirroring the orbit families,
transforming to observables, writing the block into $\mathbf{X}$, and
freeing the set before the next is read. Because the histogram archives
interleave the apertures of several kinematic sets within one file, the bulk
binary parsers were extended with a mask that keeps their sequential record
walk synchronised with the file layout while scattering values only for the
set being processed; passing anything other than the complete per-file
aperture enumeration to these parsers silently desynchronises the walk (we
caught and fixed exactly this failure during development, and pin it with a
dedicated regression test). On glibc systems, freed multi-GiB allocations
are returned to the operating system immediately, which we verified
directly, making the streaming frees reliable.

The trade-off is repeated decompression: at production scale, streamed
assembly takes 1888 s against 249 s when the library is read once, an
overhead of $\sim$27 min per model on our network-attached storage. Since
the solve itself is bounded below by hours, this is negligible per model
but worth avoiding when many models share one orbit library; the mode is
therefore optional (`stream_orblib_reads`).

### A.x.4 Single precision and validation

Because $\mathbf{A}$ no longer exists, all large arrays on the solved path
can be held in single precision (`nnls_dtype: float32`), halving every
figure quoted above. The solver's coordinate descent accumulates in double
precision internally, and the returned weights agree with the double-
precision solution's $\chi^2_{\rm kin}$ to $1.6\times10^{-5}$ relative on
omega Centauri, consistent with earlier validations of single-precision
weight solving on NGC 6278.

Our acceptance criteria distinguish three levels of agreement. (i)
*Bitwise identity* is required of constructions that claim to rearrange,
not redefine, arithmetic: verified between the fused and two-step paths,
between streamed and non-streamed assembly, and between single-run repeats
of the same binary, using synthetic libraries and digest comparison of the
production matrix. (ii) Round-off-level agreement ($\leq 10^{-11}$
relative in $\chi^2$) is achieved between old and new code paths within one
code version. (iii) Comparison against solutions recorded under earlier
versions reaches only $\sim 10^{-6}$ relative in $\chi^2$: the threaded
coordinate-descent solver does not reduce floating-point quantities in a
fixed order, so successive runs diverge at the last significant digit and
the ALM iterate selection amplifies this to a redistribution of weight at
the $\sim 10^{-6}$ level among essentially unchanged sets of orbits (6537
of ~6540 positive weights shared; maximum absolute difference
$2.8\times10^{-6}$). This is a property of the parallel solver, not of the
memory work, and sets the reproducibility floor for any threaded
implementation of this algorithm.

### A.x.5 Results

Table~\ref{tab:rss} summarises peak resident memory per weight solve, measured
by sampling `/proc` at 1 s resolution through instrumented end-to-end runs
on the real production library, together with the number of concurrent
solves this permits on our 192-core, 1416 GB node. The combination of all
three changes reduces the per-solve footprint by a factor of 3.3 and raises
the achievable concurrency from one solve to six, with the double-precision
options remaining available for verification runs.

*[LaTeX table skeleton]*

```latex
\begin{table}
  \centering
  \caption{Peak resident memory per weight solve ...}
  \label{tab:rss}
  \begin{tabular}{lccc}
    \hline
    Configuration & Peak memory & Concurrent solves & Validation \\
    \hline
    Two-step construction, float64 & 626 GiB & 1 & reference \\
    Fused assembly                 & 502 GiB & 2 & $\chi^2$ equal to $1.3\times10^{-11}$ \\
    \quad + streamed reads         & 380 GiB & 3 & bitwise-equal construction \\
    Single precision               & 251 GiB & 4 & $\chi^2_{\rm kin}$ equal to $1.6\times10^{-5}$ \\
    \quad + streamed reads         & 190 GiB & 6 & bitwise-equal weights \\
    \hline
  \end{tabular}
\end{table}
```

Reproducibility: all changes, regression tests, profiling scripts and
measurement logs are part of our public DYNAMITE fork (branch
`orblib-performance`); the measurement methodology, raw timelines and the
failure-mode post-mortem are documented in the repository alongside the
code.
