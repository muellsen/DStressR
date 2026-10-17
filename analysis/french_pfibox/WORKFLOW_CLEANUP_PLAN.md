# PFIbox integrated-activity workflow: cleanup roadmap

## Scope and protected result

The manuscript example now uses the signed integral of reconstructed unscaled
PA, not summaries of the deposited per-colony min-max-scaled PA. Keep the
current fit and selected Figure 6 as regression targets during cleanup.
At the user's request, abandoned PFIbox explorations and their generated
outputs have been removed. The exact reconstruction audit is retained.

The four manually chosen illustrations are ompC / Ampicillin-18,
msrB / Sulfamethoxazole-18, recN / Ciprofloxacin-18, and
rmf / NalidixicAcid-18. They were picked from the documented joint-score,
diversity-constrained 10-positive plus 10-negative display selection.

## Source and redistribution

- Study: French, Coutts and Brown (2018), Cell Systems,
  DOI <https://doi.org/10.1016/j.cels.2018.07.004>.
- Dataset: Mendeley Data v1, published 24 August 2018,
  <https://data.mendeley.com/datasets/fmv6vv4bsf/1>;
  DOI <https://doi.org/10.17632/fmv6vv4bsf.1>.
- Dataset licence: CC BY 4.0, verified on the repository page.
- Local source: `data-raw/french_pfibox/source/gfpData.Rdata`.
- Source size: 1,311,202,733 bytes; deliberately Git-ignored.
- SHA-256: `41b3f09bc7a8c3c435af9d03d0e33cede018ba7a5fd0ec900cf6f7d07079ed63`.
- Published cluster assignments: Table S3 / mmc2.xlsx, retained locally with
  separate article provenance. These are published untreated time-course
  cluster labels, not a ground-truth partition of perturbation effects.
- The public PFIbox processing scripts are a later protocol example. Record
  the exact downloaded version/hash; do not describe them as frozen 2018 code.

## Completed cleanup

The retained pipeline has annotated source reconstruction, scalar integration,
model fitting, saved-fit inference, diagnostics and figure assembly stages.
The metadata and scalar-response exports total 6.7 MB compressed and include
source attribution, units, observation keys and checksums. A reusable
`plot_activity_trajectories()` helper supplies both the 20-pair display and
Figure 6b. Abandoned exploratory outputs occupied approximately 1.86 GB.

Validation confirmed exact equality of all 273,792 reconstructed responses.
A refit from the public derived tables reproduces the saved total effects,
p-values and all 5,101 global BH calls. The selected 20 pairs and every
displayed trajectory value are unchanged. The figure-layer contract passes.

Optional additional scalar summaries and further adoption of package plotting
APIs remain possible extensions; they are not part of this retained analysis.

## Readable stages

1. Source README in the data folder: download instructions, expected path,
   version, checksum, source citation, licence and missing-file diagnostics.
2. Import and validate: dimensions, reporter/colony identity, condition codes,
   unique observation IDs, finite-trace eligibility, and reference structure.
3. Reconstruct PA: difference the deposited biomass-normalized fluorescence,
   retain differences 20--216, quadratic LOESS span 0.5, convert to per-hour
   derivative using the five-minute interval. Include the numerical recovery
   check against deposited min-max-scaled PA. Keep signed derivatives.
4. Scalar summaries: trapezoidal integral with configurable T; optionally
   other explicitly named summaries, all using declared windows and units.
   The primary result must reproduce the current integral exactly.
5. Export a stable observation-keyed experimental-design table and long-format
   scalar-response table, with a schema, units, exclusion report and metadata.
   Preserve original source identifiers and distinguish repeated colony
   positions from independent biological replication. Never duplicate controls
   across treatment comparisons.
6. DStressR: annotated prepare_assay(), fit_destress(), then results(); save
   fit separately from inference. Default target here is total effects;
   Figure 6 uses global BH adjustment across 59,301 hypotheses.
7. Figure assembly: computation from standardized results, then rendering.
   Reuse standard package matrix/volcano/p-value plot functions where possible;
   keep the four-pair marking, manuscript typography and composite layout in
   the paper-specific layer. Extract a reusable trajectory plotting function
   with explicit observation grouping, condition, time, units and summaries.

## Export feasibility

The stable exports are `experimental_design.tsv.gz` (3,473,238 bytes) and
`scalar_responses.tsv.gz` (3,250,079 bytes), each with 273,792 rows. Their
README and manifest document the schema, CC BY attribution, transformations
and checksums. The 1.31-GB source and generated analysis caches stay out of Git.

## Regression checks

- 1,797 reporters; 33 treatments; 273,792 finite colony responses.
- 59,301 fitted total effects; 5,101 global BH calls at 0.05.
- Correlation / complete-linkage clustering of 1,796 nonconstant profiles.
- 1,478 reporters matched to published cluster labels; ARI 0.03254586.
- The four selected pairs are identical across heatmap, time courses and volcano.
- Current raw and moderated supplementary results use the same integral fit.
- Both manuscript PDF and HTML must render, with no stale scaled-response
  figures or counts in the PFIbox application and its supplementary notes.

## Literature notes for the four examples

French et al. explicitly mention recN among ciprofloxacin-induced SOS reporters.
Choi et al. (2024), DOI 10.1080/21505594.2024.2367648, report decreased ompC
transcripts after ampicillin in wild-type UPEC. Bollenbach et al. (2009),
DOI 10.1016/j.cell.2009.10.025, report rmf upregulation during trimethoprim
DNA-synthesis stress; this is pathway context, not direct verification of
the nalidixic-acid pair. Grimaud et al. (2001), DOI 10.1074/jbc.M105509200,
establish MsrB's protein-repair function; no exact sulfamethoxazole-direction
reference was established in this pass.
