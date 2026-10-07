# Diamond Light Source I22 example

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23207188.svg)](https://doi.org/10.5281/zenodo.23207188)

This example demonstrates processing of operando SAXS, WAXS, and USAXS data
from Diamond Light Source beamline I22. Four SAXS/WAXS sample measurements are
paired with one empty-cell background and detector-specific calibration and
mask files. Three USAXS sample acquisitions share one empty-furnace background.

## Packaged inputs

- Sample masters: `i22-978003.nxs`, `i22-978008.nxs`, `i22-978013.nxs`, and
  `i22-978018.nxs`.
- Empty-cell background master: `i22-977723.nxs`.
- Glassy-carbon absolute-intensity calibration master: `i22-966700.nxs`, with
  its detector and readout sidecars. Its native thickness field is missing, so
  preprocessing records the independently known thickness of 1.0 mm.
- Absolute reference curve: `Glassy Carbon L average.dat` (Q in Å⁻¹ and
  intensity in cm⁻¹ sr⁻¹). The accompanying email was reviewed but is omitted
  because it contained no useful provenance information.
- Unobstructed-beam transmission reference: `i22-977723.nxs` for this data set.
  It is configured separately from the background in each notebook so another
  reference can be used, for example when the background is an empty capillary.
- Detector, I0, beamstop-diode, and user-tetramm HDF5 sidecars for each master.
- SAXS/WAXS calibration and mask files in `data/processing/`.
- USAXS empty-furnace scans `977724`--`977727` and sample scan groups
  `978497`--`978500`, `978502`--`978505`, and `978507`--`978510`, with their
  NeXus/HDF5 sidecars in `data/usaxs/`.
- Recommended ordinary-processing pipeline:
  `pipelines/I22_SAXS_WAXS_solids_operando.yaml`. It processes only SAXS and
  WAXS, in one graph, and uses the SAXS glassy-carbon result for both detector
  branches. The detector-specific pipelines remain available to the chunked
  transport examples during their migration: `pipelines/I22_SAXS_solids_operando.yaml`
  plus one WAXS nosecone profile,
  `pipelines/I22_WAXS_solids_operando_usaxs_saxs_waxs.yaml` for the combined
  USAXS/SAXS/WAXS configuration without an aluminium attenuator, or
  `pipelines/I22_WAXS_solids_operando_standard_saxs_waxs.yaml` for the standard
  SAXS/WAXS configuration with the 1 mm aluminium attenuator.
- Test-only DAWN comparison pipelines:
  `pipelines/I22_SAXS_DAWN_crosscheck.yaml` and
  `pipelines/I22_WAXS_DAWN_crosscheck.yaml`.

The small NeXus master files contain relative external links to their HDF5
sidecars. Keep each master and its sidecars together in `data/`.

The data payload is licensed under CC BY 4.0; see `DATA_LICENSE.json` and
`DATA_LICENSE.txt`.

## Data availability

The notebooks use version 1.0.0 of the Diamond Light Source I22 example
dataset, archived on Zenodo under the version-specific DOI
[10.5281/zenodo.23207188](https://doi.org/10.5281/zenodo.23207188). The tracked
`data-manifest.json` binds this example to the exact archive and records its
size, SHA-256 checksum, and the checksums of all extracted files.

The preferred citation is:

> Pauw, Brian Richard, Müller-Elmau, Johanna, Bleßmann, Finn Eric, & Smith,
> Andrew. (2026). *Diamond Light Source I22 SAXS/WAXS/USAXS example data for
> MoDaCor* (Version 1.0.0) [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.23207188

## Running the example

If `data/` is absent, first run `python tools/data_repository.py download
DLS/I22` from the repository root. This downloads and verifies dataset version
1.0.0 from Zenodo.

Start Jupyter from the repository root or this directory, select the prepared
MoDaCor kernel, and choose one focused notebook:

- `I22_USAXS_server.ipynb` preprocesses the four-scan USAXS acquisitions and
  runs the signed-Q correction, transmission, background-remapping, diode
  scaling, and indexed averaging pipeline. Its compact `step_blocks` document
  expands to 184 independently traceable execution steps. See
  `USAXS_PROCESSING_DESIGN.md` for the processing contract and provisional
  settings.
- `I22_solids_server.ipynb` performs ordinary non-chunked batch processing.
- `I22_solids_chunked_buffer.ipynb` uploads notebook-sliced `BufferSource`
  chunks.
- `I22_solids_chunked_hdf.ipynb` lets the runtime slice `HDFSource` inputs.
- `I22_solids_chunked_tiled.ipynb` runs the same direct-slice workload through
  a notebook-owned Tiled service.

The old `I22_solids_server_operando_preprocessed.ipynb` filename is retained as
a short index so existing links lead to these notebooks. Shared discovery,
preprocessing, registrations, and explicit I22 plan construction live in
`i22_helpers.py`; the local Tiled lifecycle lives in `i22_tiled.py`.
The packaged measurements used the special USAXS/SAXS/WAXS nosecone, so the
notebooks set `BEAMLINE_CONFIGURATION = "usaxs_saxs_waxs"`. Change it to
`"standard_saxs_waxs"` only for measurements made with the standard nosecone.

The processing notebooks write compact MoDaCor-facing files below
`work/preprocessed/`. Ordinary results go below `work/output/`; chunked results
go below `work/chunk_server/`. Generated preprocessing links are relative, so
the work tree remains movable with this example. The Tiled notebook requires
the `tiled-tests` MoDaCor extra.

Preprocessing only reshapes or summarizes incompatible frame-wise metadata.
It averages channel 1 of both `bsdiodes` and `I0`, including their standard
errors on the mean. The unobstructed-beam measurement supplies the scalar
calibration ratio `mean(bsdiodes) / mean(I0)`. Frame-wise sample transmission
is `(bsdiodes / I0) / calibration_ratio`; its propagated SEM includes both
sample readouts and both reference readouts. The results are written to
`/entry1/sample/transmission` and `/entry1/sample/transmission_sem` in each
preprocessed file. MoDaCor resolves detector geometry and corrections from the
original NeXus metadata and the packaged calibration files.

The recommended combined SAXS/WAXS pipeline is authored with four compact
`step_blocks`: frame normalization, detector geometry, detector correction,
and final sample integration. Their `for_each` items expand into independent
SAXS, WAXS, background, and glassy-carbon lanes as appropriate; the runtime
still executes the same 96 ordinary steps. Each detector is first normalized
by its own count time, followed by its own transmission and incident-`I0`
steps. This ordering keeps the non-integrating transmission and flux readouts
downstream of the integrating count-time correction. SAXS and WAXS remain separate,
visually continuous lanes through Poisson uncertainty generation, transmission,
incident-flux normalization, frame averaging, silicon-detector efficiency,
polarization, and specimen-thickness normalization. They use the same metadata
and correction parameters without repeatedly merging and splitting the graph.
Geometry, masks, background subtraction, and azimuthal integration remain
detector-specific. After applying the shared absolute scale, the two corrected
curves are concatenated and Q-sorted into one `saxs_waxs_1d` dataset and one
live 1D plot; `detector_index` records SAXS as 0 and WAXS as 1. The sample
masters do not contain a
usable thickness value, so preprocessing records 50 µm from their acquisition
titles.

Measurement 966700 follows the complete SAXS sample path. Its averaged curve
is divided by its independently known 1.0 mm thickness before a lognormal,
uncertainty-weighted `FindScaleFactor1D` fit to the supplied reference over
`0.12–2.0 nm⁻¹`. This conservative range avoids relying on the high-Q tail of
the desmeared APS USAXS reference. I0 is represented as the incident count rate
in `count/s`; dividing the count-time-normalized detector signal by I0 therefore
leaves a dimensionless detector response. The fitted factor is dimensionless
and is applied unchanged to both thickness-normalized detector curves. The 50 µm
thickness is represented as `50e-6 m`, so the resulting intensity units are
`1/(m sr)`.
The WAXS glassy-carbon data are therefore not fitted independently. The fixed
calibration branch is cached and reused when only the sample source changes.

## Current validation status

The SAXS/WAXS results below were originally validated with MoDaCor 1.8.0.
Their pipeline files have since been migrated to the current development
baseline documented in the repository root. All packaged I22 pipelines prepare
with that baseline. The shared helper discovers exactly four packaged samples,
validates matching `(1679, 1475)` calibration/mask shapes, and preprocesses the
samples, background, and glassy-carbon measurement with working relative HDF5
links. The USAXS tests validate
the compact 184-step graph and preprocessing, and all four packaged USAXS
acquisitions preprocess successfully. Complete SAXS and WAXS pipeline
equivalence and the full server processing loop remain follow-up workflow
validation tasks.

On 2026-10-07, the published Zenodo record and public file endpoint were
checked against the frozen local archive. The filename, 2,283,982,640-byte
size, and Zenodo MD5 checksum matched exactly; all 113 local payload files also
passed their manifest SHA-256 checks.

An October 2026 real-data smoke run of the combined graph fitted the
dimensionless factor `1.68143098e-17` over `0.12–2.0 nm⁻¹` (88 reference points; reduced
χ² 0.176). Both the SAXS and WAXS sample outputs used that exact same factor.
Fit-window checks from 0.5 through 3.3 nm⁻¹ changed the
factor by about 1.3%, with convergence by 2.0 nm⁻¹. An exploratory independent
WAXS fit had a pronounced Q-dependent residual and reduced χ² 7.25, supporting
the decision not to use the desmeared reference at WAXS Q for scale estimation.
A second measurement completed as a partial rerun and reproduced the cached
shared factor exactly. I0 and the beamstop detector are incident and transmitted
count rates in `count/s`; their ratio-derived transmission remains dimensionless.
A follow-up per-step unit trace with the thickness expressed as `50e-6 m`
verified `1/sr` after solid-angle correction and `1/(m sr)` after thickness
normalization and after applying the shared dimensionless factor for both
detectors. Q remains stored as `1/m`; the
live plots numerically convert it to `1/nm`.

The BufferSource notebook is configured for all four
measurements and ten ten-frame chunks per measurement. It derives independent
SAXS and WAXS schemas from pilot chunks and stores two detector-specific run
groups in one physical HDF5 file. Each 40-chunk plan is initialized, populated,
inspected, and finalized independently. Lightweight trace events are retained
under `/processing/tracer/<run_name>/chunks/<chunk_id>/`. The complete
80-pipeline-run exercise remains an interactive beamline validation. A reduced
real-data smoke run has completed two ten-frame chunks for each detector in one
shared file and verified both pilot comparisons and trace manifests.

Two further notebooks run that workload with the server reading sample
slices directly. The first uses `HDFSource`; the second exposes the same
preprocessed files through a notebook-owned read-only Tiled server and uses
`TiledSource`. Typed source bindings apply one frame selection to the detector,
I0 mean and SEM, sample transmission and SEM, and detector count time. The finalized
HDFSource output is compared with the BufferSource result, and the TiledSource
output is compared with the HDFSource result, one stored chunk at a time. A
one-chunk real SAXS smoke run has completed through both direct transports,
including pipeline execution, trace publication, finalization, and exact
cross-transport comparison; the complete 80-run variants remain interactive.

The changing sample is chunked, but the background is deliberately handled as
one reusable aggregate per detector. Each pilot loads, corrects, and reduces
the complete background; later sample chunks reuse the reduced branch. The
complete detector stack and its cached/working copies, masks, uncertainties,
and numerical temporaries must therefore fit comfortably in worker memory.
SAXS and WAXS are processed sequentially and their sessions are deleted between
runs so the two background working sets do not coexist. A background that does
not fit requires a separate mergeable weighted-reduction workflow or an
explicit frame-pairing policy; averaging chunk means blindly is not sufficient.

## Supplementary poster visualizations and testing

The supplementary performance material under `supplementary/poster_2026/` is
retained to reproduce poster graphics and exploratory timing measurements. It
is not part of the main I22 correction example or its normal validation path.
`supplementary/poster_2026/I22_performance_benchmark.ipynb` compares
fresh-session server processing, sample-only partial reruns, and direct-slice
chunked-HDF throughput; its shared implementation is in the adjacent
`i22_performance.py`.

The benchmark uses the numerical form of the recommended pipelines: all
corrections and azimuthal integration are retained, while interactive plots and
intermediate 2D file sinks are removed. Preprocessing, server startup, session
construction, and the chunk-schema pilot are outside the reported timings. On
2026-09-15, an Apple M1 Max with 64 GiB RAM running Python 3.14.5 and MoDaCor
1.8.0 gave median end-to-end times of 29.2 s (full) and 9.24 s (sample rerun)
for SAXS, and 25.2 s and 6.14 s for WAXS. Forty direct-HDF ten-frame chunks per
detector gave median throughputs of 8.59 SAXS frames/s and 8.97 WAXS frames/s,
including chunk publication.

The notebook writes raw observations, per-step durations, machine metadata, and
PNG/SVG/PDF plots below
`work/supplementary/poster_2026/performance/`. This regenerable output tree is
ignored by Git. The processing date and MoDaCor version are also recorded in
the notebook so later benchmark runs can be distinguished from these results.

## Provisional processing values and follow-up checks

- The old absolute-intensity scalar `3.8e-15` is retained only for the DAWN
  cross-check and is still carried by older chunk-source schemas for
  compatibility. The recommended operational pipeline does not use it.
- The DAWN cross-check pipelines reproduce selected recorded behaviour and are
  not the recommended physical-correction pipelines.
- The 50 µm specimen thickness is title-derived and must be confirmed against
  the sample log before a release-quality absolute calibration is claimed.
- The combined pipeline intentionally transfers the SAXS-derived scale to WAXS;
  the WAXS glassy-carbon/reference shape mismatch remains a diagnostic rather
  than a second calibration.
- A general pipeline-level output-unit conversion for intensity and Q remains
  to be added; `UnitsLabelUpdate` is not a numerical unit converter.
- A future thickness path should infer the effective thickness from measured
  X-ray absorption together with atomic composition and gravimetric density,
  rather than relying on the acquisition-title value.
- The NeXus metadata include experiment identifiers `sm43108-1` and
  `sm43533-1`, plus facility account usernames `ckn54496` and `hck38156`.
  No email address was found; the apparent email-like value was ordinary `@`
  temperature notation in sample titles. These non-secret identifiers are
  retained as acquisition provenance in the public data package.
