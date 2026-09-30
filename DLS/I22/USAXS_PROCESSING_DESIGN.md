# I22 USAXS processing design

Status: initial implementation validated 2026-09-29; revision backlog agreed
2026-09-30; priority-1 schema expansion implemented and validated

## Purpose

This document defines the first MoDaCor processing workflow for the I22
Bonse--Hart USAXS measurements. It replaces the earlier exploratory workflow
with a small instrument preprocessor and an explicit pipeline built from
reusable MoDaCor modules.

The packaged example uses one empty-furnace acquisition and three sample
acquisitions:

| Acquisition | Low-gain dark | Low-gain scan | High-gain dark | High-gain scan |
| --- | ---: | ---: | ---: | ---: |
| Empty furnace | 977724 | 977725 | 977726 | 977727 |
| Sample at 1153.2 degrees C | 978497 | 978498 | 978499 | 978500 |
| Sample at 1162.7 degrees C | 978502 | 978503 | 978504 | 978505 |
| Sample at 1175.4 degrees C | 978507 | 978508 | 978509 | 978510 |

Each acquisition contains low- and high-gain rocking scans. Each rocking scan
records a front HG-USAXS diode, a rear BS diode, an I0 monitor, analyser yaw,
and count time. The diode electronics integrate over the exposure time.

## Revision decisions and implementation order (2026-09-30)

The first end-to-end result is a starting point rather than the final pipeline.
The following decisions supersede the corresponding parts of the initial
implementation where they differ.

| Priority | Concerns | Work item | Cost | Benefit | Principal dependencies |
| ---: | --- | --- | --- | --- | --- |
| 0 | 1, 4, 6 | Correct terminology, set the agreed final `q_min`, and define how negative results are assessed | Very low | Immediate clarity and a safer example default | None |
| 1 | 7, 9 | Add schema-level `for_each` step-block expansion, including loading blocks | Medium | Very high reduction in authored YAML without weakening module contracts | None |
| 2 | 2 | Generalize `YawToQ` into an angle-convention-aware `AngleToQ` | Low--medium | Reusable for USAXS and step-scanning diffractometers | Explicit angle convention |
| 3 | 3 | Add a general BaseData-coordinate indexer and retain `IndexPixels` as a scattering convenience interface | Medium--high | Reusable binning outside detector images and Q/Psi conventions | Generic bin/ROI/periodicity schema |
| 4 | 8 | Determine one beam centre per physical front/rear scan pair | Medium | Correct physical grouping and cleaner graph lanes | Paired-profile combination rule |
| 5 | 5, 10 | Scale and merge all four series, determine and retain transmission, then subtract the merged background | High | Largest scientific improvement to the current correction chain | Common scale anchor, paired centres, merged-coordinate reduction |
| 6 | 6 | Validate post-subtraction residuals statistically and at detector/gain handoffs | Medium | Distinguishes harmless negative noise from systematic mismatch | Revised merged pipeline |

Priority 1 is now implemented: schema-level `for_each` step-block expansion.
The detailed design and implementation record is in the MoDaCor core note
[`pipeline-foreach-expansion.md`](https://github.com/BAMResearch/MoDaCor/blob/main/docs/development/design/pipeline-foreach-expansion.md).
The compact pipeline is 496 rather than 1,255 lines and expands to the same
123 ordinary steps. It reproduces the prior pooled signal, pooled Q, and
transmission scalar exactly for sample scans 978497--978500 against background
scans 977724--977727. The priority-0 decisions are recorded here; the `q_min`
configuration change should be made alongside the next scientific pipeline
revision so that the tracked pipeline and its validation results change
together.

### 1. Eight diode readouts, not twelve scientific inputs

There are four dark-corrected diode readouts for the sample and four for the
background: eight primary readouts in total. The four I0 arrays correspond to
the sample-low, sample-high, background-low, and background-high physical
scans. They are auxiliary monitors, not additional USAXS readouts.

The initial implementation temporarily represents those four monitor arrays as
helper DataBundles because the existing pairwise arithmetic modules consume
DataBundles. That creates twelve workspace bundles during loading, but it does
not change the eight-readout scientific input model. A mapped loader or monitor
normalizer may instead attach or consume I0 as auxiliary BaseData without
presenting it as four additional measurement lanes.

### 2. General angle-to-Q conversion

The public conversion step should become `AngleToQ`. Renaming alone is not
sufficient because the configured angle must have an explicit convention:

- `scattering_angle` or `two_theta` uses
  `Q = 4*pi/lambda * sin(angle/2)`;
- `bragg_angle` or `theta` uses
  `Q = 4*pi/lambda * sin(angle)`.

The step retains an optional zero/centre, measured photon energy or wavelength,
signed output, units, axes, and BaseData uncertainty propagation. The current
USAXS yaw is treated as a signed scattering-angle displacement. The numerical
kernel is already phrased in terms of scattering angle; the main work is the
public configuration contract, migration compatibility, tests, and docs.

### 3. General BaseData indexing

The bin-assignment operation should not fundamentally depend on the names `Q`
and `Psi`, nor on pixels. A new generic indexer should accept one or more
coordinate specifications referring to arbitrary BaseData keys. Each
specification defines explicit or generated bin edges, units, linear/log
spacing, and whether it is a binned coordinate or only a region-of-interest
selector. Periodic coordinates require an explicit period and wrap policy.

For multiple binned coordinates, the indexer should combine the per-coordinate
indices into one flat integer index while retaining bin shape and edge metadata.
All coordinate arrays must be broadcast-compatible with the indexed data
domain. `IndexPixels` can remain as a backwards-compatible scattering adapter
that translates Q/Psi configuration into this generic contract.

The numerical binning is straightforward. Most of the cost is defining units,
periodic ranges, multidimensional flattening, metadata, and exact dependencies.
Generalizing `IndexedAverager` is a separate decision; the generic index map can
initially continue to feed the existing averager for signal and Q reduction.

### 4. Conservative lower Q limit

The agreed lower limit for the final I22 USAXS example is
`q_min = 2e-3 1/nm`. This replaces the provisional `8e-4 1/nm` value in both
final binning and any final-curve scaling interval whose lower boundary is
intended to match the usable Q range. Full signed scans remain available for
centering, scaling studies, transmission, and wing diagnostics.

### 5. Transmission from the merged four-series measurement

The target transmission calculation uses all four normalized series for an
acquisition: front and rear diodes at low and high gain. The four sample series
are brought onto one response scale and merged; the four background series are
treated equivalently. Transmission is then calculated over the full valid
merged scans:

```text
T = integral(merged_sample) / integral(merged_background)
```

Merging means an uncertainty-aware average after response scaling, not a raw
sum that would count the same photons two or four times. Sample and background
must share a common response anchor, for example low-gain front fixed to one,
so independent arbitrary rescaling cannot bias `T`. Saturated or otherwise
invalid points contribute no weight.

The sample merged curve is divided by `T`; only then is the merged background
remapped and subtracted. Background subtraction already follows transmission
normalization in the initial pipeline. The required change is therefore the
four-series merged input to transmission and subtraction, rather than their
relative ordering.

This revision moves diode/gain response scaling before transmission and
background subtraction. It consequently requires a reviewed rule for deriving
response factors without relying on the already background-subtracted curve.

### 6. Negative background-subtracted intensities

Negative intensities are retained. They are acceptable when their distribution
is consistent with propagated noise. Validation should examine standardized
residuals, sign balance in background-dominated regions, contiguous runs versus
Q, and detector/gain handoff regions. A persistent offset, long same-sign run,
or Q-dependent trend is evidence of systematic background, transmission, or
scale mismatch; isolated values of order one propagated standard uncertainty
are not a correction failure.

### 7 and 9. Generic repeated-step expansion and loading

Modules should keep their current specific task rather than acquiring lists of
unrelated targets or pairs. Repetition belongs at the pipeline-composition
layer. Investigation of the runtime and the 123-step USAXS graph showed that a
single-step map is too narrow: it reduces repeated module declarations but does
not represent repeated chains such as “for sample and background, do x, y, and
z.” The preferred feature is therefore schema-level `for_each` expansion of a
block containing one or more ordinary step templates. A one-step block covers
the original mapped-step case.

An illustrative one-step block is:

```yaml
step_blocks:
  load_usaxs_readouts:
    for_each:
      PD_SLF:
        output_key: SLF
        location: sample::/entry1/low_gain_front/signal
      PD_SLR:
        output_key: SLR
        location: sample::/entry1/low_gain_rear/signal
    steps:
      load:
        module: AppendProcessingData
        configuration:
          processing_key: "${output_key}"
          databundle_output_key: signal
          signal_location: "${location}"
```

Each local step is expanded once per item into an ordinary ProcessStep instance
before graph validation. It therefore retains its own step id, configuration
validation, dependency contract, tracing, failure reporting, and partial-rerun
behavior. Local prerequisites are resolved within each item instance. The
implementation must instantiate one child module per item; it must not
repeatedly mutate and reuse one stateful ProcessStep instance.

This mechanism applies equally to loaders, masks, conversions, arithmetic, and
other arbitrary registered modules. It also handles multi-step repeated blocks
for sample/background processing. Source-scoped blocks make the 28 current
loading nodes concise without creating an I22-specific bulk loader. The
runtime DAG remains flat; graph export carries block/item/local-step metadata
so a renderer can group children into physical-scan lanes or show a collapsed
summary.

### 8. Beam centre belongs to one physical scan pair

There are four physical front/rear pairs:

- sample low gain: `SLF + SLR`;
- sample high gain: `SHF + SHR`;
- background low gain: `BLF + BLR`;
- background high gain: `BHF + BHR`.

Each pair must yield one shared beam centre used by both readouts. Low/high and
sample/background profiles must not be combined merely because they use the
same centering module. A repeated step block can invoke the selected paired
estimator four times, but the estimator's scientific combination rule belongs
to the module, not to pipeline expansion.

When this item is implemented, the paired estimator must specify whether it
centres a response-scaled combined profile or combines two independently
estimated centres. Saturation masks and uncertainties must participate in that
choice. The current rear-derived centre copied to the front diode remains the
fallback until the paired rule is selected.

### 10. Normalization factors are outputs

A monitor-normalization step must retain the exact applied normalization factor
as a named BaseData output, with units, propagated uncertainty, weights/mask
semantics, and provenance. A pointwise I0 normalization factor is an array; the
merged transmission is a scalar dimensionless factor. They should not be
conflated merely because both are divisors.

For transmission, the revised pipeline should expose at least:

- merged sample integral;
- merged background integral;
- scalar `transmission_factor` with uncertainty;
- the transmission-normalized sample curve.

These values should use explicit diagnostic names rather than overwriting an
object still named `sample_integral`. Whether factor determination and factor
application remain separate steps or are exposed by one normalizer will be
settled with the normalization-module contract; in either case the factor must
remain addressable ProcessingData.

## Package boundary

The I22 preprocessor, acquisition discovery, example pipeline, and notebook
belong in `MoDaCor_examples/DLS/I22`. Format-independent mathematics and
generally useful process steps belong in the main MoDaCor package.

The preprocessor is deliberately lean. It may decode the facility layout,
reduce electronic subreads, perform the exposure-aware dark subtraction, and
write provenance. Scientific corrections remain visible pipeline operations.

## Preprocessed data contract

One compact NeXus file is written for each four-scan acquisition. It exposes
four ordinary one-dimensional diode readouts:

- low-gain front diode;
- low-gain rear diode;
- high-gain front diode;
- high-gain rear diode.

Loading one sample and one background therefore gives the pipeline eight
dark-adjusted readouts.

For every readout the file stores:

- integrated, exposure-adjusted dark-subtracted signal;
- subread SEM and dark-offset SEM as separate uncertainty components;
- analyser yaw;
- I0 mean and SEM;
- scan count time;
- measured photon energy;
- scan identifiers, titles, source files, and preprocessing version.

For scan exposure `t_scan` and matching dark exposure `t_dark`, the stored
signal is

```text
adjusted = scan - mean(dark / t_dark) * t_scan
```

The pipeline later divides by `t_scan`. No I0 dark subtraction is applied:
the acquisitions do not provide a clean I0-dark measurement. Non-positive I0
samples are dynamically masked before division as invalid readouts; this is a
validity check, not a dark correction.

## Beam-centre determination in the initial implementation

The beam centre is an intensity centroid rather than a Gaussian-fit centre.
This avoids assuming a symmetric or structureless analyser rocking curve.

For each physical low- or high-gain scan:

1. Create axis-independent validity masks for non-finite data, invalid I0,
   diode saturation, and configured signal-to-uncertainty limits.
2. Normalize to count time and I0.
3. Use the unsaturated rear diode to locate a coarse maximum.
4. Select the contiguous valid peak region within a configurable physical
   half-width around that maximum.
5. Subtract an optional scalar baseline. Zero is the default because the
   electronic dark has already been removed.
6. Clip remaining negative intensity weights to zero and compute

   `centre = sum(yaw * intensity) / sum(intensity)`.
7. Recenter the window and repeat until convergence or a small iteration
   limit is reached.
8. Store the centroid, propagated uncertainty, contributing-point count,
   window, and convergence diagnostics.

The centre found from the rear diode is applied to both diodes from that
physical scan. Low- and high-gain scans are centred independently.

## Signed momentum transfer and wing policy in the initial implementation

The yaw-to-Q module uses measured photon energy and

```text
Q = 4*pi/wavelength * sin((yaw - yaw_zero)/2)
```

It retains signed Q so the two analyser wings remain distinguishable. The
positive wing may contain additional analyser-crystal scattering. It remains
part of the full-scan transmission integral because those photons were not
absorbed by the sample. The initial final scattering product uses the negative
wing and converts it to positive `abs(Q)` before logarithmic binning. Positive,
negative, and asymmetry diagnostics remain available.

## Transmission and background subtraction in the initial implementation

Transmission is evaluated after count-time and I0 normalization and before
background subtraction. The primary estimate uses the unsaturated low-gain
front diode; the low-gain rear diode is a consistency check.

Both sample and empty-furnace curves are integrated over the full valid scan.
The yaw samples are stably sorted for quadrature. Exact repeated encoder
coordinates are consolidated by arithmetic-mean intensity, with independent
uncertainties propagated to the mean; this avoids assigning an arbitrary
finite integration interval to either repeated observation.
The scalar transmission is

```text
T = integral(sample) / integral(empty_furnace)
```

The sample readouts are divided by `T`, then each matching empty-furnace
readout is remapped and subtracted.

The reusable subtraction module supports two explicit remapping modes:

- `nearest`: choose the nearest background coordinate for every sample
  coordinate. This is the simple, setting-free mode formerly described as
  "no interpolation".
- `linear`: linearly interpolate between bracketing background coordinates.

Neither mode extrapolates. Points outside the background domain are
invalidated. Nearest mode intentionally has no distance-tolerance setting;
linear mode can optionally invalidate points whose bracketing interval exceeds
a configured maximum width. Linear interpolation propagates each uncertainty
component using its interpolation coefficients, not by interpolating standard
deviations directly. The sample coordinate is preserved as the output
coordinate.

## Scaling, concatenation, and averaging in the initial implementation

After transmission normalization and background subtraction, the four sample
readouts are scaled to one selected reference readout. `FindScaleFactor1D`
gains a lognormal mode based on the uncertainty-weighted mean log ratio over
positive, valid overlap points. The user explicitly selects the propagated
uncertainty component used for fitting weights; the pipeline does not invent a
combined uncertainty for this purpose. Before interpolation, exact repeated
coordinates are consolidated by inverse-variance averaging using that same
selected uncertainty component.

The scaled curves are pooled with a generic `ConcatenateDatabundles` step.
It concatenates matching `BaseData` entries, converts compatible units, and
preserves named uncertainties, weights, masks, and provenance. Its
`sort_by` option defaults to `None`. When set to a concatenated BaseData key,
all concatenated entries are reordered together; sorting is ascending unless
`descending: true` is configured. Sorting is optional because `IndexPixels`
does not require monotonic input.

The pooled data are reduced using the existing scattering modules:

1. `IndexPixels` assigns bins from actual Q values. For pure one-dimensional Q
   binning, `Psi` becomes optional when no azimuthal ROI is requested.
2. `IndexedAverager` computes the weighted mean signal and the weighted mean of
   the actual Q values that entered each bin. Nominal bin centres are not used
   as output coordinates.
3. With `stats_keys: [signal, Q]`, it reports signal and Q SEM/STD from the
   actual within-bin values.
4. Additional diagnostics report raw bin count, positive-weight count, summed
   weight, and effective sample count.

The user selects `uncertainty_weight_key` when inverse-variance weighting is
enabled. No uncertainty components are combined implicitly.

## Initial pipeline order

```text
reduce subreads and subtract exposure-adjusted diode darks
  -> load four sample and four background readouts
  -> build finite/I0/saturation/SNR masks
  -> normalize to count time and I0
  -> determine four rear-diode centroids
  -> convert yaw to signed Q for all eight readouts
  -> integrate full scans and determine transmission
  -> divide sample readouts by transmission
  -> nearest or linear remapped background subtraction per readout
  -> select negative wing and use abs(Q)
  -> determine and apply lognormal readout scales
  -> ConcatenateDatabundles
  -> IndexPixels
  -> IndexedAverager
  -> write final I(Q), actual mean Q, Q scatter, masks, and diagnostics
```

`Integrate1D` gains an opt-in coordinate-sorting mode for jittered acquisition
axes. Its existing strict-monotonic behavior remains the default.

## Core MoDaCor work

- Add a reusable intensity-centroid model and `FindCenterOfMass1D` step.
- Add a signed scattering-angle-to-Q kernel and process step.
- Extend `FindScaleFactor1D` with lognormal fitting and an explicit uncertainty
  weighting key.
- Extend `Integrate1D` with opt-in coordinate sorting and duplicate-coordinate
  consolidation.
- Add `SubtractInterpolated1D` with nearest and linear remapping.
- Add `ConcatenateDatabundles` with optional coordinated sorting.
- Make `Psi` optional in `IndexPixels` for Q-only one-dimensional binning.
- Add bin-count, summed-weight, positive-weight-count, and effective-count
  diagnostics to `IndexedAverager`.

Every public step requires focused model/module tests, exact dependency-contract
assertions, export through `modacor.modules`, and regenerated reference docs.

## Example work

- Add the lean I22 USAXS preprocessor and scan-set discovery helper.
- Add focused synthetic preprocessing tests and validate all four packaged
  acquisitions end to end.
- Add the tracked USAXS pipeline YAML.
- Copy the structure of `I22_solids_server.ipynb` into a focused USAXS server
  notebook that preprocesses the four acquisition sets, previews the graph,
  runs the three samples against the reusable background, and displays peak,
  wing, transmission, scaling, and final-I(Q) diagnostics.
- Document the measured 14 keV energy and any provisional electronic
  saturation thresholds next to their YAML configuration.

## Validation gates

- Exactly four dark-adjusted curves are produced per acquisition and eight are
  loaded for a sample/background run.
- Dark subtraction is equivalent to subtracting dark count rate after time
  normalization, including uncertainty propagation.
- Rear-diode centroids are stable against reasonable changes in centroid
  window and SNR threshold.
- Full-scan transmission estimates from the two unsaturated low-gain diodes
  agree within their reviewed uncertainty budget.
- Nearest and linear background remapping are tested on mismatched, descending,
  and partially overlapping coordinates without silent extrapolation.
- Scale fits use only positive valid overlap data and remain stable under a
  reviewed fit-range perturbation.
- Final Q values equal the weighted means of actual contributing Q values, and
  reported Q STD/SEM and contribution diagnostics match direct calculations.
- The final curve has no unexplained steps at diode or gain handoffs.

Resolution desmearing, sample-thickness normalization, and absolute-intensity
calibration are deliberately later phases.
