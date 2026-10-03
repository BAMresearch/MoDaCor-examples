# I22 USAXS processing design

Status: initial implementation validated 2026-09-29; revision backlog agreed
2026-09-30; priorities 0--7 implemented and validated 2026-10-03

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
| 2 | 2 | Generalize angle-to-Q conversion with explicit conventions and uncertainty-aware energy/wavelength input | Low--medium | Reusable for USAXS and step-scanning diffractometers | Explicit angle convention |
| 3 | 3 | Replace `IndexPixels` with a general one-dimensional BaseData-coordinate indexer and make indexed reduction coordinate-agnostic | Medium | Reusable binning outside detector images and Q/Psi conventions | One-dimensional bin/edge/output schema |
| 4 | 8 | Determine one beam centre per physical front/rear scan pair | Medium | Correct physical grouping and cleaner graph lanes | Paired-profile combination rule |
| 5 | 5, 10 | Scale and merge all four series, determine and retain transmission, then subtract the merged background | High | Largest scientific improvement to the current correction chain | Common scale anchor, paired centres, merged-coordinate reduction |
| 6 | 6 | Validate post-subtraction residuals statistically and at detector/gain handoffs | Medium | Distinguishes harmless negative noise from systematic mismatch | Revised merged pipeline |
| 7 | Aside | Allow `XSGeometryFromPixelCoordinates` to accept photon energy as an alternative to wavelength | Low | Reuses the common uncertainty-aware energy/wavelength conversion and simplifies source adaptation | Existing photon conversion helper |

Priority 1 is now implemented: schema-level `for_each` step-block expansion.
The detailed design and implementation record is in the MoDaCor core note
[`pipeline-foreach-expansion.md`](https://github.com/BAMResearch/MoDaCor/blob/main/docs/development/design/pipeline-foreach-expansion.md).
The compact pipeline is 504 rather than 1,255 lines and expands to the same
123 ordinary steps. It reproduces the prior pooled signal, pooled Q, and
transmission scalar exactly for sample scans 978497--978500 against background
scans 977724--977727. The subsequent `AngleToQ` migration reproduced those
values within floating-point roundoff before the independently agreed Q-range
change. The pipeline now uses `q_min = 2e-3 1/nm` for final binning and the
matching final-curve scaling interval.

Priority 3 is now implemented. `IndexByCoordinate` replaces `IndexPixels`, and
`IndexedAverager` consumes only the resulting index map plus its independently
configured value, optional measured axis, and optional mask. At that revision
stage all three real sample acquisitions still produced 386 populated bins
(IDs 0--385), exactly matching the preceding implementation. The later
dark-noise mask deliberately removes one sparsely populated top-Q bin.

Priorities 4 and 5 are now implemented. Each simultaneous front/rear pair is
response-scaled and averaged before centroid determination. The low-gain
rear/front response factor is then refitted near the centred direct beam and
applied to both gains. The resulting low- and high-gain pair curves are fitted,
pooled without a hard handoff, and averaged onto a shared fine signed-Q grid.
Transmission and background subtraction operate on those merged acquisition
curves. The expanded DAG has 184 ordinary steps; mapped blocks keep the
sample/background and gain lanes visible in the authored YAML and graph.

Across the three example samples, the low-gain rear/front response factor is
`29.760`--`29.813`, compared with `29.851` for the reused background. The
sample high-to-low factor is `0.9807`--`1.0066`, while the background factor is
`1.0344`. The merged full-scan transmission values are `0.58102`, `0.58290`,
and `0.58459`. All runs produce 385 populated final bins with first mean Q near
`2.008e-3 1/nm`.

Priority 6 exposed and corrected a missing dynamic-range condition. Positive
rear-diode tails below the measured dark-noise floor had remained eligible for
pair averaging. The preprocessor now retains the dark-rate standard deviation
as a separate operating array. Explicit `CopyDataBundleKeys` and
`DivideDatabundles` steps preserve and divide a copy of the signal to create a
dimensionless `signal_to_dark_noise` diagnostic; `ThresholdMask` then only
compares that ratio with the configured bound. The example uses the original
I22 dark-noise multipliers: 5 for low-gain front, 7 for both rear scans, and no
lower noise cutoff for high-gain front. This leaves roughly 450 rear-diode
points around each direct beam while the front diode supplies the wings. The
upper high-gain-front saturation threshold remains provisional.

Priority 7 is now implemented. `AngleToQ` and
`XSGeometryFromPixelCoordinates` share the same `photon_source`,
`photon_units_source`, and `photon_uncertainties_sources` interface. A shared
`BaseData` helper infers photon energy versus wavelength from Pint
dimensionality and produces wavelength with propagated uncertainties. Because
the photon description is acquisition metadata rather than a dynamically
derived pipeline value, `AngleToQ` reads it directly from its `IoSource`; the
USAXS graph therefore no longer needs eight energy-loading steps.

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

The public conversion step is `AngleToQ`, with an explicit configured
convention:

- `scattering_angle` or `two_theta` uses
  `Q = 4*pi/lambda * sin(angle/2)`;
- `bragg_angle` or `theta` uses
  `Q = 4*pi/lambda * sin(angle)`.

The step retains an optional zero/centre, signed output, units, axes, and
BaseData uncertainty propagation. It loads photon energy or wavelength through
the shared `photon_*` IoSource interface and infers the representation from its
units. The current USAXS yaw is treated as a signed scattering-angle
displacement. No `YawToQ` compatibility alias is retained because the earlier
step had no external consumer. Reusable `BaseData` helpers convert photon
energy and wavelength through `h*c`, preserving units, named uncertainty
components, axes, rank, and weights. The material-attenuation path uses the
same conversion helper.

### 3. General BaseData indexing

The bin-assignment operation does not depend on the names `Q` and `Psi`, nor on
pixels. `IndexByCoordinate` bins exactly one arbitrary `BaseData` coordinate
key. Its configuration defines explicit or generated bin
edges, units, and linear or logarithmic spacing. Non-finite and out-of-range
coordinates receive index `-1`.

Region selection does not belong in the indexer. Static or dynamic coordinate
regions are expressed through ordinary mask-producing steps such as
`ThresholdMask`, and the downstream averager applies the resulting mask. This
keeps selection independently inspectable and avoids embedding periodic or
technique-specific ROI policy in a general bin-assignment module.

The indexer emits a dimensionless index map and the actual bin edges as
`BaseData`. The edges document the assignment boundaries but are not an input
to indexed reduction. Coordinate arrays must be broadcast-compatible with the
indexed data domain. Because no production compatibility is required,
`IndexPixels` is retired rather than retained as a scattering-specific adapter;
the tracked example pipelines are migrated to the generic contract.

Two-coordinate binning is deferred. Useful cases include Q/Psi cake plots and
Qx/Qy reciprocal-space maps, but they require a multidimensional reduction
contract, retained output shape, and multiple physical output axes. Treating
them as a flattened one-dimensional index would hide that materially different
result model.

The agreed separation of concerns is:

1. `IndexByCoordinate` reads one configured coordinate, constructs or accepts
   bin edges, and writes `bin_index` plus diagnostic `bin_edges`. It does not
   inspect signal values or apply masks. Bins are left-inclusive and
   right-exclusive except that the final right edge is included.
2. `IndexedAverager` reads a value, an index map, and an optional mask. It never
   reads bin edges and never reassigns points. An optional `axis_key` identifies
   a second value, such as Q, that is averaged over exactly the same accepted
   observations and weights to provide the measured output axis.
3. Only bins with positive total weight are emitted, ordered by their original
   integer bin ID. The original IDs are retained as a diagnostic. Empty bins
   are omitted because they have neither a measured signal nor an actual mean
   coordinate.

This keeps the index map reusable by a future indexed sum or integrator without
adding reduction modes to `IndexedAverager`. Two-coordinate binning remains a
separate future result model.

### 4. Conservative lower Q limit

The agreed lower limit for the final I22 USAXS example is
`q_min = 2e-3 1/nm`. This replaces the provisional `8e-4 1/nm` value in both
final binning and any final-curve scaling interval whose lower boundary is
intended to match the usable Q range. Full signed scans remain available for
centering, scaling studies, transmission, and wing diagnostics.

All three example acquisitions execute the current 184-step pipeline with this
limit. Their final curves contain 385 points, with the first mean bin Q at about
`2.008e-3 1/nm`.

### 5. Transmission from the merged four-series measurement

The implemented transmission calculation uses all four normalized series for
an acquisition: front and rear diodes at low and high gain. The four sample
series are brought onto one response scale and merged; the four background
series are treated equivalently. Transmission is then calculated over the full
valid merged scans:

```text
T = integral(merged_sample) / integral(merged_background)
```

Merging means an uncertainty-aware average after response scaling, not a raw
sum that would count the same photons two or four times. The low-gain front
diode is the response anchor within both sample and background. Saturated or
otherwise invalid points contribute no weight. There is no fitted scale between
the sample and background acquisitions, because their relative magnitude
contains the unknown transmission.

The sample merged curve is divided by `T`; only then is the merged background
remapped by nearest neighbour and subtracted. The reusable remapper preserves
the sample Q values and invalidates points beyond the background domain.

Response calibration follows two levels:

1. Fit the low-gain rear diode onto the simultaneous low-gain front diode over
   `-8e-4 <= Q <= 8e-4 1/nm`, where the direct beam provides strong overlap.
   Apply that same rear/front factor to the rear diode at both gains. This is
   provisional pending beamline confirmation that both rear-amplifier settings
   have the same gain ratio, but the measured rear low/high direct-beam
   integrals support that assumption to about one percent.
2. Average the aligned front and scaled-rear values separately for low and high
   gain. Fit the merged high-gain curve onto the merged low-gain reference over
   the negative-wing interval `-5e-3 <= Q <= -2e-3 1/nm`.

The low/high curves are not joined at a threshold. Both valid curves are
concatenated, assigned to a shared fine signed-Q grid, and averaged within each
occupied bin. This makes the overlap statistically useful and avoids an
arbitrary handoff discontinuity.

Both fit levels retain the fitted scalar, its named uncertainty component,
point count, actual fit range, and reduced chi-square. Reduced chi-square is
`1.93`--`3.80` for the sample diode fits and `5.97`--`12.20` for the sample gain
fits. These values warn that the selected `subread_sem` does not fully explain
the pointwise curve mismatch; the formal factor uncertainty is therefore not
silently enlarged or presented as a goodness-of-fit substitute.

The named factor uncertainties propagate through BaseData arithmetic and are
preserved when curves with different uncertainty components are concatenated.
They are scalar calibration uncertainties, however, and hence correlated
across every point to which a factor was applied. The present diagonal
uncertainty representation carries their marginal contribution but indexed
averaging and integration cannot preserve that covariance. Final uncertainties
must not yet be interpreted as a complete covariance-aware budget.

### 6. Negative background-subtracted intensities

Negative intensities are retained. They are acceptable when their distribution
is consistent with propagated noise. Validation should examine standardized
residuals, sign balance in background-dominated regions, contiguous runs versus
Q, and detector/gain handoff regions. A persistent offset, long same-sign run,
or Q-dependent trend is evidence of systematic background, transmission, or
scale mismatch; isolated values of order one propagated standard uncertainty
are not a correction failure.

All three packaged samples currently have zero negative final bins after the
dark-noise mask correction. That observation is descriptive, not evidence that
the high-Q signal should equal zero: the corrected sample curves remain
strongly positive and may contain real scattering. The notebook therefore
does not calculate a generic final-signal z-score against zero. It instead
plots standardized low/high differences over their actual scaling interval,
where agreement is expected, using the explicitly selected `subread_sem`
component. It also reports negative-bin counts and longest contiguous negative
runs should later datasets contain them.

The standardized overlap diagnostics and reduced chi-square values show no
hard handoff—none exists—but do show more scatter than `subread_sem` alone
predicts. Over the gain-fit interval the standardized sample differences have
means between `-0.04` and `-0.13`, but standard deviations of `2.45`--`3.50`;
the reused background has mean `-0.22` and standard deviation `1.82`. Thus the
fitted levels are not systematically displaced, while the selected pointwise
uncertainty component is too narrow. This remains a calibration/uncertainty-
model diagnostic rather than a reason to switch curves at a threshold.

### 7 and 9. Generic repeated-step expansion and loading

Modules should keep their current specific task rather than acquiring lists of
unrelated targets or pairs. Repetition belongs at the pipeline-composition
layer. Investigation of the runtime and the then-current 123-step USAXS graph showed that a
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

Each pair yields one shared beam centre used by both readouts. Low/high and
sample/background profiles must not be combined merely because they use the
same centering module. The `prepare_acquisition_center` repeated block invokes
the chain once for the sample acquisition and once for the background, while
keeping independent low- and high-gain centroid branches inside each item:

1. copy the rear signal so centre preparation cannot mutate a scientific
   readout;
2. fit the low-gain rear copy to the low-gain front response in yaw space with
   the configured `subread_sem` weighting;
3. apply the same fitted rear/front scalar to the low- and high-gain rear
   copies;
4. concatenate each front and scaled-rear profile while recording both diode
   identity and source-local point position;
5. use the source-local position directly as the `IndexedAverager` group map,
   averaging the actual yaw and signal with the selected uncertainty weights;
6. determine an independent iterative centroid for the low- and high-gain
   merged profiles; and
7. copy each shared centre to the two original readouts from that gain scan
   before `AngleToQ`.

Front and rear yaw arrays are identical within each physical scan. The
source-position index is therefore an exact aligned-observation grouping, not
a coordinate binning approximation. `alignment_key: yaw` verifies this
pointwise at runtime before concatenation. Each diode's dynamic mask is applied
before fitting and pooling. In particular, saturated high-gain front points
are absent and the rear diode supplies the direct-beam region.

The early full-scan factor remains isolated on temporary centre-estimation
bundles. After centering and Q conversion, the scientific rear/front factor is
refitted over the configured direct-beam Q interval for response scaling and
transmission.

### 10. Normalization factors are outputs

A monitor-normalization step must retain the exact applied normalization factor
as a named BaseData output, with units, propagated uncertainty, weights/mask
semantics, and provenance. A pointwise I0 normalization factor is an array; the
merged transmission is a scalar dimensionless factor. They should not be
conflated merely because both are divisors.

For transmission, the pipeline exposes:

- merged sample integral;
- merged background integral;
- scalar `transmission_factor` with uncertainty;
- the transmission-normalized sample curve.

Factor determination and application are separate graph steps. The numerator
is copied to `transmission_factor`, divided by `background_integral`, and then
used to normalize `sample_merged`. The original `sample_integral` and
`background_integral` remain addressable diagnostics rather than being
overwritten.

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
- exposure-adjusted dark-rate standard deviation as a separate masking
  operating array, not as an automatically propagated signal uncertainty;
- analyser yaw;
- I0 mean and SEM;
- scan count time;
- measured photon energy;
- scan identifiers, titles, source files, and preprocessing version.

For scan exposure `t_scan` and matching dark exposure `t_dark`, the stored
signal is

```text
adjusted = scan - mean(dark / t_dark) * t_scan
dark_noise_std = std(dark / t_dark) * t_scan
```

`dark_noise_std` supplies the divisor for the retained dimensionless
`signal_to_dark_noise` diagnostic. A separate `ThresholdMask` compares that
ratio with the dynamic-range threshold; the mask step performs no
normalization. The dark-noise array is not attached as a propagated
uncertainty component. The pipeline later divides valid signal points by
`t_scan`. No I0 dark subtraction is applied:
the acquisitions do not provide a clean I0-dark measurement. Non-positive I0
samples are dynamically masked before division as invalid readouts; this is a
validity check, not a dark correction.

## Beam-centre determination

The beam centre is an intensity centroid rather than a Gaussian-fit centre.
This avoids assuming a symmetric or structureless analyser rocking curve.

For each physical low- or high-gain scan:

1. Create axis-independent validity masks for non-finite data, invalid I0,
   diode saturation, and configured signal-to-dark-noise limits.
2. Normalize to count time and I0.
3. Scale and uncertainty-average the aligned front/rear pair as described
   above.
4. Use the merged profile to locate a coarse maximum.
5. Select the contiguous valid peak region within a configurable physical
   half-width around that maximum.
6. Subtract an optional scalar baseline. Zero is the default because the
   electronic dark has already been removed.
7. Clip remaining negative intensity weights to zero and compute

   `centre = sum(yaw * intensity) / sum(intensity)`.
8. Recenter the window and repeat until convergence or a small iteration
   limit is reached.
9. Store the centroid, propagated uncertainty, contributing-point count,
   window, and convergence diagnostics.

The centre found from the merged pair is applied to both original diodes from
that physical scan. Low- and high-gain scans are centred independently.

## Signed momentum transfer and wing policy

The `AngleToQ` module treats analyser yaw as a signed scattering angle, uses
the measured photon energy through the shared uncertainty-aware
energy-to-wavelength helper, and evaluates

```text
Q = 4*pi/wavelength * sin((yaw - yaw_zero)/2)
```

It retains signed Q so the two analyser wings remain distinguishable. The
positive wing may contain additional analyser-crystal scattering. It remains
part of the full-scan transmission integral because those photons were not
absorbed by the sample. The final scattering product uses the negative wing
and negates Q before logarithmic binning. Both wings are retained before this
selection for diagnostics.

## Transmission and background subtraction

Transmission is evaluated after count-time and I0 normalization and before
background subtraction. Both the sample and empty-furnace inputs are the
merged front/rear and low/high acquisition curves described above. Their Q
samples are stably sorted for quadrature.
The scalar transmission is

```text
T = integral(sample) / integral(empty_furnace)
```

The merged sample is divided by `T`, then the merged empty-furnace curve is
remapped and subtracted.

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

## Scaling, concatenation, and averaging

Response scaling occurs before transmission normalization and background
subtraction. `FindScaleFactor1D` uses a lognormal fit based on the
uncertainty-weighted mean log ratio over positive, valid overlap points. The
user explicitly selects the propagated uncertainty component used for fitting
weights; the pipeline does not invent a combined uncertainty for this purpose.
Before interpolation, exact repeated coordinates are consolidated by
inverse-variance averaging using that same selected uncertainty component.

The scaled curves are pooled with a generic `ConcatenateDatabundles` step.
It concatenates matching `BaseData` entries, converts compatible units, and
preserves named uncertainties, weights, masks, and provenance. Its `sort_by`
option defaults to `None`. When set to a concatenated BaseData key, all
concatenated entries are reordered together; sorting is ascending unless
`descending: true` is configured. Sorting is optional because
`IndexByCoordinate` does not require monotonic input. The opt-in
`uncertainty_key_policy: fill_zero` takes the union of named uncertainty
components and supplies zero for a component absent from one input. This is
used when only a scaled branch carries the corresponding fit uncertainty.

The pooled data are reduced using the generic indexing modules:

1. `IndexByCoordinate` assigns bins from actual Q values and retains the edges
   as a diagnostic. Selection masks remain independent of bin assignment.
2. `IndexedAverager` computes the weighted mean signal and the weighted mean of
   the actual Q values that entered each bin. Nominal bin centres are not used
   as output coordinates.
3. With `stats_keys: [signal, Q]`, it reports signal and Q SEM/STD from the
   actual within-bin values.
4. Additional diagnostics report raw bin count, positive-weight count, summed
   weight, and effective sample count.

The user selects `uncertainty_weight_key` when inverse-variance weighting is
enabled. No uncertainty components are combined implicitly.

## Current pipeline order

```text
reduce subreads and subtract exposure-adjusted diode darks
  -> load four sample and four background readouts
  -> build finite/I0/saturation/SNR masks
  -> normalize to count time and I0
  -> provisionally scale and uncertainty-average aligned front/rear pairs on copies
  -> determine four merged-pair centroids and apply each to its two readouts
  -> convert yaw as a scattering angle to signed Q for all eight readouts
  -> fit the low-gain rear/front response near the direct beam
  -> apply each acquisition's rear/front factor to both gains
  -> uncertainty-average the aligned front/rear points at each gain
  -> fit each merged high-gain curve onto its low-gain reference
  -> concatenate low/high points and average on a common fine signed-Q grid
  -> integrate full merged scans and determine transmission
  -> divide the merged sample by transmission
  -> nearest-neighbour remap and subtract the merged background
  -> select negative wing and use abs(Q)
  -> IndexByCoordinate
  -> IndexedAverager
  -> write final I(Q), actual mean Q, Q scatter, masks, and diagnostics
```

`Integrate1D` uses its opt-in coordinate-sorting mode for the merged acquisition
curves; strict monotonicity remains the module default.

## Core MoDaCor work

- Add a reusable intensity-centroid model and `FindCenterOfMass1D` step.
- Add a signed scattering-angle-to-Q kernel and process step.
- Extend `FindScaleFactor1D` with lognormal fitting and an explicit uncertainty
  weighting key.
- Extend `Integrate1D` with opt-in coordinate sorting and duplicate-coordinate
  consolidation.
- Add `SubtractInterpolated1D` with nearest and linear remapping.
- Add `ConcatenateDatabundles` with optional coordinated sorting.
- Add `IndexByCoordinate` for one-coordinate bin assignment and retire
  `IndexPixels`.
- Add bin-count, summed-weight, positive-weight-count, and effective-count
  diagnostics to `IndexedAverager`.
- Let `AngleToQ` and `XSGeometryFromPixelCoordinates` accept photon energy or
  wavelength through one source interface and the shared uncertainty-aware
  `BaseData` conversion helper.

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
- Merged-pair centroids are stable against reasonable changes in centroid
  window and SNR threshold.
- The assumption that the rear-diode amplifier changes by the same factor as
  the front path is confirmed by the beamline scientist or replaced by an
  explicit calibration.
- Full-scan merged transmission remains stable under reviewed diode-response,
  gain-overlap, and fine-grid perturbations.
- Nearest and linear background remapping are tested on mismatched, descending,
  and partially overlapping coordinates without silent extrapolation.
- Scale fits use only positive valid overlap data and remain stable under a
  reviewed fit-range perturbation.
- Final Q values equal the weighted means of actual contributing Q values, and
  reported Q STD/SEM and contribution diagnostics match direct calculations.
- Front/rear and low/high overlap residuals show no unexplained systematic
  discontinuity; there is deliberately no hard handoff threshold.

Resolution desmearing, sample-thickness normalization, and absolute-intensity
calibration are deliberately later phases.
