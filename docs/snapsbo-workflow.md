# SnapSBO workflow for TrialEmulationj 1.1.1

This guide maps the SnapSBO non-operative-management analysis to the guided
landmark workflow. It contains variable names and analysis rules only; no
patient-level data are distributed with TrialEmulationj.

## Recommended analysis sequence

1. Open the analysis-ready SnapSBO OMV file in jamovi.
2. Complete **1. Target-Trial Protocol** and prespecify the 24-hour primary and
   48-hour sensitivity windows.
3. Use **2. Causal DAG** to confirm the baseline adjustment set.
4. Open **3. Landmark Trial (24/48 h)** and assign the variables below.
5. Review the timing audit and cohort flow before interpreting any estimate.
6. Review propensity overlap, effective sample size, and the Love plot.
7. Compare AIPW with Hajek IPTW and repeat the prespecified missingness and
   timing-policy sensitivity analyses.

## Landmark variable mapping

| Role | Analysis-ready variable | Rule |
|---|---|---|
| Treatment A / stratum | `ngt` | NGT present versus absent |
| Treatment B / exposure | `wsec` | WSEC challenge present versus absent |
| Index time | `admit_datetime` | Admission is time zero |
| Treatment B time | `wsec_datetime` | Date-time of WSEC administration |
| Failure component 1 | `rescue_surgery` | Surgery after failed NOM |
| Failure component 1 time | `surgery_decision_datetime` | Decision for surgery |
| Failure component 2 | `death` | In-hospital death |
| Failure component 2 time | `death_datetime` | Date-time of death |
| Discharge/censoring time | `discharge_datetime` | Date-time or date of discharge |
| Cluster | `center_analysis` | Missing centres assigned to one `Unknown` level |
| Person ID | `record_id` | Unique patient record |

Choose composite outcome construction, date-time mode, hours, 24 and 48 hours,
and clustered inference. The default v1.1.1 discharge policy excludes a row at a
landmark when continued observation cannot be verified because discharge time
is missing, unparseable, or negative.

## Baseline covariates

The final set should follow the prespecified DAG. Candidate baseline variables
used in the prior analysis include age, sex, comorbidities, prior abdominal
surgery, previous WSEC, raw RCRI, and counts or prespecified categories of prior
operative and non-operative SBO episodes.

Do not use `RCRI_adj`, because subtracting a point according to subsequent
surgery makes the covariate outcome-informed. Avoid adjusting for mediators,
colliders, treatment consequences, or post-landmark information.

## Coding and missingness

Blank values may be recoded as absence only when the data-collection process
supports that interpretation. The historical SnapSBO code treated blank NGT,
WSEC, death, comorbidity, previous WSEC, and prior SBO counts as absence or zero.
TrialEmulationj records every requested recode in its audit table.

Missing cluster identifiers should be resolved before analysis or explicitly
assigned to one shared `Unknown` cluster. Missing baseline covariates should be
handled using the prespecified complete-case or multiple-imputation analysis;
neither choice repairs unmeasured confounding.

## Known data limitations

The source dataset does not contain an NGT-placement time suitable for the
landmark analysis. Leaving Treatment A time blank therefore assumes that the
recorded NGT status describes stratum membership at each landmark. This matches
the historical subgroup analysis but cannot verify 24/48-hour NGT timing.

The timing audit must be reported. In particular, investigate WSEC, failure,
death, and discharge dates occurring before admission as well as recorded
events or treatments with missing times.
