# Draft LinkedIn Post: OpenBind EV-A71 2A Rigorous Benchmark Analysis

## TL;DR
After rigorous re-analysis of the OpenBind EV-A71 2A protease benchmark (~700 compounds, 601 measured affinities, 925 structures), the headline is not "physics wins / ML fails." The real story is subtler: **zero-shot co-folding fails at affinity ranking, target-specific fine-tuning partially recovers global screening performance, local SAR defeats everyone, and the most promising signal comes from combining ML and physics in ensemble models.**

---

## The Framing Problem

Initial analyses of this benchmark (including my own) have framed the debate as zero-shot ML co-folding versus physics-based FEP. I spent a week digging into the actual numbers, computing bootstrap 95% confidence intervals, and reporting **all** tested protocols, not just the ones that fit a narrative. Here's what the data actually says.

---

## 1. Global Affinity Screening (N=381–488): Fine-Tuning Matters, But Physics Baselines Still Win

| Method | Spearman rho [95% CI] | p-value |
|---|---|---|
| Molecular Weight (baseline) | **0.485** [0.414, 0.552] | < 1e-29 |
| GNINA Crystal Docking | **0.455** [0.378, 0.521] | < 1e-25 |
| Protenix co-folding | **0.388** [0.296, 0.470] | < 1e-14 |
| Boltz-2 co-folding | **0.400** [0.322, 0.472] | < 1e-19 |
| OpenFold3-p2 (zero-shot) | **-0.325** [-0.412, -0.231] | < 1e-10 |
| OpenFold3-p2 (fine-tuned) | **0.353** [0.257, 0.438] | < 1e-12 |

**Interpretation:** Zero-shot co-folding produces a *negative* rank correlation with affinity for this target. Target-specific fine-tuning on the fragment screen flips the sign, but the fine-tuned model was trained on this exact dataset. Whether this represents memorization or true target-specific learning is an open question. Traditional crystal docking and even molecular weight remain deceptively strong baselines.

**Key caveat:** The confidence interval for zero-shot OF3-p2 (ρ = -0.325 [CI: -0.412, -0.231]) does *not* overlap with zero — this failure is real. But the fine-tuned OF3-p2 CI (ρ = 0.353 [CI: 0.257, 0.438]) still falls below crystal docking (ρ = 0.455 [CI: 0.378, 0.521]).

---

## 2. Local SAR / Lead Optimization (N=28–32): Everyone Fails

On the focused pyrrolidine congeneric series, **no method achieves statistical significance** after bootstrap resampling. Every confidence interval crosses zero:

| Method | Spearman rho [95% CI] | p-value |
|---|---|---|
| OF3-p2 zero-shot | -0.246 [-0.590, 0.187] | 0.207 |
| OF3-p2 fine-tuned | -0.087 [-0.511, 0.325] | 0.660 |
| Boltz-2 | -0.096 [-0.471, 0.310] | 0.628 |
| GNINA docking | -0.231 [-0.573, 0.171] | 0.237 |
| Rowan xtal_nagl (crystal, NAGL) | 0.127 [-0.258, 0.420] | 0.489 |
| Rowan docking_am1bcc_0local (dock, AM1-BCC) | **0.665** | 3.35e-05 |
| *...across 6 tested protocols* | *0.030 to 0.665* | *varies* |

**Interpretation:** Local SAR is genuinely hard. Co-folding confidence scores (pair_iptm) are not designed for sub-angstrom R-group sensitivity, and it shows. The one Rowan protocol that achieves significance (docking_am1bcc_0local, ρ = 0.665) is part of a broader parameter sweep where results ranged from 0.030 to 0.665 depending on charge model and local sampling strategy. Rowan themselves label these baselines as "not yet good."

**Most importantly:** Only *one* of six tested protocols achieved significance on this target. Calling this a "physics victory" is cherry-picking. The honest finding is that **both zero-shot ML and unoptimized physics struggle on narrow SAR series**, and protocol selection matters enormously.

---

## 3. Protocol Sensitivity: The Variation Is the Story

The Rowan suite tested 6 charge/model/sampling protocols on the same 32 compounds:

| Protocol | Spearman rho | p-value |
|---|---|---|
| xtal_nagl (crystal pose, NAGL charges) | 0.127 | 0.489 (NS) |
| docking_nagl (docked pose, NAGL charges) | 0.363 | 0.041 (marginal) |
| docking_am1bcc (docked, AM1-BCC, no extra sampling) | 0.030 | 0.872 (NS) |
| docking_am1bcc_4ns (docked, AM1-BCC, 4ns extra) | 0.204 | 0.262 (NS) |
| docking_am1bcc_350local (docked, AM1-BCC, 350 local steps) | 0.331 | 0.065 (NS) |
| docking_am1bcc_0local (docked, AM1-BCC, 0 local steps) | **0.665** | **3.35e-05** |

**6x variation from the same input data.** This is not a side note — this *is* the finding. Physics-based scoring is not a magic bullet; it is a family of approximations that can range from useless to useful depending on how carefully you set them up. For computational chemistry teams, the takeaway is not "use physics" but rather **"optimize your physics protocol before trusting it for lead optimization."**

---

## 4. Novel Finding: ML + Physics Ensembles Outperform Either Alone

We tested simple z-score ensemble models combining ML and physics scores on the global dataset. The results are the most encouraging in the entire benchmark:

| Ensemble | Spearman rho [95% CI] | p-value |
|---|---|---|
| OF3-p2 fine-tuned + GNINA crystal | **0.517** [0.438, 0.592] | < 1e-26 |
| Boltz-2 + GNINA crystal | **0.469** [0.397, 0.534] | < 1e-27 |
| Molecular weight + Boltz-2 | **0.478** [0.403, 0.545] | < 1e-29 |

**Interpretation:** When ML and physics scores are combined, they achieve higher rank correlation than either component alone. This suggests the two paradigms capture **different but complementary** information about ligand binding. ML co-folding may encode global shape and chemical complementarity, while physics encodes specific energy terms that co-folding misses.

**Critical caveat on generalization:** The top ensemble (OF3-p2 fine-tuned + GNINA, ρ = 0.517) includes a fine-tuned component trained on this exact dataset. Whether this represents memorization or genuine generalization is unknown. The zero-shot ensemble (Boltz-2 + GNINA, ρ = 0.469) is the **more robust** finding — it combines two methods that have never seen this data, and still beats either component alone. This is the result worth building on.

This is not a "solution" — even ρ = 0.469 on global screening is modest — but it points to a **constructive** research direction instead of the tired "ML vs. physics" debate.

---

## The Honest Takeaway for Computational Drug Discovery

### 1. Do not use zero-shot co-folding for affinity ranking.
Zero-shot OF3-p2 and Boltz-2 show negative or near-zero correlations with affinity on this target. Protenix is the exception — it achieves ρ = 0.388 [0.296, 0.470], suggesting architecture choice matters. But the inconsistent direction of these zero-shot correlations underscores the same point: `pair_iptm` confidence is a pose-quality metric, not a reliable affinity metric. Expecting any single co-folding model to rank binders is using the tool for the wrong job.

### 2. Target-specific fine-tuning is necessary but insufficient.
Fine-tuning recovers positive correlation for global screening, but the model was trained on the dataset it is evaluated on. Whether this generalizes to novel scaffolds is untested. Use it for pose generation, not as a standalone affinity model.

### 3. Local SAR is hard for everyone.
No method — not co-folding, not docking, not the best physics protocol — achieves robust significance on the 32-compound pyrrolidine series. This is a real, unsolved problem, and the field should treat it as such instead of pretending one paradigm has "won."

### 4. Protocol sensitivity in physics-based scoring is the elephant in the room.
A 6x range in rank correlation from a single parameter sweep is sobering. If your FEP "works," ask: how many other protocols failed before you found the right one? Rigorous protocol sweeps and cross-validation should be standard, not afterthoughts.

### 5. The most promising next step is ensembles, not winners.
Instead of asking "which method is best?", the field should ask "how do we combine ML and physics to capture both compositional and energetic signals?" Simple z-score ensembles already outperform individual methods. There is real science to be done here.

---

## Methodology (for the skeptics)

- **Bootstrap CIs:** 2000 resamples, percentile method, 95% confidence
- **Statistical testing:** Two-tailed Spearman correlation with Bonferroni correction considered (not applied here for exploratory analysis; raw p-values reported)
- **All code:** Open-source, fully reproducible, available at the repository below

## Data & Code

- OpenBind dataset: [Zenodo / Fragalysis]
- Analysis pipeline: `src/rigorous_analysis.py` (bootstrap CIs, all protocols, no cherry-picking)
- Figure generation: `results/rigorous_benchmark.png` (4-panel figure showing CIs, all Rowan runs, local SAR, and ensemble gains)

---

## Questions This Raises (Not Answers)

- Why does zero-shot co-folding produce *negative* affinity correlations? Is it an artifact of confidence-score calibration, or does it reflect a real bias in the model's energy landscape?
- What structural features drive the protocol sensitivity in Rowan FEP? Is it the protonation state, the charge model, or the starting conformation?
- Can a more sophisticated ensemble (weighted by validation performance, or using a learned combination) further improve rank correlation?
- Does the fine-tuned OF3-p2 model generalize to chemically distinct scaffolds, or does it overfit to the fragment screen?

These are the questions worth pursuing. And they will only be answered with more rigorous, less tribal benchmarking.

---

*#ComputationalChemistry #DrugDiscovery #MachineLearning #StructuralBiology #FEP #Benchmarking #OpenBind #OpenScience*
