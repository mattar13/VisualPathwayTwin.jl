# RetinalTwin.jl: Revised Fellowship Application

## Project Overview

RetinalTwin.jl is a biophysically grounded mechanistic simulation of the vertebrate retinal neural circuit, designed as the foundational digital organ component of an autonomous experimental biology platform. The retina is the ideal starting point: it is the best-characterized neural circuit in existence, with complete cell-type catalogs, single-cell transcriptomic atlases, and a non-invasive quantitative functional readout—the electroretinogram (ERG)—that maps directly onto specific cell populations.

Unlike machine learning black boxes, this model is a white box where every parameter is biologically interpretable: ion channel conductances, receptor kinetics, synaptic time constants, and enzymatic rates. This interpretability enables direct mapping between model parameters and genetic/pharmacological perturbations, forming the basis for mechanistic disease modeling and precision medicine.

---

## Platform Architecture

### Layer 1: The Digital Organ (RetinalTwin.jl)

A mechanistic simulation implemented in Julia comprising **49 coupled ordinary differential equations (ODEs)** organized as a single retinal column—the minimal vertical pathway from photoreceptors to ganglion cells. The model includes **10 retinal cell types** governed by **approximately 197 biophysical parameters**, each with direct biological meaning.

#### Cell Types and State Variables

**1. Rod Photoreceptor (21 ODEs, 59 parameters)**
- **Phototransduction cascade** (4 ODEs): Light → rhodopsin (R*) → transducin (T*) → phosphodiesterase (PDE*) → cGMP hydrolysis
  - Parameters: quantum efficiency (λ = 0.67), rhodopsin inactivation rate (kR1 = 400 s⁻¹), PDE activity (kHYDRO = 9000 s⁻¹)
- **Cyclic nucleotide-gated (CNG) channels** with Hill coefficient n = 3
- **Hyperpolarization-activated current (Ih)** with 5-state Markov model for temporal filtering
- **Voltage-gated conductances** (4 ODEs): Kv (delayed rectifier), L-type Ca²⁺, KCa (Ca²⁺-activated K⁺), ICl(Ca)
- **Inner segment Ca²⁺ dynamics** (6 ODEs): Fast/slow compartments with dual-affinity buffers (low: Kd = 0.2 μM; high: Kd = 90 μM)
- **Graded glutamate release**: Voltage-dependent sigmoid with 10 ms time constant

**Interpretability example**: Setting the dark current conductance `iDARK → 0` simulates photoreceptor degeneration and eliminates the ERG a-wave, matching retinitis pigmentosa phenotypes.

**2. ON Bipolar Cell (4 ODEs, 20 parameters)**
- **mGluR6 metabotropic cascade** (1 ODE): Glutamate activates Gαo → inhibits TRPM1 channels
  - **Sign inversion**: High glutamate (dark) → TRPM1 closed → hyperpolarized; Low glutamate (light) → TRPM1 open → depolarized
  - Parameters: TRPM1 conductance (g_TRPM1 = 2.0 nS), cascade time constant (tau_mGluR6 = 50 ms)
- **Morris-Lecar excitability**: Ca²⁺ (depolarizing) and K⁺ (repolarizing) conductances with voltage-dependent gating
- **Graded glutamate output** to downstream amacrine and ganglion cells

**Disease modeling**: Setting `g_TRPM1 → 0` simulates complete congenital stationary night blindness (cCSNB) and produces an **electronegative ERG** (no b-wave, preserved a-wave), matching the clinical phenotype of patients with TRPM1 null mutations.

**3. OFF Bipolar Cell (4 ODEs, 19 parameters)**
- **Ionotropic glutamate receptors (iGluR)** with fast kinetics (tau = 3 ms)
- **Sign-preserving pathway**: High glutamate → depolarization (opposite of ON bipolar)

**4. Horizontal Cell (3 ODEs, 16 parameters)**
- **Lateral inhibition**: Provides negative feedback to photoreceptors via GABA or pH modulation
- **Gap junction coupling** (parameter: g_gap) forms spatially extended network

**5-7. Amacrine Cells (9 ODEs total)**
- **A2 (AII) Amacrine** (3 ODEs, 16 parameters): Glycinergic, narrow-field
  - Critical for **oscillatory potential (OP) generation** via reciprocal inhibition with GABAergic amacrines
  - Fast dynamics (low capacitance C_m = 5 pF, high phi = 0.5) enable 100-160 Hz oscillations
- **GABAergic Amacrine** (3 ODEs, 16 parameters): Forms reciprocal network with A2
- **Dopaminergic Amacrine** (3 ODEs, 16 parameters): Slow modulatory signaling (tau_DA = 200 ms)

**8. Ganglion Cell (2 ODEs, 16 parameters)**
- **Action potential generation**: Morris-Lecar bistable dynamics
- **Output neuron**: Integrates excitatory (bipolar) and inhibitory (amacrine) inputs

**9. Müller Glial Cell (4 ODEs, 9 parameters)**
- **K⁺ spatial buffering**: Inward-rectifying K⁺ channels (Kir4.1) at endfoot and stalk compartments
  - Parameters: g_Kir_end = 15 nS, g_Kir_stalk = 8 nS, diffusion time constant tau_K = 100 ms
- **Glutamate uptake**: Excitatory amino acid transporters (EAAT1/2) with Michaelis-Menten kinetics (Km = 10 μM, Vmax = 0.5 μM/ms)
- **ERG contribution**: Generates P3 component via transmembrane K⁺ currents

**10. Retinal Pigment Epithelium (RPE) (2 ODEs, 10 parameters)**
- **Subretinal K⁺ sensing**: Apical Kir channels detect photoreceptor dark current modulation
- **c-wave generation**: Light → reduced photoreceptor K⁺ efflux → [K⁺]subretinal drops → RPE hyperpolarizes
  - Very slow dynamics: tau_RPE = 1000 ms (c-wave peak at 2-5 seconds)

#### Electroretinogram (ERG) Output

The ERG is computed as a weighted sum of transmembrane currents from photoreceptors, bipolar cells, Müller glia, and RPE:

```
ERG(t) = w_PR · I_PR(t) + w_BC · I_BC(t) + w_Muller · I_Muller(t) + w_RPE · I_RPE(t)
```

**Components**:
- **a-wave** (0-30 ms): Photoreceptor hyperpolarization (negative deflection)
- **b-wave** (30-100 ms): ON bipolar depolarization (positive deflection)
- **Oscillatory potentials (OPs)**: Amacrine network oscillations (100-160 Hz, 75-300 Hz bandpass)
- **c-wave** (2-5 s): RPE hyperpolarization

#### Parameter Interpretability

All **197 parameters** are stored in CSV files with biologically meaningful names, default values, and prior bounds for Bayesian inference:

**Example (photoreceptor_params.csv)**:
```csv
Key,Value,LowerBounds,UpperBounds,DEFAULT
kHYRDO,9000,0,10000,9000      # PDE hydrolysis rate (s⁻¹)
gKV,2.0,0,10,2.0               # Voltage-gated K⁺ conductance (nS)
tau_Glu,10.0,1,50,10.0         # Glutamate release time constant (ms)
V_Glu_half,-40.0,-50,-30,-40   # Glutamate release half-activation (mV)
```

**Disease perturbation mapping**:
| Parameter Perturbation | Biological Mutation | ERG Phenotype |
|------------------------|---------------------|---------------|
| `g_TRPM1 = 0` | TRPM1 null (cCSNB) | Electronegative ERG (no b-wave) |
| `kHYDRO ↓ 50%` | PDE6 mutation | Reduced a-wave amplitude |
| `iDARK ↓ 90%` | Photoreceptor degeneration (RP) | Flat ERG |
| `g_Kir ↓ 50%` | Müller dysfunction | Reduced b-wave, abnormal OPs |

---

### Layer 2: Parameter Estimation via Staged Bayesian Inference

**Challenge**: Fitting 197 parameters simultaneously from a single ERG waveform is ill-posed (underdetermined system).

**Solution**: Exploit the known temporal and cellular origins of ERG components to fit parameters sequentially.

#### Staged Fitting Pipeline

**Stage 1: Photoreceptor Parameters (59 parameters)**
- **Data**: Isolated a-wave (first 20-30 ms, before b-wave onset)
- **Method**: Bayesian inference via No-U-Turn Sampler (NUTS) using isolated photoreceptor model
- **Outputs**:
  - Point estimates for phototransduction kinetics (kHYDRO, kREC, lambda)
  - Full posterior distributions for uncertainty quantification
  - Identifiable parameters: dark current (iDARK), Ca²⁺ buffer capacities (Bl, Bh)

**Stage 2: Bipolar Parameters (39 parameters)**
- **Data**: b-wave amplitude and time-to-peak (30-100 ms)
- **Method**: Fix photoreceptor parameters from Stage 1, fit ON/OFF bipolar parameters
- **Key parameters**:
  - mGluR6 cascade gain (alpha_mGluR6 = 2.0)
  - TRPM1 conductance (g_TRPM1 = 2.0 nS)
  - Bipolar membrane time constant (C_m / g_L)

**Stage 3: Amacrine Parameters (48 parameters)**
- **Data**: Oscillatory potentials (100-160 Hz bandpass filtered ERG)
- **Method**: Fit reciprocal A2-GABA network parameters to OP frequency and amplitude
- **Key parameters**:
  - Glycine/GABA release rates (alpha_Gly, alpha_GABA)
  - Network coupling strengths (g_Gly_inh, g_GABA_inh)
  - Fast membrane dynamics (phi = 0.5 for high-frequency oscillations)

**Stage 4: Müller and RPE Parameters (19 parameters)**
- **Data**: Slow components (c-wave peak at 2-5 s, P3 component)
- **Method**: Fit K⁺ buffering dynamics and slow RPE time constants
- **Key parameters**:
  - Kir conductances (g_Kir_end, g_Kir_stalk)
  - RPE time constant (tau_RPE = 1000 ms)

#### Prior Construction

Priors are informed by:
1. **Literature values**: Ion channel conductances from patch-clamp studies
2. **Transcriptomic constraints**: Gene expression levels constrain max conductances
3. **Cross-species conservation**: Parameters consistent across mouse/human/primate

**Example prior specification**:
```julia
using Turing

@model function fit_photoreceptor(erg_data, time)
    # Informative priors from literature
    iDARK ~ truncated(Normal(5040, 500), 3000, 6000)  # pA
    kHYDRO ~ truncated(Normal(9000, 1000), 5000, 12000)  # s⁻¹

    # Simulate model
    prob = ODEProblem(rod_model!, u0, tspan, [iDARK, kHYDRO, ...])
    sol = solve(prob, Rodas5())

    # Likelihood: match ERG a-wave
    a_wave_predicted = extract_a_wave(sol)
    erg_data ~ MvNormal(a_wave_predicted, σ_noise^2 * I)
end

# Sample posterior
chain = sample(fit_photoreceptor(erg_data, t), NUTS(), 2000)
```

**Outputs**:
- **Point estimates**: Maximum a posteriori (MAP) or posterior mean
- **Uncertainty quantification**: 95% credible intervals for each parameter
- **Correlations**: Posterior correlation matrix identifies co-varying parameters

---

### Layer 3: AI Experimental Designer (Proposed Integration with FutureHouse)

An LLM-powered reasoning agent that synthesizes three evidence channels to perform differential diagnosis, disease staging, and experiment design:

#### Evidence Channels

**1. Fitted Model Parameters**
- Input: Posterior distributions from Layer 2 Bayesian inference
- Example: "Patient ERG shows `g_TRPM1_posterior ~ Normal(0.5, 0.2)` vs control `Normal(2.0, 0.1)` → 75% TRPM1 conductance loss"

**2. Clinical and Genomic Data**
- Input: Patient metadata (age, family history, visual acuity, fundus images)
- Genomic variants from whole-exome sequencing (WES) or targeted panel
- Example: "TRPM1 p.R684C variant (missense) + reduced b-wave → likely pathogenic"

**3. Structured Disease Knowledge Base (RAG over Literature)**
- Indexed corpus: ~5,000 ERG and retinal disease publications
- Query: "Retrieve ERG phenotypes associated with TRPM1 mutations"
- Returns: "Electronegative ERG (b/a ratio < 1.0), normal scotopic threshold, no progressive degeneration"

#### Core Capabilities

**Differential Diagnosis**:
```
Input: Patient ERG shows reduced b-wave (40% of normal), normal a-wave, no OPs
Agent reasoning:
1. Query disease KB: "reduced b-wave + normal a-wave" → candidates: cCSNB, XLRS, rod-cone dystrophy
2. Check model fits: TRPM1 conductance posterior: 0.8 nS (60% reduction)
3. Check genomic data: TRPM1 c.2051G>A (p.R684H) variant
4. Integrate evidence: High confidence for TRPM1-mediated cCSNB
Output: Diagnosis = "Complete congenital stationary night blindness (TRPM1 mutation)",
        Confidence = 0.92
```

**Disease Staging and Prognosis**:
- Track parameter evolution over time (longitudinal ERGs)
- Example: "Photoreceptor dark current declined 15%/year → predict legal blindness in 8 years"

**Active Learning for Experiment Design**:
The agent identifies parameters with high posterior uncertainty and proposes experiments to resolve ambiguity:

**Example**:
```
Current uncertainty: Ca²⁺ buffer parameters (Hb1, Hb2) have wide posteriors (CV > 50%)
Proposed experiment: "Pharmacological Ca²⁺ buffer block (BAPTA-AM) + ERG recording"
Expected information gain: Reduces Ca²⁺ parameter posterior variance by 70%
Justification: BAPTA specifically blocks high-affinity buffers (Hb) → disambiguates Hb1 vs Lb1
```

**Experiment Prioritization**:
The agent ranks experiments by:
1. **Expected information gain** (KL divergence between prior and predicted posterior)
2. **Feasibility** (cost, ethical constraints, available models)
3. **Clinical relevance** (does it change treatment decisions?)

#### Required AI Infrastructure (FutureHouse Support)

1. **LLM with extended context** (e.g., Claude 3.5 with 200K context) for:
   - Parsing complex ERG waveforms and clinical reports
   - Reasoning over disease knowledge base (RAG queries)
   - Generating natural language explanations of model predictions

2. **Vector database** for semantic search over literature:
   - Embed 5,000 ERG/retinal disease papers
   - Query: "What ERG phenotypes are associated with NYX gene mutations?"
   - Return: Ranked passages from Bech-Hansen et al. (2000), "Complete CSNB with negative ERG"

3. **Experiment simulation engine**:
   - Agent generates proposed perturbations (e.g., "Set g_TRPM1 = [0, 0.5, 1.0, 2.0] nS")
   - Simulate ERG outputs for each condition
   - Compute Fisher information to estimate parameter resolution

4. **Multi-modal integration**:
   - Combine ERG data (time series) + genomic variants (structured data) + fundus images (vision)
   - Example: "Fundus shows bone spicules (pigmentary degeneration) + flat ERG → advanced RP"

---

### Layer 4: Closed-Loop Wet-Lab Validation (Blackshaw & Peachey Labs)

#### Validation Data Sources

**Blackshaw Lab (Johns Hopkins)**:
- Single-cell RNA-seq atlases of retinal cell types (mouse, human)
- CRISPR knockouts: Targeted deletion of ion channels, receptors, signaling proteins
- Transcription factor perturbations: Cell-fate specification experiments
- **Use case**: Validate parameter priors from gene expression levels

**Peachey Lab (Case Western Reserve University)**:
- Clinical ERG database:
  - X-linked retinoschisis (XLRS): n ≈ 20 patients
  - Retinitis pigmentosa (RP): n ≈ 30 patients (RHO, PDE6, RPGR mutations)
  - Congenital stationary night blindness (CSNB): n ≈ 10 patients (TRPM1, GRM6, NYX mutations)
- Longitudinal recordings for disease progression tracking
- **Use case**: Fit model to patient ERGs, validate predictions against genotype

#### Iterative Refinement Loop

1. **Initial fit**: Use literature priors, fit to control ERG dataset (n = 50 healthy subjects)
2. **Validate on knockouts**: Predict ERG phenotype of TRPM1⁻/⁻ mouse → compare to Peachey lab recording
3. **Identify discrepancies**: Model predicts OP amplitude 120 μV, data shows 80 μV → refine amacrine parameters
4. **Update priors**: Incorporate knockout data as additional constraints on parameter space
5. **Re-fit and test**: Improved model now predicts CSNB patient ERGs with R² = 0.91 (vs 0.74 before)

#### Uncertainty-Guided Validation

The model's posterior uncertainty directly guides which experiments to prioritize:

**High-uncertainty parameters → targeted validation**:
- If `tau_mGluR6` posterior is wide (25-75 ms range), propose:
  - **Experiment**: Voltage-clamp bipolar cells, measure mGluR6 cascade kinetics via puff application of glutamate
  - **Expected outcome**: Resolves time constant to ±5 ms
  - **Impact**: Improves b-wave time-to-peak predictions by 30%

---

## Current Implementation Status

### Completed (Phase 1)
- ✅ All 10 cell types implemented with biophysically grounded ODEs
- ✅ 197 parameters stored in version-controlled CSV files
- ✅ Photoreceptor + ON bipolar working example (25 ODEs)
- ✅ Glutamate release dynamics with voltage-dependent kinetics
- ✅ mGluR6 sign-inverting synapse validated (produces correct depolarization)
- ✅ Morris-Lecar framework for all neurons (horizontal, bipolar, amacrine, ganglion)
- ✅ Müller K⁺ buffering model with spatial compartments
- ✅ RPE c-wave generation mechanism

### In Progress (Phase 2)
- ⏳ Full retinal column assembly (49 ODEs with all cell types coupled)
- ⏳ ERG field potential calculation (weighted sum of transmembrane currents)
- ⏳ Oscillatory potential extraction (bandpass filtering + peak detection)
- ⏳ Parameter fitting pipeline (Bayesian inference via Turing.jl)

### Planned (Phase 3+)
- 🔲 Cone photoreceptor pathway (requires separate phototransduction cascade)
- 🔲 Spatial grid expansion (2D array of columns for receptive field dynamics)
- 🔲 Gap junction coupling (A2-cone bipolar, horizontal network)
- 🔲 LLM experimental designer (requires FutureHouse AI infrastructure)
- 🔲 Integration with wet-lab validation data (Blackshaw/Peachey labs)

---

## Why the Retina?

The retina is uniquely positioned as the **ideal first neural circuit** for an autonomous experimental biology platform:

1. **Complete cell-type catalog**: All ~100 retinal neuron types characterized (Macosko et al. 2015)
2. **Non-invasive readout**: ERG clinically accessible, maps directly to cellular activity
3. **Genetic tractability**: >250 Mendelian retinal diseases with known causative genes
4. **Translational relevance**: 10 million patients with inherited retinal dystrophies worldwide
5. **Conserved across species**: Mouse models directly inform human disease mechanisms

**Contrast with other neural circuits**:
- **Cortex**: 1000+ cell types, no direct non-invasive readout, limited disease models
- **Spinal cord**: Difficult genetic access, requires invasive recordings
- **Hippocampus**: Behavioral readouts confounded by distributed networks

**Scalability**: Once validated on retina, the platform can extend to:
- Other sensory systems (cochlea → auditory brainstem response)
- Motor circuits (spinal cord → electromyography)
- Central circuits (cortex → EEG/MEG), though interpretability will decrease

---

## Technical Innovations

1. **Hierarchical Bayesian fitting**: Exploit known ERG-to-cell-type mapping to decouple parameter estimation across cell types (avoids 197-dimensional optimization)

2. **Interpretable parameterization**: Every parameter is a physical quantity, enabling direct mapping to genetic mutations and pharmacological interventions

3. **Uncertainty-guided experiment design**: Posterior variance identifies which parameters are poorly constrained, enabling active learning for biology

4. **Cross-scale integration**: Molecular (ion channels) → cellular (membrane dynamics) → tissue (ERG) → clinical (diagnosis) within a single unified model

5. **Modular architecture**: Cell models can be validated in isolation, then assembled into circuit (reduces debugging complexity)

---

## Expected Impact

### Scientific Impact
- **Mechanistic disease models**: Replace phenomenological clinical diagnostic criteria with quantitative biophysical parameter perturbations
- **Personalized medicine**: Individual patient ERGs → fitted parameter sets → tailored treatment (e.g., gene therapy dosing)
- **Accelerate target discovery**: Model identifies which ion channels/receptors are most sensitive to perturbation → prioritize drug development

### Clinical Impact
- **Early diagnosis**: Detect parameter deviations before symptomatic vision loss (e.g., 20% photoreceptor loss detectable via ERG a-wave)
- **Prognosis**: Predict disease progression rates from longitudinal parameter trajectories
- **Treatment monitoring**: Quantify gene therapy efficacy via parameter recovery (e.g., "TRPM1 conductance restored from 0.5 nS → 1.8 nS post-treatment")

### Methodological Impact
- **Blueprint for other organ systems**: Retinal digital twin provides template for heart (ECG), brain (EEG), muscle (EMG)
- **AI-driven biology**: Demonstrates closed-loop integration of mechanistic models + LLM reasoning + wet-lab validation
- **Open science**: All code, parameters, and fitted models publicly available for community validation and extension

---

## Timeline and Milestones

**Months 1-3**: Complete full retinal column coupling (49 ODEs), validate ERG waveform morphology
**Months 4-6**: Implement staged Bayesian fitting pipeline, fit to control ERG database (n=50)
**Months 7-9**: Validate disease models (CSNB, RP, XLRS) against clinical ERG data
**Months 10-12**: Integrate LLM experimental designer (FutureHouse collaboration), demonstrate active learning loop

**Year 2**: Spatial grid expansion, cone pathway, wet-lab validation with Blackshaw/Peachey data

---

## Requested Support from FutureHouse

1. **AI engineering**: LLM infrastructure for experimental reasoning agent (Claude 3.5 API, RAG pipeline)
2. **Compute resources**: Bayesian sampling (10⁴-10⁵ ODE solves per fit), parallelized across patient cohorts
3. **Vector database**: Literature indexing (5,000 papers → embeddings → semantic search)
4. **Collaboration network**: Connections to clinical ERG databases, genetic consortia
5. **Mentorship**: Guidance on productionizing research code → deployable diagnostic tool

---

## References

1. Robson & Frishman (2014), "The origin of the electroretinogram," *Prog Retin Eye Res*
2. Macosko et al. (2015), "Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets," *Cell*
3. Lamb & Pugh (1992), "A quantitative account of the activation steps involved in phototransduction," *J Physiol*
4. Shiells & Falk (1990), "Glutamate receptors of rod bipolar cells are linked to a cyclic GMP cascade," *Proc R Soc Lond B*
5. Newman et al. (1984), "Control of extracellular potassium levels by retinal glial cell K⁺ siphoning," *Science*

---

**Contact**: Dr. Matt Tarchick
**Repository**: https://github.com/mattar13/RetinalTwin.jl
**Documentation**: See `docs/retinal_digital_twin_spec.md` for full mathematical specification
