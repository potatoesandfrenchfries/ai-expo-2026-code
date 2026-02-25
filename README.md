# ai-expo-2026-code

We should hopefully be implementing a novel approach to phase one of drug discovery trials, and integrate it with existing systems. Not gonna be exceptionally complex. We should preferrably be having a working rocq system, and a barebones GNN that can utilise those files (can jump over to LLM based infrastructure if we dont have the time to).

This should be the general template for the README. Ill create a separate folder for the actual research to go in. I'll also create a google docs to make copy pasting easier for citations and references. 

## Knotworking- Inverse Problem for Drug Discovery

## Concept

## Pipeline
```
┌─────────────────────────────────────────────────────────┐
│          DRUG DISCOVERY ARCHITECTURE                    │
└─────────────────────────────────────────────────────────┘

╔═════════════════════════════════════════════════════════╗
║  1. SYMPTOMS → DIAGNOSIS                                ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ Patient reports symptoms
    │
    ├─→ Standard diagnosis (Doctor + Tests)
    │
    └─→ [Alternative] ML-based diagnosis
                ↓
╔═════════════════════════════════════════════════════════╗
║  2. DIAGNOSIS → MECHANISM                               ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ Identify genes/proteins involved
    │
    ├─→ Simulate genome → Predict dysfunctional proteins
    │
    └─→ Identify causal pathway:
        Gene mutation → Protein dysfunction → 
        Pathway disruption → Symptoms
                ↓
╔═════════════════════════════════════════════════════════╗
║  3. MECHANISM → TARGET IDENTIFICATION                   ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ Which protein to target?
    │
    ├─→ Identify rate-limiting proteins
    │
    ├─→ Select druggable targets:
    │   • Has binding pocket
    │   • Not essential for survival
    │
    └─→ [CHOSEN THERAPEUTIC TARGET]
                ↓
╔═════════════════════════════════════════════════════════╗
║  4. TARGET VALIDATION & HIT IDENTIFICATION              ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ ╔════════════════════════════════════════╗
    │   ║ COMPUTATIONAL TARGET CONFIDENCE      ║
    │   ║ Probabilistic causal modeling        ║
    │   ║ Estimate likelihood of reversal      ║
    │   ║ Cross-reference literature & omics   ║
    │   ╚════════════════════════════════════════╝
    │            ↓
    ├─→ Generate billions of candidate molecules
    │   (AI generative models)
    │            ↓
    ├─→ Predict binding affinity (GNN/ML)
    │            ↓
    ├─→ ╔════════════════════════════════════════╗
    │   ║ COMPUTATIONAL DRUG VALIDATION         ║
    │   ║ Filter molecules that satisfy:        ║
    │   ║ • Valid chemistry (bonds, valence)    ║
    │   ║ • Drug-like properties                ║
    │   ║ • Synthesizable                       ║
    │   ║ • No toxic groups                     ║
    │   ╚════════════════════════════════════════╝
    │            ↓
    └─→ [FILTERED VALID MOLECULES + CONFIDENCE SCORES]
                ↓
╔═════════════════════════════════════════════════════════╗
║  5. LEAD OPTIMIZATION                                   ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ Optimize molecular properties:
    │   • Side chains
    │   • Functional groups
    │
    ├─→ Improve drug efficacy:
    │   • Potency
    │   • Solubility
    │   • Safety
    │   • Stability
    │
    ├─→ Preclinical testing:
    │   • Cell assays
    │   • Animal trials
    │   (Preferred over pure ML due to risk)
    │
    └─→ [OPTIMIZED DRUG CANDIDATE]
                ↓
╔═════════════════════════════════════════════════════════╗
║  6. CLINICAL TRIALS & REGULATORY APPROVAL               ║
╚═════════════════════════════════════════════════════════╝
    │
    ├─→ Phase I: Safety (100 people)
    │
    ├─→ Phase II: Efficacy (1000 people)
    │
    ├─→ Phase III: Confirmation (10,000 people)
    │
    ├─→ FDA Approval
    │
    └─→ ★★★ APPROVED DRUG ★★★

═══════════════════════════════════════════════════════════
LEGEND:
═══════════════════════════════════════════════════════════
╔════╗  Main Stage
║    ║  
╚════╝  

╔════╗  Validation Step
║    ║  (Computational / Heuristic checks)
╚════╝

[    ]  Key Milestone/Output

───→   Process Flow

```
## Architecture

```
╔═══════════════════════════════════════════════════════════════════════════════╗
║                    DRUG DISCOVERY ARCHITECTURE LAYERS                         ║
╚═══════════════════════════════════════════════════════════════════════════════╝

┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 1: DATA INGESTION                                       [DATA LAYER]    │
│                                                                               │
│ Collects patient symptoms, medical history, genome sequences,                 │
│ clinical lab results, imaging data → standardized format                      │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 2: DISEASE PREDICTION (ML)                             [ML LAYER]       │
│                                                                               │
│ ML models analyze patient data to predict disease diagnosis                   │
│ and classify disease subtype/severity                                         │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 3: GENOTYPE-TO-PHENOTYPE ASSOCIATION              [DATA/ML LAYER]       │
│                                                                               │
│ Maps genetic variants to protein function using GWAS, QTLs,                   │
│ and sequence-based variant effect predictors (e.g., ESM-variants)             │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 4: CAUSAL PATHWAY MODELING                        [SIMULATION LAYER]    │
│                                                                               │
│ Builds directed graph of protein interactions and biochemical pathways        │
│ showing how mutations lead to disease phenotype                               │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 5: TARGET IDENTIFICATION                          [ML LAYER]            │
│                                                                               │
│ Analyzes pathway graph to identify rate-limiting proteins that are            │
│ druggable (binding pockets, not essential for survival)                       │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
╔═══════════════════════════════════════════════════════════════════════════════╗
║ LAYER 6: TARGET CONFIDENCE SCORING                    [VALIDATION LAYER]      ║
║                                                                               ║
║ Aggregates multi-omics, literature, and causal graph metrics to score target  ║
║ viability; flags high-confidence targets for downstream generation            ║
╚═══════════════════════════════════════┬═══════════════════════════════════════╝
                                        │
                                        ↓
┌──────────────────────────────────────────────────────────────────���────────────┐
│ LAYER 7: MOLECULAR GENERATION (AI)                      [ML LAYER]            │
│                                                                               │
│ Generative models (VAE/diffusion/RL) create billions of candidate             │
│ drug molecules likely to bind target                                          │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 8: PROPERTY PREDICTION (ML)                       [ML LAYER]            │
│                                                                               │
│ Fast neural networks predict binding affinity, solubility, toxicity           │
│ for each candidate (millions per minute)                                      │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
╔═══════════════════════════════════════════════════════════════════════════════╗
║ LAYER 9: COMPUTATIONAL DRUG VALIDATION               [VALIDATION LAYER]       ║
║                                                                               ║
║ Cheminformatics pipeline (RDKit) filters molecules to ensure:                 ║
║  • Valid chemistry/valency       • Drug-likeness (Lipinski rules)             ║
║  • Substructure toxicity alerts  • Basic synthesizability (SA score)          ║
╚═══════════════════════════════════════┬═══════════════════════════════════════╝
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 10: MULTI-FIDELITY SCREENING PIPELINE         [OPTIMIZATION LAYER]      │
│                                                                               │
│ Hierarchical filtering: ML → MD → QM → Experiments                            │
│ to validate surviving candidates                                              │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 11: LEAD OPTIMIZATION ENGINE                   [OPTIMIZATION LAYER]     │
│ If you read this, send me a message on insta saying "Potato". If not, I'll genuinely
crash out considering you haven't read what I have written                                                                              │
│ Generates molecular modifications, predicts property changes,                 │
│ explores optimization landscape (potency/safety/solubility)                   │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 12: SYNTHESIS PLANNING                         [OPTIMIZATION LAYER]     │
│                                                                               │
│ Retrosynthesis algorithms design step-by-step chemical reactions              │
│ to manufacture drug from starting materials                                   │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 13: TRANSCRIPTOMIC PERTURBATION PREDICTION     [ML LAYER]               │
│                                                                               │
│ Predicts gene expression changes (e.g., using Graph-based models)             │
│ in response to drug binding to estimate cellular-level effects                │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 14: PK/PD MODELING                             [SIMULATION LAYER]       │
│                                                                               │
│ Predicts how drug moves through body (ADME) and relates                       │
│ concentration to effect over time                                             │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 15: PATIENT STRATIFICATION MODELING            [ML LAYER]               │
│                                                                               │
│ Clusters patient multi-omics profiles to identify responsive subpopulations   │
│ and predict variable drug efficacy across different demographic groups        │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 16: TOXICITY PREDICTION                        [ML LAYER]               │
│                                                                               │
│ ML models predict organ-specific toxicity (liver, kidney, heart),             │
│ drug-drug interactions, adverse effects                                       │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 17: VIRTUAL COHORT SCREENING                   [STATISTICAL LAYER]      │
│                                                                               │
│ Uses historical trial data and statistical models to optimize trial design,   │
│ predict patient dropout rates, and estimate optimal sample sizes              │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 18: ACTIVE LEARNING CONTROLLER                 [ML LAYER]               │
│                                                                               │
│ Decides which experiments to run next based on model uncertainty              │
│ to maximize information gain, minimize costs                                  │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 19: SEMI-AUTOMATED ASSAY INTEGRATION           [INTEGRATION LAYER]      │
│                                                                               │
│ Formats candidate lists for external CROs / automated HTS systems;            │
│ standardizes returned assay results to update model priors                    │
└────────────────┬──────────────────────┬───────────────────────────────────────┘
                 │                      │
                 │                      │ (Feedback loops)
                 │                      ├────────────────┐
                 │                      │                │
                 ↓                      ↓                ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 20: EVIDENCE & CONFIDENCE MANAGER              [OUTPUT LAYER]           │
│                                                                               │
│ Aggregates all ML predictions, literature cross-references, and assay data    │
│ into a unified confidence score for each lead candidate                       │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 21: MULTI-SCALE ERROR PROPAGATION              [OUTPUT LAYER]           │
│                                                                               │
│ Tracks how uncertainties/approximations compound across scales                │
│ (molecular → organism) to provide confidence bounds                           │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 22: IND PRE-FILING ASSISTANT                   [OUTPUT LAYER]           │
│                                                                               │
│ Uses LLMs to draft initial sections of Investigational New Drug (IND) reports │
│ based on the collected computational and pre-clinical assay evidence          │
└───────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════
                              FEEDBACK LOOPS
═══════════════════════════════════════════════════════════════════════════════

    Layer 19 ─────┐ (Update ML models)
                  │
                  └──────→ Layer 8

    Layer 19 ─────┐ (Refine stratification models)
                  │
                  └──────→ Layer 15

    Layer 9 ──────┐ (Guide generation away from invalid)
                  │
                  └──────→ Layer 7
═══════════════════════════════════════════════════════════════════════════════
                          LAYER DEPENDENCIES
═══════════════════════════════════════════════════════════════════════════════

CRITICAL PATH (must be sequential):
  L1 → L2 → L3 → L4 → L5 → L6 → L7 → L8 → L9 → ... → L22

PARALLELIZABLE:
  • L7 (generation) can produce while L8 (prediction) scores previous batch
  • L13-L17 (simulations & predictions) can run independently on different candidates
  • L15 (stratification) runs across patient populations parallel to drug design
FEEDBACK-DEPENDENT:
  • L18 (active learning) requires L19 (experiments) results
  • L21 (error propagation) requires all prediction layers

═══════════════════════════════════════════════════════════════════════════════
