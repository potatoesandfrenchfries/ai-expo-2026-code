# ai-expo-2026-code


This should be the general template for the README. Ill create a separate folder for the actual research to go in. I'll also create a google docs to make copy pasting easier for citations and references. 

## Knotworking- Inverse Problem for Drug Discovery

## Concept

## Pipeline

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
    │   ║ FORMAL VERIFICATION (Rocq/Lean)      ║
    │   ║ Build causal model of disease        ║
    │   ║ Prove: modulating target reverses    ║
    │   ║        disease phenotype              ║
    │   ╚════════════════════════════════════════╝
    │            ↓
    ├─→ Generate billions of candidate molecules
    │   (AI generative models)
    │            ↓
    ├─→ Predict binding affinity (GNN/ML)
    │            ↓
    ├─→ ╔════════════════════════════════════════╗
    │   ║ FORMAL VERIFICATION                   ║
    │   ║ Prove molecules satisfy:              ║
    │   ║ • Valid chemistry (bonds, valence)    ║
    │   ║ • Drug-like properties                ║
    │   ║ • Synthesizable                       ║
    │   ║ • No toxic groups                     ║
    │   ╚════════════════════════════════════════╝
    │            ↓
    └─→ [FILTERED VALID MOLECULES + PROOF CERTIFICATES]
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

╔════╗  Formal Verification Step
║    ║  (Mathematical proof in Rocq/Lean)
╚════╝

[    ]  Key Milestone/Output

───→   Process Flow
## Architecture

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
│ LAYER 3: GENOTYPE-TO-PHENOTYPE SIMULATION               [SIMULATION LAYER]    │
│                                                                               │
│ Simulates how genetic variants affect protein function                        │
│ and cellular behavior to identify molecular disease mechanism                 │
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
║ LAYER 6: TARGET VALIDATION WITH FORMAL VERIFICATION   [VERIFICATION LAYER]    ║
║                                                                               ║
║ Proves mathematically (Rocq/Lean) that modulating target reverses disease;    ║
║ outputs proof certificate guaranteeing target validity                        ║
╚═══════════════════════════════════════┬═══════════════════════════════════════╝
                                        │
                                        ↓
┌───────────────────────────────────────────────────────────────────────────────┐
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
║ LAYER 9: FORMAL CHEMICAL VERIFICATION                [VERIFICATION LAYER]     ║
║                                                                               ║
║ Theorem prover checks molecules satisfy:                                      ║
║  • Valid bonds/angles/valence  • Drug-likeness (Lipinski)                     ║
║  • No toxic substructures      • Synthesizable                                ║
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
│                                                                               │
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
│ LAYER 13: CELLULAR SIMULATION                        [SIMULATION LAYER]       │
│                                                                               │
│ Simulates drug effects on cells (signaling, gene expression, phenotype)       │
│ to predict efficacy and toxicity                                              │
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
│ LAYER 15: DIGITAL TWIN GENERATION                    [SIMULATION LAYER]       │
│                                                                               │
│ Creates personalized computational models of individual patients              │
│ (genome, physiology) to simulate individual responses                         │
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
│ LAYER 17: CLINICAL TRIAL SIMULATION                  [SIMULATION LAYER]       │
│                                                                               │
│ Simulates thousands of virtual patients to predict trial outcomes,            │
│ dosing, patient stratification, sample sizes                                  │
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
│ LAYER 19: EXPERIMENTAL INTEGRATION                   [INTEGRATION LAYER]      │
│                                                                               │
│ Interfaces with robotic labs to execute selected experiments;                 │
│ feeds results back to update all models                                       │
└────────────────┬──────────────────────┬───────────────────────────────────────┘
                 │                      │
                 │                      │ (Feedback loops)
                 │                      ├────────────────┐
                 │                      │                │
                 ↓                      ↓                ↓
┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 20: PROOF CERTIFICATE MANAGER                  [OUTPUT LAYER]           │
│                                                                               │
│ Stores formal proofs from verification; generates reports showing             │
│ which properties are mathematically guaranteed vs predicted                   │
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
│ LAYER 22: REGULATORY DOCUMENTATION GENERATOR         [OUTPUT LAYER]           │
│                                                                               │
│ Compiles computational evidence, experimental data, proof certificates,       │
│ simulations into FDA-compliant submission                                     │
└───────────────────────────────────────────────────────────────────────────────┘

═══════════════════════════════════════════════════════════════════════════════
                              FEEDBACK LOOPS
═══════════════════════════════════════════════════════════════════════════════

    Layer 19 ─────┐ (Update ML models)
                  │
                  └──────→ Layer 8

    Layer 19 ─────┐ (Refine digital twins)
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
  • L13-L17 (simulations) can run independently on different candidates
  • L15 (digital twins) one per patient, fully parallel

FEEDBACK-DEPENDENT:
  • L18 (active learning) requires L19 (experiments) results
  • L21 (error propagation) requires all simulation layers

═══════════════════════════════════════════════════════════════════════════════

## Existing research

## Hypothetical Advantages

## Requirements

## References

