# ai-expo-2026-code

We should hopefully be implementing a novel approach to phase one of drug discovery trials, and integrate it with existing systems. Not gonna be exceptionally complex. We should preferrably be having a working rocq system, and a barebones GNN that can utilise those files (can jump over to LLM based infrastructure if we dont have the time to).

This should be the general template for the README. Ill create a separate folder for the actual research to go in. I'll also create a google docs to make copy pasting easier for citations and references. 

## Knotworking- Inverse Problem for Drug Discovery

## Concept

## Architecture

```
╔═══════════════════════════════════════════════════════════════════════════════╗
║                    DRUG DISCOVERY ARCHITECTURE LAYERS                         ║
╚═══════════════════════════════════════════════════════════════════════════════╝

┌───────────────────────────────────────────────────────────────────────────────┐
│ LAYER 1: DATA INGESTION                                       [DATA LAYER]    │
│                                                                               │
│ Collects patient symptoms, medical history, genome sequences,                 │
│ clinical lab results, imaging data → standardized format                      |
│ Will be using doctor data, due to risks associated (schema provided)          |
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
┌──────────────────────────────────────────────────────────────────────────────┐
│ LAYER 7: MOLECULAR GENERATION (AI)                      [ML LAYER]            │
│                                                                               │
│ Generative models (VAE/diffusion/RL) create billions of candidate             │
│ drug molecules likely to bind target                                          │
└───────────────────────────────────────┬───────────────────────────────────────┘
                                        │
------------------------------------------------------------------------------------------------------------------------
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
║ Cheminformatics pipeline (RDKit) filters molecules to ensure:                 ║         ====================> Novel layers for Phase 1 of Drug Discovery
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
-------------------------------------------------------------------------------------------------------------------------
═══════════════════════════════════════════════════════════════════════════════
                              FEEDBACK LOOPS
═══════════════════════════════════════════════════════════════════════════════


    Layer 9 ──────┐ (Guide generation away from invalid)
                  │
                  └──────→ Layer 7
═══════════════════════════════════════════════════════════════════════════════
                          LAYER DEPENDENCIES
═══════════════════════════════════════════════════════════════════════════════


PARALLELIZABLE:
  • L7 (generation) can produce while L8 (prediction) scores previous batch

═══════════════════════════════════════════════════════════════════════════════
