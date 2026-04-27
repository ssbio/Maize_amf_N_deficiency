# A Genome-Scale Model Tracks How Arbuscular Mycorrhizal Symbiosis Mitigates Maize N Deficiency Stress

This repository contains genome-scale metabolic models (GSMs) and associated analysis files supporting the study:

> **Decouard B., Chowdhury N.B., Saou A., Ampimah N., et al.**
> *"A genome-scale model tracks how arbuscular mycorrhizal symbiosis mitigates maize N deficiency stress."*

## Overview

Arbuscular mycorrhizal fungi (AMF) hold significant potential to reduce nitrogen (N) fertilizer use in agriculture. This study investigates the transcriptomic and metabolomic reprogramming in two maize (*Zea mays*) inbred lines — **B73** and **Lo3** — during AM symbiosis with *Rhizophagus irregularis* under varying nitrogen levels.

Maize transcriptomic data were integrated into the genome-scale metabolic model **iZMA6517** (6,515 genes, 5,228 reactions, 3,007 metabolites) to characterize condition-specific metabolic states, predict metabolite accumulation, and uncover key metabolic shifts — particularly in the **pyrimidine metabolic pathway** — under nitrogen-limited AM symbiosis conditions.

## Experimental Conditions

| Label | Description |
|-------|-------------|
| **HN_CONTROL** | High nitrogen (10 mM NO₃⁻), no AMF inoculation |
| **HN_RHIZO** | High nitrogen (10 mM NO₃⁻), inoculated with *R. irregularis* |
| **LN_CONTROL** | Low nitrogen (1 mM NO₃⁻), no AMF inoculation |
| **LN_RHIZO** | Low nitrogen (1 mM NO₃⁻), inoculated with *R. irregularis* |

## Repository Structure

```
├── Codes/     # GAMS scripts for FBA, Flux Sum Analysis, and Metabolic Bottleneck Analysis
├── GAMS/      # Condition-specific stoichiometric matrices (sij), reaction/metabolite lists,
│              # and flux bounds derived from transcriptomic integration into iZMA6517
│              ├── B73/
│              │   ├── HN_CONTROL/
│              │   ├── HN_RHIZO/
│              │   ├── LN_CONTROL/
│              │   └── LN_RHIZO/
│              └── LO3/
│                  ├── HN_CONTROL/
│                  ├── HN_RHIZO/
│                  ├── LN_CONTROL/
│                  └── LN_RHIZO/
└── Python/    # Condition-specific genome-scale metabolic models in SBML (.xml) format
```

## Computational Methods

- **Flux Balance Analysis (FBA)** — Predicts optimal metabolic flux distributions under each condition
- **Flux Sum Analysis** — Quantifies the total metabolic activity associated with each metabolite
- **Metabolic Bottleneck Analysis** — Identifies rate-limiting reactions constraining metabolic performance

## Dependencies

- [GAMS](https://www.gams.com/) — General Algebraic Modeling System (for FBA and flux analysis scripts)
- [COBRApy](https://opencobra.github.io/cobrapy/) — Python package for constraint-based metabolic modeling

## Correspondence

- **Alia Dellagi** — alia.dellagi@agroparistech.fr (Université Paris-Saclay, INRAE, AgroParisTech, IJPB)
- **Rajib Saha** — rsaha2@unl.edu (University of Nebraska-Lincoln)
