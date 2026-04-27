# Maize AMF Nitrogen Deficiency — Metabolic Modeling

This repository contains genome-scale metabolic models and flux balance analysis (FBA) files for studying the effect of arbuscular mycorrhizal fungi (AMF) colonization on maize (*Zea mays*) metabolism under nitrogen-deficient conditions.

## Genotypes
- **B73** — A widely used maize reference inbred line
- **LO3** — A maize inbred line used for comparative analysis

## Experimental Conditions
| Condition | Description |
|-----------|-------------|
| HN_CONTROL | High nitrogen, no AMF inoculation |
| HN_RHIZO | High nitrogen, with AMF/rhizosphere inoculation |
| LN_CONTROL | Low nitrogen, no AMF inoculation |
| LN_RHIZO | Low nitrogen, with AMF/rhizosphere inoculation |

## Repository Structure

```
├── Codes/          # GAMS scripts for FBA, flux sum, and metabolic bottleneck analysis
├── GAMS/           # Stoichiometric matrices, reaction/metabolite lists, and flux bounds per condition
└── Python/         # Genome-scale metabolic models in SBML (.xml) format
```

## Analysis
- **Flux Balance Analysis (FBA)** — Predicts metabolic flux distributions under each condition
- **Flux Sum Analysis** — Quantifies overall metabolic activity per metabolite
- **Metabolic Bottleneck Analysis** — Identifies limiting reactions in the metabolic network

## Dependencies
- [GAMS](https://www.gams.com/) (General Algebraic Modeling System)
- [COBRApy](https://opencobra.github.io/cobrapy/) — Python package for constraint-based metabolic modeling
