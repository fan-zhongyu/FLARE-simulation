# FLARE-simulation

Part of the simulation code of the FLARE manuscript.

| Driver | Simulation |
|---|---|
| `01_pure.R` | simulated factors, FLARE with and without the graph |
| `02_semi_real.R` | semi-real data, FLARE and comparison methods |
| `03_null_gwas.R` | null GWAS, number of factors returned |
| `04_false_discovery.R` | false-discovery designs |
| `05_gfa.R`, `06_flash.R` | GFA and FLASH on the same data sets |

`R/` holds the scoring functions and the adapters of the comparison methods, `slurm/` the batch scripts.
The drivers call the `FLARE` R package. The simulated data sets are not included; their location is set
through the environment variables listed at the top of `R/common.R`.
