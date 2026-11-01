# High-throughput sequencing analysis code

This repository organizes code for two analyses of PE150 high-throughput sequencing data: TAM sequence characterization and mutation-frequency monitoring in directed protein evolution libraries. It contains analysis scripts, example metadata and reference sequences. It does not claim to include every analysis or experiment in the associated manuscript.

## Workflows

| Analysis | Purpose | Main steps | Instructions |
| --- | --- | --- | --- |
| TAM sequence analysis | Find and count sequences downstream of a specified target | Demultiplex pooled R1 reads with Ultraplex and a barcode table; count TAM sequences in each sample | [TAM analysis README](01_TAM_analysis/README.md) |
| Protein library mutation analysis | Monitor amino-acid mutations and compare library rounds or biological replicates | Dual-barcode demultiplexing; fastp quality control; minibwa mapping and mutation counting; optional visualization | [Protein library analysis README](02_Protein_library_mutation_analysis/README.md) |

Run each workflow from its own directory, using separate writable working directories for data and results. Each README documents input formats, required software, configuration, outputs and limitations. Sample names in the provided metadata are examples; barcode sequences, reference sequences, PCR products, primers and experimental grouping must be checked and replaced where appropriate for each dataset.

## Visualization and interpretation

The protein workflow includes two distinct visualization paths. The minibwa-table path plots mutation frequencies, replicate correlations and coverage after aggregating counts. The local fragment-table path produces log2 fold-change heatmaps and amino-acid coverage figures; its coverage thresholds are applied per fragment before combining detections. These approaches use different inputs and definitions and should not be treated as interchangeable. See the [protein workflow instructions](02_Protein_library_mutation_analysis/README.md) for details, including reference-coordinate assumptions and a known limitation in the historical enrichment calculation.

## Data and reproducibility

Raw FASTQ reads and large intermediate files are not included in the two workflow directories. Sequencing accessions, figure assignments and the license remain to be confirmed before public release. The local `I325Y2/` directory is supporting historical visualization material; inspect its contents and publication suitability separately before including it in a public repository.

Selected preparation and plotting steps have been checked against example data, as described in the workflow READMEs. The complete workflows have not been rerun from raw sequencing data, and reproduction of final manuscript figures has not been established.
