# PE150 sequencing analysis code

Code and example files for TAM sequence characterization and directed protein evolution library mutation analysis.

## Workflows

| Folder | Purpose | Main processing | Guide |
| --- | --- | --- | --- |
| `01_TAM_analysis/` | Count TAM sequences per sample | Ultraplex demultiplexing followed by target-sequence search and counting | [README](01_TAM_analysis/README.md) |
| `02_Protein_library_mutation_analysis/` | Track amino-acid mutations across library rounds and replicates | Dual-barcode demultiplexing, fastp QC, minibwa mapping and mutation counting | [README](02_Protein_library_mutation_analysis/README.md) |

Each workflow README lists input formats, software, commands, outputs and calculation rules. Example barcode tables, sample metadata and reference sequences must be matched to the dataset being analyzed.

The protein analysis folder contains minibwa visualizations and a separate local workflow for log2 fold-change heatmaps and fragment-level coverage. Their input tables and coverage rules differ; the READMEs describe those calculation paths separately.

The folders contain scripts, example metadata and reference sequences. Raw FASTQ reads and large intermediate files are stored separately. 
