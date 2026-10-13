# Antibody Humanization Pipeline

HumAbY v1.0 is a computational platform for antibody humanization through complementarity-determining region (CDR) grafting onto human germline-derived frameworks, followed by the application of predefined rules for selective residue modifications in framework regions. The platform provides a graphical user interface, command-line execution, batch processing, and structured report generation.

## Contents

- [Overview](#overview)

- [Installation](#installation)

- [Project structure](#project-structure)

- [Databases](#databases)

- [Scientific background and implementation](#scientific-background-and-implementation)

- [Input format](#input-format)

- [Usage](#usage)

- [Pipeline modes and framework backmutation rules](#pipeline-modes-and-framework-backmutation-rules)

- [Output files](#output-files)

- [Reports and comparison](#reports-and-comparison)

- [Troubleshooting](#troubleshooting)

- [License](#license)

- [References](#references)

## Overview

HumAbY v1.0 processes the amino acid sequences of antibody heavy-chain (VH) and light-chain (VL) variable domains. The computational workflow uses ANARCII to identify and delineate complementarity-determining regions (CDRs) and framework regions (FRs) according to the Kabat numbering scheme. Human germline-derived frameworks are then selected through BLAST+ sequence alignments, followed by grafting of the parental CDRs onto the selected human frameworks and the application of predefined rules for selective framework residue modifications, according to the specified processing parameters. Finally, the platform generates humanized antibody candidate sequences, framework alignment and selection metrics, detailed execution logs, and optional structured reports documenting the backmutations performed.

## Installation

### Prerequisites

HumAbY v1.0 requires Python 3.10 or later and the following software dependencies:

- **Python 3.10+** — Programming environment with virtual environment support.

- **ANARCII** — Antibody variable-domain numbering and CDR/framework identification according to the Kabat numbering scheme.

- **PyTorch** — Computational dependency required by ANARCII.

- **Biopython** — Biological sequence processing and analysis.

- **NCBI BLAST+** — Sequence alignment and reference database management (blastp, makeblastdb, and blastdbcmd).

- **Dash** — Web-based graphical user interface.

### Python environment setup

All commands should be executed from the project's root directory.

**1. Create and activate a virtual environment**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

If a virtual environment already exists, activate it using:

```bash
source .venv/bin/activate
```

**2. Install Python dependencies**

```bash
python -m pip install biopython anarcii
python -m pip install -r requirements-dashboard.txt
```

ANARCII requires PyTorch and its associated dependencies. For installation instructions and platform-specific requirements, consult the [official ANARCII repository](https://github.com/oxpig/ANARCII).

### NCBI BLAST+ installation

NCBI BLAST+ must be installed separately and its executables must be available in the system PATH.

**Ubuntu/Debian:**

```bash
sudo apt-get update
sudo apt-get install ncbi-blast+
```

For other operating systems, download the appropriate distribution from the [official NCBI BLAST+ repository](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/).

### Installation verification

After installation, verify that the required Python packages and BLAST+ executables are accessible.

**Python dependencies:**

```bash
python -c "import Bio, anarcii, torch, dash"
```

**ANARCII:**

```bash
anarcii --help
```

**NCBI BLAST+:**

```bash
blastp -version
makeblastdb -version
blastdbcmd -version
```

Successful execution of these commands confirms that the principal software dependencies are accessible in the configured environment.

### Reference computational environment

The following software versions are documented for the reference computational environment:

| **Software** | **Version** |
| --- | --- |
| Python | 3.10.12 |
| Biopython | 1.88 |
| ANARCII | 2.0.7 |
| PyTorch | 2.14.1 |
| Dash | 4.4.1 |

These versions describe the reference environment. For reproducibility of the published analyses, the software versions, reference databases, and execution parameters should correspond to those documented in the associated manuscript and supplementary materials.

## Project structure

The HumAbY v1.0 repository is organized into modules for antibody humanization, workflow execution, reference database management, and result processing.

```text
.
├── dash_app_integrado.py    # Dash-based graphical user interface
├── run_pipeline.py         # Command-line execution entry point
├── Makefile                # Shortcuts for pipeline execution and the dashboard
├── core/                   # Core antibody humanization and backmutation algorithms
├── integration/            # Workflow orchestration, batch processing, and reporting
├── data/
│   ├── inputs/             # Input files containing parental VH/VL sequences
│   └── databases/          # Human and murine reference sequence databases
├── results/
│   ├── analysis/           # Results from graphical interface and batch executions
│   └── manual/             # Results from command-line executions
└── tools/
    └── databases/          # Reference database construction utilities
```

The computational architecture separates the core antibody humanization algorithms from workflow execution and result management. The core/ directory contains the scientific procedures responsible for candidate generation and framework residue modifications, whereas the integration/ directory coordinates individual and batch executions, manages working directories, and organizes computational outputs.

Each execution is performed in a dedicated working directory, preserving the input sequences, processing parameters, execution logs, and generated results. This organization supports traceability, systematic comparison of humanized candidates, and reproducibility of computational analyses.

## Databases

### Human germline reference library (IMGT/GENE-DB)

HumAbY v1.0 uses a human germline reference library constructed from functional immunoglobulin gene alleles obtained from IMGT/GENE-DB. The library was designed to provide human germline-derived reference sequences for selecting acceptor frameworks during antibody humanization.

The reference library was generated by combining human immunoglobulin variable (V), diversity (D), and joining (J) gene segments according to the immunoglobulin chain type.

For heavy chains (VH), functional IGHV, IGHD, and IGHJ gene segments were combined to generate theoretical V(D)J sequences. For light chains (VL), IGKV and IGKJ gene segments were combined to generate kappa VJ sequences, while IGLV and IGLJ gene segments were combined to generate lambda VJ sequences.

The reference library was constructed using the following gene-segment combinations:

| **Chain** | **Gene segments** | **Combination strategy** |
| --- | --- | --- |
| Heavy chain (VH) | IGHV, IGHD, IGHJ | V × D × J |
| Kappa light chain (Vκ) | IGKV, IGKJ | V × J |
| Lambda light chain (Vλ) | IGLV, IGLJ | V × J |

The resulting human germline reference library comprised:

| **Reference sequences** | **Number** |
| --- | --- |
| Heavy-chain V(D)J combinations (VH) | 241,920 |
| Light-chain VJ combinations (VL) | 1,851 |
| **Total germline-derived combinations** | **243,771** |

These sequences represent computationally generated combinations of human germline gene segments and should not be interpreted as naturally observed rearranged antibody sequences.

The resulting reference library was used for BLASTP-based sequence comparison and selection of human acceptor frameworks for subsequent CDR grafting and selective framework residue modification.

### Murine germline reference database

HumAbY v1.0 uses a murine immunoglobulin germline V-gene reference database to identify framework residues potentially associated with somatic hypermutation (SHM) in parental antibodies.

Parental heavy-chain (VH) and light-chain (VL) variable-region amino acid sequences are independently aligned against murine germline V-gene references using NCBI BLASTP. The highest-ranking germline reference is retained for each chain, and framework positions at which parental residues differ from their corresponding inferred germline residues are classified as putative SHM-derived sites.

These positions are subsequently considered during selective framework backmutation. When the human acceptor framework contains a residue different from the parental antibody at a position selected for modification, the parental residue is restored according to the predefined backmutation rules.

#### NCBI murine germline database

The backmutation workflow uses a murine V-gene reference database originating from the NCBI IgBLAST distribution.

In the current documented repository structure, the database is located at:

```text
data/databases/mouse_ncbi/mouse_V
```

The database contains **512 sequences and 50,181 amino acid residues**. According to the repository's provenance records, it originates from the NCBI IgBLAST 1.22.0 distribution and was built on August 24, 2022.

The original public implementation also references the database at:

```text
back/ncbi-igblast/internal_data/mouse/mouse_V
```

Although the reference database originates from the NCBI IgBLAST distribution, sequence alignment in the documented backmutation implementation is performed using the BLAST+ blastp executable.

#### Database provenance and verification

Information regarding the database source, location, sequence count, and SHA-256 checksum is recorded in:

```text
data/databases/catalogs.json
```

The installed database can be inspected using:

```bash
blastdbcmd -db data/databases/mouse_ncbi/mouse_V -info
```

This command reports the indexed sequence count, total amino acid residues, and database build information.

#### Reproducibility considerations

The murine germline reference database influences germline assignment and the identification of putative SHM-derived framework positions, which may affect subsequent backmutation decisions.

For reproducibility, the reference database version, alignment settings, and corresponding software implementation should be documented alongside each computational analysis.

Differences between parental and inferred germline sequences are interpreted as putative SHM-associated positions rather than definitive evidence of somatic hypermutation.

### Human germline reference database construction

HumAbY v1.0 provides a computational utility for constructing human germline-derived reference databases from immunoglobulin gene-segment sequences obtained from IMGT/GENE-DB.

The construction procedure processes human immunoglobulin variable (V), diversity (D), and joining (J) gene segments to generate theoretical heavy-chain V(D)J and light-chain VJ combinations. These sequences serve as references for identifying human acceptor frameworks during antibody humanization.

The repository includes a database construction script that can be executed using the preserved source files.

**1. Activate the Python environment**

```bash
source .venv/bin/activate
```

**2. Construct the human germline reference database**

```bash
python tools/databases/build_germline_database.py \
  data/databases/human_imgt/sources \
  /tmp/humanizer-human-imgt-rebuild
```

The script processes the supplied nucleotide sequences, translates gene segments into amino acid sequences, applies sequence-quality filters, and reconstructs heavy-chain and light-chain reference sequences.

The procedure also performs light-chain FR4 consistency checks, generates FASTA files and database statistics, and builds searchable BLAST indices using makeblastdb.

**3. Verify the reconstructed database**

```bash
blastdbcmd \
  -db /tmp/humanizer-human-imgt-rebuild/human_germline_frameworks \
  -info
```

This command reports the number of indexed sequences, the total number of amino acid residues, and database build information.

#### Reproducibility considerations

The human germline reference library described in the associated HumAbY v1.0 manuscript comprised 241,920 heavy-chain V(D)J combinations and 1,851 light-chain VJ combinations, totaling **243,771 theoretical germline-derived combinations**.

Reproducing this exact library requires the original IMGT/GENE-DB source files, accessed in September 2025, together with the corresponding database construction procedures and filtering criteria.

The number of reconstructed sequences may vary depending on the source dataset, gene-segment composition, and sequence-processing filters.

#### Updating the reference database

To construct a database from a different IMGT/GENE-DB release, obtain the corresponding human immunoglobulin gene-segment FASTA files from [IMGT/GENE-DB](https://www.imgt.org/genedb/).

The required source files are:

- Heavy chain: IGHV.fasta, IGHD.fasta, and IGHJ.fasta.

- Kappa light chain: IGKV.fasta and IGKJ.fasta.

- Lambda light chain: IGLV.fasta and IGLJ.fasta.

New databases should be constructed in separate output directories to preserve existing references and their associated computational results.

After successful construction and verification, the generated database can be selected using the --database command-line argument or the HUMANIZER_DB_PATH environment variable.

## Scientific background and implementation

HumAbY v1.0 implements a rule-based computational workflow for antibody humanization through complementarity-determining region (CDR) grafting onto human germline-derived acceptor frameworks, followed by selective framework residue modifications.

The workflow processes parental heavy-chain (VH) and light-chain (VL) variable-domain amino acid sequences independently, preserving their original association for subsequent generation of humanized VH–VL combinations.

### CDR annotation and framework identification

Parental and human germline-derived variable-region sequences are numbered using ANARCII according to the Kabat numbering scheme.

The variable domains are divided into framework regions (FR1–FR4) and complementarity-determining regions (CDR1–CDR3), following the sequence organization:

```text
FR1 — CDR1 — FR2 — CDR2 — FR3 — CDR3 — FR4
```

CDR boundaries are defined as follows:

| **Region** | **Heavy chain (VH)** | **Light chain (VL)** |
| --- | --- | --- |
| CDR1 | Kabat H31–H35 | Kabat L24–L34 |
| CDR2 | Kabat H50–H65 | Kabat L50–L56 |
| CDR3 | Defined by conserved flanking motifs | Kabat L89–L97 |

Insertion-coded positions are included within their corresponding CDR boundaries. For VH, CDR3 is identified between the conserved cysteine-containing motif at the end of FR3 and the WG motif marking the beginning of FR4.

All remaining residues within the numbered variable domains are assigned to their corresponding framework regions.

### Human acceptor framework selection

Parental VH and VL sequences are independently compared against human germline-derived reference libraries using BLASTP.

The reference libraries contain computationally reconstructed heavy-chain V(D)J and light-chain VJ combinations derived from human immunoglobulin germline gene segments.

BLASTP searches are performed using an E-value threshold of 1 × 10⁻¹⁰, retaining up to 50 target sequences per parental chain.

Candidate human acceptor frameworks are identified from the ranked alignments. Framework sequences are compared to eliminate redundant candidates, and up to five highest-ranking nonredundant frameworks are retained for each chain.

### CDR grafting and sequence validation

For each selected human acceptor framework, parental CDR1, CDR2, and CDR3 sequences are grafted onto the corresponding human framework regions (FR1–FR4).

The resulting reconstructed variable domains preserve the parental CDR sequences while incorporating human germline-derived framework residues.

Following CDR grafting, candidate sequences are reprocessed using ANARCII to verify successful variable-domain numbering and appropriate chain classification.

Candidates are retained only when the expected framework–CDR organization is preserved, without unintended alterations in the grafted CDR sequences or their junctions with adjacent framework regions.

### Selective framework modifications

Following CDR grafting, predefined sequence- and position-based rules are applied to identify framework residues selected for modification.

These rules include the evaluation of putative somatic hypermutation-derived positions, predefined Vernier-zone residues, selected cysteine and proline substitutions, and heavy-chain Kabat position H71.

When the human acceptor framework differs from the parental antibody at a position selected for backmutation, the corresponding parental residue is restored according to the implemented modification criteria.

Following framework modification, FR4 sequences are normalized to contain the WG-X-G motif in VH and FG-X-G motif in VL, preserving the X residue derived from the corresponding human J-gene sequence.

FR4 normalization is treated as a sequence-standardization procedure and is not classified as a parental framework backmutation.

The complete residue-level modification criteria are described in the associated manuscript and supplementary materials.

### Generation of humanized VH–VL candidates

Following CDR grafting, framework modifications, and sequence validation, the retained humanized VH and VL sequences are combinatorially paired within each parental antibody.

For *nH* heavy-chain candidates and *nL* light-chain candidates, the workflow generates *nH × nL* humanized VH–VL combinations.

Pairing between sequences originating from different parental antibodies is not performed.

In the associated HumAbY v1.0 study, this approach generated 60 humanized VH–VL combinations from three parental anti-IL-7Rα antibodies, of which 13 were selected for experimental characterization.

## Input format

HumAbY v1.0 accepts amino acid sequences corresponding to the variable domains of parental antibody heavy (VH) and light (VL) chains.

For individual antibody humanization, the input must be provided as a plain-text file (.aa) containing exactly two lines, without FASTA headers:

- **First line:** Heavy-chain variable-domain (VH) amino acid sequence.

- **Second line:** Light-chain variable-domain (VL) amino acid sequence.

Both sequences must correspond to the same parental antibody and preserve their original VH–VL association.

### Example input: 4A10 antibody

The repository includes an example input file, data/inputs/4A10.aa, containing the parental VH and VL amino acid sequences of the anti-IL-7Rα antibody 4A10.

```text
QVQLQQPGAELVMPGASVKLSCKASGYTFTSYWMHWVKQRPGEGLEWIGEIDPSDSYTNDNQKFKGKATLTVDKSSSTAYMQLSSLTSEDSAVYYCARRLYSNSYYYAMDYWGQGTSVTVSS
DIQMTQSPSSLSASLGGKVTITCKASQDIKKYIAWYQHKPGKGPRLLIHYTSTLQPGIPSRFSGSGSGRDYSFSISNLEPVDIATYYCLQYDNLLTFGAGTKLELK
```

These sequences provide an example input for the antibody humanization workflow.

### Sequence requirements

Input sequences must satisfy the following criteria:

- Contain complete VH and VL variable-domain amino acid sequences.

- Include at least 50 amino acid residues per chain.

- Use the standard amino acid alphabet: ACDEFGHIKLMNPQRSTVWY.

- Maintain the required sequence order: VH followed by VL.

- Exclude FASTA headers and additional annotation lines.

### Input validation

The graphical interface validates and normalizes input sequences before processing. The command-line workflow also performs sequence-cleaning and validation procedures, including minimum sequence-length and amino acid diversity checks.

Valid input sequences are subsequently processed for variable-domain numbering, CDR and framework identification, and human acceptor framework selection.

The 4A10 reference input provides a reproducible example for running the computational workflow without requiring disclosure of proprietary antibody sequences.

## Usage

### Graphical user interface (GUI)

HumAbY v1.0 provides a web-based graphical user interface developed with Dash, allowing users to perform antibody humanization analyses without directly interacting with the command-line workflow.

The interface supports individual antibody processing and batch execution, providing access to computational settings, execution status, generated humanized candidates, and structured analysis reports.

#### Launching the interface

From the project's root directory, activate the Python virtual environment and start the Dash application:

```bash
source .venv/bin/activate
python dash_app_integrado.py
```

Once the application is running, open the following address in a web browser:

[http://127.0.0.1:8055](http://127.0.0.1:8055/)

The interface allows users to:

- Submit parental VH and VL amino acid sequences for humanization.

- Process multiple antibody sequence pairs in batch mode.

- Select available processing options and computational settings.

- Monitor execution status and inspect processing logs.

- Access generated humanized VH and VL candidate sequences.

- Review computational outputs and structured backmutation reports, when generated.

To stop the application, return to the terminal and press Ctrl+C.

#### Environment configuration

The graphical interface can be configured using the following environment variables:

| **Environment variable** | **Default value** | **Description** |
| --- | --- | --- |
| DASH_HOST | 127.0.0.1 | Network address used to host the Dash application |
| DASH_PORT | 8055 | Port used to access the graphical interface |
| HUMANIZER_RUNS_DIR | results/analysis | Directory for storing execution results |
| HUMANIZER_DB_PATH | Human IMGT database prefix | Path to the human germline-derived reference database used for framework selection |

By default, the application runs locally and can be accessed through the specified browser address.

Each analysis is stored in a dedicated execution directory, allowing users to retrieve input sequences, processing settings, execution logs, and generated computational results.

### Command-line interface (CLI)

HumAbY v1.0 provides a command-line interface (CLI) for executing antibody humanization workflows, allowing users to configure processing modes, select reference databases, apply framework modifications, and specify output directories.

The CLI supports individual antibody analyses and can be integrated into automated computational workflows.

#### Running the pipeline

All commands should be executed from the project's root directory with the Python virtual environment activated.

**1. Display available commands and options**

```bash
python run_pipeline.py --help
```

This command displays the available execution parameters and their corresponding descriptions.

**2. Perform CDR grafting**

```bash
python run_pipeline.py data/inputs/4A10.aa --mode graft
```

This command processes the parental VH and VL sequences of the 4A10 antibody, selects human germline-derived acceptor frameworks, and generates humanized candidates through CDR grafting.

**3. Execute the workflow with framework backmutation**

```bash
python run_pipeline.py data/inputs/4A10.aa \
  --mode complete \
  --optimization 1 \
  --backmutation \
  --output results/manual/4A10_level1
```

This example executes the complete processing mode at optimization level 1, followed by the independent backmutation procedure. Generated results are stored in the specified output directory.

The selected optimization level determines which processing operations are performed before the optional backmutation stage.

#### Command-line arguments

| **Argument** | **Description** |
| --- | --- |
| input | Input file containing parental VH and VL amino acid sequences in two-line format |
| --mode | Selects the pipeline processing mode; default: complete |
| --optimization | Specifies the processing level (1–4) when using complete mode; default: 4 |
| --backmutation | Enables the independent framework backmutation procedure following complete or graft processing |
| --database | Specifies the BLAST database prefix for the human germline-derived framework reference library |
| --output | Specifies a new directory for storing execution results; existing output directories are rejected |

Additional information about available processing modes and optimization levels is provided in the [Pipeline modes and framework backmutation rules](#pipeline-modes-and-framework-backmutation-rules) section.

#### Makefile shortcuts

The repository includes a Makefile with predefined commands for executing the pipeline and launching the graphical interface.

**Run the 4A10 example with backmutation:**

```bash
make NAME=4A10 ARGS="--mode complete --optimization 1 --backmutation"
```

**Launch the Dash graphical interface:**

```bash
make dashboard
```

These shortcuts use the Python environment configured for the project.

#### Direct scientific execution

The core antibody humanization workflow can also be accessed through:

```text
core/humanization/humanize.py
```

Available execution parameters can be inspected using the script's --help option.

Direct execution of the scientific core requires appropriate configuration of the working directory and reference database paths. For routine analyses, the run_pipeline.py interface is recommended because it coordinates workflow execution and organizes generated results.

### Batch processing

HumAbY v1.0 supports batch processing of multiple parental antibody VH–VL sequence pairs, allowing users to perform humanization analyses using consistent computational settings.

Batch processing is available through the graphical user interface and can also be performed from the command line.

#### Batch input through the graphical interface

The graphical interface accepts tabular input containing three fields:

- **name** — Unique identifier for each parental antibody.

- **VH** — Heavy-chain variable-domain amino acid sequence.

- **VL** — Light-chain variable-domain amino acid sequence.

Input fields may be separated by commas, semicolons, or tabs. Column headers are optional, and sample identifiers must be nonempty and unique.

An example of comma-separated batch input is provided below:

```csv
name,VH,VL
sample_1,<VH amino acid sequence>,<VL amino acid sequence>
sample_2,<VH amino acid sequence>,<VL amino acid sequence>
```

Replace the placeholders with the corresponding complete VH and VL amino acid sequences. Each row must contain sequences derived from the same parental antibody.

All submitted samples are processed sequentially using the selected computational settings. Each analysis receives an individual execution directory containing its corresponding results, logs, and processing status.

#### Batch execution through the command line

Multiple antibody input files can also be processed sequentially using a shell loop.

For example, to execute the workflow for all .aa files located in data/inputs/:

```bash
for input_file in data/inputs/*.aa; do
  python run_pipeline.py "$input_file" \
    --mode complete \
    --optimization 1
done
```

Each input file must follow the two-line VH/VL format described in the [Input format](#input-format) section.

This approach allows multiple parental antibodies to be processed using the same execution parameters without submitting each file individually.

#### Reproducibility considerations

Each antibody is processed as an independent analysis. The original VH–VL association is preserved within each parental antibody, and candidate generation is performed separately for each input pair.

Batch processing does not introduce cross-pairing between heavy and light chains originating from different parental antibodies.

Generated results should be interpreted according to the selected processing mode, optimization level, and reference database configuration.

## Pipeline modes and framework backmutation rules

HumAbY v1.0 implements a rule-based antibody humanization workflow based on CDR grafting onto human germline-derived acceptor frameworks, followed by selective framework residue modifications.

### Pipeline modes

| **Mode** | **Description** |
| --- | --- |
| structures | Identifies and extracts CDR1–CDR3 and FR1–FR4 regions from parental variable-domain sequences using ANARCII and Kabat numbering. |
| scores | Performs BLASTP alignments against human germline-derived reference sequences and reports sequence identity, E-value, and bit score. |
| graft | Generates humanized variable-domain candidates by grafting parental CDR sequences onto selected human acceptor frameworks. |
| graft_backmutation | Performs CDR grafting followed by selective framework backmutation according to predefined residue-level rules. |
| complete | Executes the integrated humanization workflow, including framework selection, CDR grafting, sequence validation, and selective framework modifications. |

### Selective framework backmutation rules

Following CDR grafting, HumAbY v1.0 applies predefined sequence- and position-based rules to selectively restore parental framework residues in humanized antibody candidates.

The framework modification strategy includes the following criteria:

| **Rule** | **Description** |
| --- | --- |
| **Putative SHM-derived positions** | Framework positions at which parental residues differ from their corresponding inferred murine germline residues are identified. When the human acceptor residue differs from the parental residue at a selected position, the parental residue is restored. |
| **Cysteine and proline safeguards** | Cysteine or proline residues introduced by the human acceptor framework at positions containing different parental residues are replaced with the corresponding parental residues. |
| **Heavy-chain Kabat position H71** | The parental residue at VH Kabat position H71 is restored whenever it differs from the corresponding human acceptor residue, considering its established structural relevance to CDR-H2 conformation. |

### FR4 sequence normalization

Following selective framework backmutation, the FR4 regions are normalized to preserve the following conserved motifs:

| **Chain** | **FR4 motif** |
| --- | --- |
| Heavy chain (VH) | WG-X-G |
| Light chain (VL) | FG-X-G |

The X residue is retained from the corresponding human J-gene-derived sequence.

FR4 normalization is treated as a sequence-standardization step and is not classified or counted as a parental framework backmutation.

### Framework modification considerations

Framework modifications are restricted to positions selected according to the predefined residue-level rules. The parental amino acid is restored only when it differs from the corresponding human acceptor residue.

The strategy aims to preserve framework residues potentially relevant to antigen recognition while retaining human germline-derived framework sequences.

The complete residue-level rules, including predefined Vernier-zone positions, triggering conditions, and resulting modifications, are described in **Supplementary Table S1** of the associated manuscript.

The contribution of individual framework modification rules was not independently evaluated in the study. Therefore, these modifications should not be interpreted as experimentally validated improvements in antigen-binding affinity, stability, or developability.

## Output files

HumAbY v1.0 organizes computational outputs into dedicated directories for each execution, allowing users to access humanized antibody candidate sequences, reference framework information, alignment results, and framework backmutation reports.

By default, analyses performed through the graphical interface are stored in results/analysis/, whereas command-line executions are stored in results/manual/. Alternatively, users can specify a custom output directory using the --output argument.

### Execution directory structure

Each execution generates a dedicated directory containing the following components:

```text
<run_directory>/
├── input/                  # Parental VH/VL input sequences
├── output/                 # Generated candidates and analysis results
├── logs/
│   └── pipeline.log        # Computational execution log
└── metadata.json           # Execution metadata and processing settings
```

### Generated output files

Depending on the selected processing mode and backmutation configuration, the output/ directory may contain the following files:

| **Output file** | **Description** |
| --- | --- |
| humanized_vh.fasta | Humanized VH candidate sequences, including candidate identifiers, human framework references, and BLASTP sequence identity information. |
| humanized_vl.fasta | Humanized VL candidate sequences, including candidate identifiers, human framework references, and BLASTP sequence identity information. |
| humanization_summary.txt | Summary of generated candidates, selected human reference sequences, alignment identity values, and variable-domain sequence lengths. |
| FRfileVH, FRfileVL | Framework sequences extracted from the parental VH and VL domains for subsequent backmutation analysis. |
| hGermVH, hGermVL | Selected human germline-derived framework reference sequences arranged according to candidate order. |
| VH.xml, VL.xml | BLASTP alignment results comparing parental VH and VL sequences with murine germline V-gene references, used during the backmutation procedure. |
| &lt;sample&gt;.backmutated | Raw output generated by the independent framework backmutation procedure. |
| backmutation_organized.txt | Human-readable summary of framework residue modifications and associated candidate information. |
| backmutation_organized.json | Structured report containing framework modification annotations and corresponding candidate data. |
| backmutation_rules_observed.txt | Human-readable records of residue-level modification rules identified through execution logs or verified reconstruction of the generated results. |
| backmutation_rules_observed.json | Structured representation of the observed residue-level modification rules and their supporting evidence. |
| execution_conditions.json | Records input sequences, computational settings, code and database hashes, and software environment information. |

### Output organization and traceability

The generated files depend on the selected execution mode and whether independent backmutation is enabled. Modes restricted to sequence annotation or alignment analysis may generate fewer output files than complete humanization runs.

Execution logs and metadata record the computational settings associated with each analysis. When batch processing is performed, metadata also tracks the corresponding individual executions.

The structured backmutation reports facilitate the inspection of framework residue modifications and provide supporting information for comparing humanized candidates.

### Reproducibility considerations

To support reproducibility, the original input sequences, reference databases, processing parameters, software versions, and generated outputs should be preserved together.

The recorded execution conditions and file hashes facilitate the identification of computational environments and reference datasets associated with individual analyses.

These records support the traceability of computational results but do not independently establish the functional or structural validity of the generated antibody candidates.

## Reports and comparison

HumAbY v1.0 provides structured reports to support the inspection, documentation, and comparison of humanized antibody candidates and their framework residue modifications.

The reporting system organizes information generated during CDR grafting and selective framework backmutation, facilitating the identification of modified residues, their corresponding reference sequences, and the computational criteria associated with each modification.

### Report generation

For an existing completed backmutation run, organized reports can be generated using:

```bash
python -m integration.reports.service results/manual/4A10_level1/output
```

The command processes the available candidate sequences, reference alignments, and backmutation results to generate structured reports within the corresponding output directory.

### Report formats

The reporting system generates two complementary formats:

| **Format** | **Description** |
| --- | --- |
| **TXT** | Human-readable summaries of humanized antibody candidates, framework residue modifications, and associated computational information. |
| **JSON** | Structured records containing candidate identifiers, reference sequences, residue-level annotations, supporting rule evidence, validation status, and comparison information. |

Reports are stored alongside the corresponding FASTA sequences, BLASTP alignment files, and raw backmutation outputs.

### Framework modification traceability

The reports document framework modifications associated with the predefined criteria described in the HumAbY v1.0 methodology:

- **Putative SHM-derived positions:** Framework residues selected based on differences between parental sequences and inferred murine germline references.

- **Cysteine and proline safeguards:** Restoration of parental residues at selected framework positions where cysteine or proline residues were introduced by the human acceptor framework.

- **VH Kabat position H71:** Selective restoration of the parental residue at this structurally relevant heavy-chain framework position.

FR4 motif normalization (WG-X-G for VH and FG-X-G for VL) is documented separately from parental framework backmutation.

Rule evidence may originate from explicit execution logs or from a verified computational replay that reproduces the saved candidate sequences.

Modification counts represent distinct residue positions associated with each rule. Because multiple rules may apply to the same position, rule-specific counts are not necessarily additive.

### Result identification and reproducibility

The reporting system uses SHA-256 identifiers to support computational traceability:

| **Identifier** | **Purpose** |
| --- | --- |
| result_sha256 | Identifies the structured computational result independently of file paths and timestamps. |
| execution_sha256 | Identifies recorded execution conditions, including software code and reference database hashes. |

For historical analyses without recorded execution conditions, the reporting system assigns:

```text
execution_status: not_recorded
```

These identifiers facilitate comparisons between computational runs and help users determine whether results were generated under equivalent recorded conditions.

### Interpretation and limitations

The reports document computationally identified framework modifications and the rules associated with their application. They do not independently establish the structural or functional contribution of individual backmutations.

Experimental characterization remains necessary to determine whether the generated antibody candidates preserve antigen recognition and exhibit suitable functional properties.

## Troubleshooting

The following table summarizes common issues that may occur during the installation, execution, and analysis of results generated by HumAbY v1.0.

| **Component** | **Recommended action** |
| --- | --- |
| **Python environment** | Activate the project's .venv environment and verify that all required Python packages and external executables are correctly installed, as described in  [Installation](#installation). |
| **BLAST+ installation** | Confirm that blastp, makeblastdb, and blastdbcmd are installed and accessible through the system PATH. |
| **Human reference database** | Verify the human germline-derived BLAST database path and use its prefix without extensions such as .fasta, .pin, or .phr. Check the index using blastdbcmd -info. |
| **Murine reference database** | Confirm that the murine germline V-gene database required for framework backmutation is accessible and properly indexed. |
| **Input sequences** | Provide complete VH and VL amino acid sequences in the required order, using the standard amino acid alphabet and the documented input format. |
| **Candidate generation** | Inspect BLASTP template matches, CDR and framework annotations, and processing messages in logs/pipeline.log. |
| **Backmutation processing** | Verify the input and candidate sequence lengths, framework reference files, and corresponding murine BLASTP XML alignments. |
| **Structured reports** | Ensure that FASTA sequences, framework and germline reference files, XML alignments, and raw backmutation outputs originate from the same execution. |
| **Output directory** | Specify a new output directory when necessary, as existing destinations are not overwritten. |
| **Graphical interface** | Confirm that the configured interface port is available. Stop any previous application instance or modify the DASH_PORT environment variable. |
| **Result reproducibility** | Verify the reference database version, software dependencies, processing parameters, and recorded execution conditions when comparing results from different runs. |

If an error persists, inspect the corresponding execution log (`logs/pipeline.log`).

## License

HumAbY v1.0 is distributed under the **Apache License, Version 2.0**.

The official English license text is available in  [LICENSE-APACHE-2.0.txt](core/humanization/LICENSE-APACHE-2.0.txt).

An unofficial Portuguese translation is provided in  [LICENSE](core/humanization/LICENSE) for informational purposes only. The official English version takes precedence in the event of any discrepancy.

## References

### Core Antibody Humanization Literature

- Jones, P.T., Dear, P.H., Foote, J., Neuberger, M.S. & Winter, G. Replacing the complementarity-determining regions in a human antibody with those from a mouse. *Nature*  **321**, 522–525 (1986). [DOI: 10.1038/321522a0](https://doi.org/10.1038/321522a0)

- Verhoeyen, M., Milstein, C. & Winter, G. Reshaping human antibodies: grafting an antilysozyme activity. *Science*  **239**, 1534–1536 (1988). [DOI: 10.1126/science.2451287](https://doi.org/10.1126/science.2451287)

- Winter, G. & Milstein, C. Man-made antibodies. *Nature*  **349**, 293–299 (1991). [DOI: 10.1038/349293a0](https://doi.org/10.1038/349293a0)

### Antibody Structure and Numbering Systems

- Kabat, E.A., Wu, T.T., Perry, H.M., Gottesman, K.S. & Foeller, C. *Sequences of Proteins of Immunological Interest*, 5th edn. (U.S. Department of Health and Human Services, 1991).

- Kabat, E.A. & Wu, T.T. Identical V region amino acid sequences and segments of sequences in antibodies of different specificities. *J. Immunol.*  **147**, 1709–1719 (1991). [PMID: 1880416](https://pubmed.ncbi.nlm.nih.gov/1880416/)

**ANARCI Numbering Tool (predecessor):**

- Dunbar, J. et al. ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics*  **32**, 298–300 (2016). [DOI: 10.1093/bioinformatics/btv552](https://doi.org/10.1093/bioinformatics/btv552)

### V(D)J Recombination and Germline Databases

**V(D)J Recombination Mechanisms:**

- Tonegawa, S. Somatic generation of antibody diversity. *Nature*  **302**, 575–581 (1983). [DOI: 10.1038/302575a0](https://doi.org/10.1038/302575a0)

- Schatz, D.G., Oettinger, M.A. & Baltimore, D. The V(D)J recombination activating gene, RAG-1. *Cell*  **59**, 1035–1048 (1989). [DOI: 10.1016/0092-8674(89)90760-5](<https://doi.org/10.1016/0092-8674(89)90760-5>)

**IMGT Database System:**

- Lefranc, M.-P. et al. IMGT, the international ImMunoGeneTics information system. *Nucleic Acids Res.*  **33**, D593–D597 (2005). [DOI: 10.1093/nar/gki065](https://doi.org/10.1093/nar/gki065)

- Giudicelli, V., Chaume, D. & Lefranc, M.-P. IMGT/GENE-DB: a comprehensive database for human and mouse immunoglobulin and T cell receptor genes. *Nucleic Acids Res.*  **33**, D256–D261 (2005). [DOI: 10.1093/nar/gki010](https://doi.org/10.1093/nar/gki010)

### Computational Methods and Algorithms

**Sequence Analysis:**

- Altschul, S.F. et al. Basic local alignment search tool. *J. Mol. Biol.*  **215**, 403–410 (1990). [DOI: 10.1016/S0022-2836(05)80360-2](<https://doi.org/10.1016/S0022-2836(05)80360-2>)

- Camacho, C. et al. BLAST+: architecture and applications. *BMC Bioinformatics*  **10**, 421 (2009). [DOI: 10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)

**BioPython Framework:**

- Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*  **25**, 1422–1423 (2009). [DOI: 10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

### Framework Engineering and Optimization

**Framework Region Analysis:**

- Foote, J. & Winter, G. Antibody framework residues affecting the conformation of the hypervariable loops. *J. Mol. Biol.*  **224**, 487–499 (1992). [DOI: 10.1016/0022-2836(92)91010-M](<https://doi.org/10.1016/0022-2836(92)91010-M>)

- Tramontano, A., Chothia, C. & Lesk, A.M. Framework residue 71 is a major determinant of the position and conformation of the second hypervariable region in the VH domains of immunoglobulins. *J. Mol. Biol.*  **215**, 175–182 (1990). [DOI: 10.1016/S0022-2836(05)80102-0](<https://doi.org/10.1016/S0022-2836(05)80102-0>)

**Structural Analysis:**

- Chothia, C. & Lesk, A.M. Canonical structures for the hypervariable regions of immunoglobulins. *J. Mol. Biol.*  **196**, 901–917 (1987). [DOI: 10.1016/0022-2836(87)90412-8](<https://doi.org/10.1016/0022-2836(87)90412-8>)

- Al-Lazikani, B., Lesk, A.M. & Chothia, C. Standard conformations for the canonical structures of immunoglobulins. *J. Mol. Biol.*  **273**, 927–948 (1997). [DOI: 10.1006/jmbi.1997.1354](https://doi.org/10.1006/jmbi.1997.1354)

### Clinical Applications and Therapeutic Development

**Immunogenicity and HAMA Response:**

- Schroff, R.W. et al. Human anti-murine immunoglobulin responses in patients receiving monoclonal antibody therapy. *Cancer Res.*  **45**, 879–885 (1985). [PMID: 3871353](https://pubmed.ncbi.nlm.nih.gov/3871353/)

- Hwang, W.Y.K. & Foote, J. Immunogenicity of engineered antibodies. *Methods*  **36**, 3–10 (2005). [DOI: 10.1016/j.ymeth.2005.01.001](https://doi.org/10.1016/j.ymeth.2005.01.001)

**Therapeutic Antibody Landscape:**

- Reichert, J.M. Antibodies to watch in 2021. *mAbs*  **13**, 1860476 (2021). [DOI: 10.1080/19420862.2020.1860476](https://doi.org/10.1080/19420862.2020.1860476)

- Kaplon, H. & Reichert, J.M. Antibodies to watch in 2019. *mAbs*  **11**, 219–238 (2019). [DOI: 10.1080/19420862.2018.1556465](https://doi.org/10.1080/19420862.2018.1556465)

### Germline Gene Repertoires

**Human Immunoglobulin Loci:**

- Matsuda, F. et al. The complete nucleotide sequence of the human immunoglobulin heavy chain variable region locus. *J. Exp. Med.*  **188**, 2151–2162 (1998). [DOI: 10.1084/jem.188.11.2151](https://doi.org/10.1084/jem.188.11.2151)

- Kawasaki, K. et al. One-megabase sequence analysis of the human immunoglobulin λ locus. *Genome Res.*  **7**, 250–261 (1997). [DOI: 10.1101/gr.7.3.250](https://doi.org/10.1101/gr.7.3.250)

### Computational Antibody Research

**Computational Antibody Design:**

- Ruffolo, J.A., Sulam, J. & Gray, J.J. Antibody structure prediction using interpretable deep learning. *Patterns*  **3**, 100406 (2022). [DOI: 10.1016/j.patter.2021.100406](https://doi.org/10.1016/j.patter.2021.100406)

- Kovaltsuk, A. et al. Observed Antibody Space: A Resource for Data Mining Next-Generation Sequencing of Antibody Repertoires. *J. Immunol.*  **201**, 2502–2509 (2018). [DOI: 10.4049/jimmunol.1800708](https://doi.org/10.4049/jimmunol.1800708)

### ANARCII Software and Method

- Greenshields-Watson, A. et al. ANARCII: A Generalised Language Model for Antigen Receptor Numbering. *bioRxiv* (2025). [DOI: 10.1101/2025.04.16.648720](https://doi.org/10.1101/2025.04.16.648720). See the [official ANARCII repository](https://github.com/oxpig/ANARCII) for installation and usage.

The literature provides methodological context for CDR grafting, antibody numbering, germline databases, and sequence analysis.
