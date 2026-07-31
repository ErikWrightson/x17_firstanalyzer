# X17 First Analysis Script Framework

## Overview

This repository contains the analysis utilities used for the X17 analysis. Together, these scripts automate the processing of reconstructed X17 experimental data, produce merged analysis ROOT files, and generate comprehensive PDF reports for validating candidate selection and monitoring the effects of each stage of the analysis.

The workflow is designed to efficiently process large data sets while providing standardized visual summaries of the reconstructed physics observables for X17 candidate searches and Møller scattering events.

The analysis chain consists of:

1. Running the `x17_firstanalyzer` over all reconstructed files for a given run.
2. Combining the individual analysis ROOT files into a single merged ROOT file.
3. Generating detailed PDF reports summarizing the reconstructed distributions for each analysis channel.
4. Reviewing the cut-by-cut evolution of candidate selection, detector response, and invariant-mass spectra.

---

# Workflow Overview

```text
Reconstructed PRad Data
          |
          v
+------------------------+
| Run_AnalysisBatch.py   |
| Batch Analysis Runner  |
+------------------------+
          |
          | Launches one instance of
          | ../x17_firstanalyzer
          | for every reconstructed file
          v
+------------------------+
| x17_firstanalyzer      |
| Main Analysis Program  |
+------------------------+
          |
          v
 Individual Analysis ROOT Files
(x17_Ana_<RUN>_<FILE>.root)
          |
          | Automatic merge (hadd)
          v
x17_Ana_<RUN>_All.root
          |
          v
+------------------------+
| PrintHistograms.py     |
| Analysis PDF Generator |
+------------------------+
          |
          v
+-----------------------------------------------+
| 3-Particle (e-e+e-) Candidate PDF             |
| Reduced 3-Particle (e-e+e-) Summary PDF       |
| (γγe) Candidate Summary PDF                   |
| Møller (e-e-) Candidate PDF                   |
+-----------------------------------------------+
```

---

# General Notes

* The plotting utility operates on merged ROOT files produced by `Run_AnalysisBatch.py`.
* Histogram discovery is performed automatically by recursively scanning the ROOT file, allowing the internal directory structure to evolve without requiring changes to the plotting code.
* Histogram naming conventions must remain consistent so that plots are correctly grouped into the three-particle, γγe, and Møller analysis categories.
* Optional live-charge normalization may be applied to invariant-mass spectra using a JSON file containing the measured run live charge.
* The generated PDF reports document the progression of the analysis through each event-selection cut and provide standardized diagnostics for detector performance and candidate reconstruction.

This README serves as the central reference for the X17 analysis framework.

Below is a list of the individual README sections describing each script used in the X17 script driven analysis workflow.

1. [Batch X17 First Analyzer Runner](#1-batch-x17-first-analyzer-runner)
2. [Analysis PDF Plot Generator](#2-analysis-pdf-plot-generator)

---
# 1. Batch X17 First Analyzer Runner

## Overview

`Run_AnalysisBatch.py` is a Python utility that automates the execution of the `x17_firstanalyzer` analysis over all reconstructed X17 experimental data files for a given run. The script processes files in parallel, merges the resulting ROOT files into a single output, and removes the intermediate per-file ROOT files upon successful completion.

This utility is intended to simplify large-scale X17 analyses while reducing total processing time through parallel execution.

## Workflow

Given a run number, the script performs the following steps:

1. Searches the input data directory for all reconstructed ROOT files corresponding to the specified run.
2. Creates the output directory if it does not already exist.
3. Launches one `x17_firstanalyzer` job for each input file.
4. Executes multiple jobs simultaneously (up to a configurable maximum number of parallel processes).
5. Waits for all jobs to finish.
6. Merges the individual output ROOT files into a single ROOT file using `hadd`.
7. Deletes the intermediate per-file ROOT files after a successful merge.
8. Prints a summary of successful and failed jobs.

## Input

The script expects reconstructed ROOT files following the naming convention:

```text
prad_<RUN>_recon_<FILE>.root
```

For example:

```text
prad_024014_recon_015.root
```

where:

* `RUN` is a six-digit run number.
* `FILE` is a three-digit file number.

The default input directory is:

```text
/volatile/hallb/prad/x17_replay/3cl_full_data/
```

## Output

Individual analysis outputs are written to:

```text
/work/hallb/prad/wrightso/x17_firstanalyzer/outfiles/<RUN>/
```

Each processed file produces:

```text
x17_Ana_<RUN>_<FILE>.root
```

After all jobs complete successfully, the script creates the merged output:

```text
x17_Ana_<RUN>_All.root
```

The individual ROOT files are then removed to conserve disk space.

## Requirements

The script requires:

* Python 3
* ROOT
* The ROOT `hadd` utility available in the system path
* The `x17_firstanalyzer` executable

By default, the executable is expected to be located at:

```text
/work/hallb/prad/wrightso/x17_firstanalyzer/x17_firstanalyzer
```

## Usage

Run the script by specifying the run number:

```bash
python3 Run_AnalysisBatch.py <RUN>
```

For example:

```bash
python3 Run_AnalysisBatch.py 24014
```

The maximum number of concurrent jobs can be changed using:

```bash
python3 Run_AnalysisBatch.py <RUN> --max-parallel <N>
```

For example:

```bash
python3 Run_AnalysisBatch.py 24014 --max-parallel 24
```

## Configuration

Several configuration parameters are defined near the top of the script:

| Parameter      | Description                                      |
| -------------- | ------------------------------------------------ |
| `MAX_PARALLEL` | Maximum number of concurrent analysis jobs.      |
| `EXECUTABLE`   | Path to the `x17_firstanalyzer` executable.      |
| `DATA_DIR`     | Directory containing reconstructed input files.  |
| `OUT_BASE`     | Base directory for output ROOT files.            |
| `RUN_WIDTH`    | Zero-padding width used for run numbers.         |
| `FILE_WIDTH`   | Zero-padding width used for output file numbers. |

These values can be modified to match the local analysis environment.

## Parallel Processing

The script uses Python's `ThreadPoolExecutor` to execute multiple instances of `x17_firstanalyzer` concurrently. Since each analysis job runs as an independent external process, multiple files can be analyzed simultaneously, substantially reducing the overall processing time for large data sets.

## Error Handling

Before processing begins, the script verifies that:

* The input directory exists.
* The `x17_firstanalyzer` executable is present and executable.
* Input ROOT files matching the requested run are found.

During execution, the script records the exit status and execution time for every analysis job. If any jobs fail, a summary table and the final lines of each job's error output are printed to assist with debugging.

If the `hadd` merge fails, the intermediate ROOT files are preserved to prevent data loss.

## Summary Output

Upon completion, the script reports:

* Run number
* Number of input files processed
* Number of successful and failed jobs
* Total wall-clock execution time
* Location of the merged ROOT file

This summary provides a quick overview of the batch processing results and highlights any files that require further investigation.

## Author Information

This script was originally written by **Rafayel Paremuzyan** ([rafopar@jlab.org](mailto:rafopar@jlab.org)) for automating batch execution of PRad-II Raw Sum Trigger Validation analysis from random trigger runs.

It was subsequently adapted and extended by **Erik Wrightson** ([wrightso@jlab.org](mailto:wrightso@jlab.org)) for use in this project.

---
# 2. Analysis PDF Plot Generator

## Overview

`PrintHistograms.py` is a Python utility that reads the output ROOT file produced by the `x17_firstanalyzer` and automatically generates comprehensive PDF reports summarizing the reconstructed event distributions. The script extracts all histograms from the input ROOT file, groups them by analysis type, and formats them into organized multi-page PDF documents.

The generated reports provide a convenient way to inspect the effects of each analysis cut, monitor detector performance, and evaluate candidate event distributions throughout the analysis chain.

## Workflow

Given an `x17_firstanalyzer` analysis ROOT file, the script performs the following steps:

1. Opens the ROOT file and recursively scans all directories for histograms.
2. Copies each histogram into memory to allow the input file to be safely closed.
3. Categorizes histograms according to their analysis type (e<sup>-</sup>X&rarr;e<sup>-</sup>(e<sup>+</sup>e<sup>-</sup>) final state candidates, γγe<sup>-</sup> final state candidates, and Møller candidates).
4. Optionally reads a JSON file containing the live charge (livetime*collected charge) for normalization.
5. Generates a detailed multi-page PDF for each analysis category.
6. Produces reduced PDF summaries containing only the most important diagnostic plots.
7. Calculates invariant-mass integrals around the 17 MeV region and annotates the corresponding plots.
8. Applies optional live-charge normalization to invariant-mass spectra when live-charge information is available.

## Input

The script expects a ROOT file containing the histograms produced by the X17 analysis framework.

For example:

```text
x17_Ana_024014_All.root
```

Optionally, a JSON file containing the measured live charge may also be supplied. The script expects the JSON file to contain the field

```json
{
    "value_nC": <live_charge_in_nanocoulombs>
}
```

This value is converted internally to mC and used to normalize invariant-mass distributions.

An optional fiducial geometry flag may also be specified to adjust the cut descriptions shown in the generated reports.

## Output

The script generates up to four PDF documents:

```text
<basename>_3Particle.pdf
```

A complete cut-by-cut report for the three-particle (e⁺e⁻e⁻) final state analysis.

```text
<basename>_3Particle_reduced.pdf
```

A condensed version of the three-particle report highlighting the primary diagnostic plots.

```text
<basename>_3Particle_Gamma_reduced.pdf
```

A reduced report for (γγe<sup>-</sup>) final state candidate reconstruction.

```text
<basename>_Moller.pdf
```

A complete report for the Møller scattering analysis.

Each report is organized into multiple pages documenting the progression of the analysis through successive event-selection cuts.

## Requirements

The script requires:

* Python 3
* CERN ROOT with Python bindings (PyROOT)
* A ROOT file generated by the X17 analysis framework

Optional:

* A JSON file containing the run live charge for normalization.

## Usage

Generate reports from an analysis ROOT file:

```bash
python3 PrintHistograms.py <ROOT_FILE>
```

For example:

```bash
python3 PrintHistograms.py x17_Ana_024014_All.root
```

To include live-charge normalization:

```bash
python3 PrintHistograms.py <ROOT_FILE> <LIVE_CHARGE_JSON>
```

For example:

```bash
python3 PrintHistograms.py x17_Ana_024014_All.root liveCharge_024014.json
```

To specify the fiducial geometry option:

```bash
python3 PrintHistograms.py <ROOT_FILE> <LIVE_CHARGE_JSON> <FIDUCIAL_OPTION>
```

For example:

```bash
python3 PrintHistograms.py x17_Ana_024014_All.root liveCharge_024014.json 1
```

## Generated Reports

### Three-Particle Candidate Report

The full three-particle report includes diagnostic plots for every analysis cut, including:

* HyCal cluster positions
* Energy v. polar angle (&theta;)
* Momentum distributions
* Total and individual cluster energies
* Timing distributions
* Coplanarity (Δ&phi;)
* Vertex position (Distance of closest appraoch [DCA] z)
* Invariant-mass spectra
* Live-charge-normalized invariant-mass spectra (when available)

Each cut stage is preceded by a description page explaining the applied event-selection criteria.

### Reduced Three-Particle Report

The reduced report presents a condensed version of the analysis, focusing on the primary monitoring histograms while minimizing redundant plots.

### γγe Candidate Report

The γγe report summarizes the analysis of candidates reconstructed as two photons and one electron. It follows the same overall structure as the three-particle report while documenting the γγ-specific event-selection requirements.

### Møller Candidate Report

The Møller report documents the event selection used for elastic Møller scattering, including:

* HyCal position distributions
* Momentum distributions
* Energy and timing spectra
* Coplanarity
* Vertex reconstruction
* Invariant-mass distributions

For later cut stages, the invariant-mass distributions are fit using a Crystal Ball function to characterize the reconstructed peak.

## Histogram Organization

The script automatically categorizes histograms based on their names within the ROOT file.

The primary categories are:

* `_X_` — Three-particle (e⁺e⁻e⁻) final state candidates
* `_Xgg_` — (γγe) final state candidate analysis
* `_M_` — Møller candidate analysis

All ROOT directories are scanned recursively, allowing the script to operate independently of the internal directory structure.

## Live-Charge Normalization

When a live-charge JSON file is supplied, invariant-mass histograms are scaled to units of counts per mC.

The script also computes and displays event yields integrated within several invariant-mass windows centered at 17 MeV:

* ±0.5 MeV
* ±1.0 MeV
* ±1.5 MeV
* ±2.5 MeV

These values are printed directly on the invariant-mass plots to facilitate comparisons between runs.

## Configuration

Several plotting options are defined within the script, including:

* Canvas layouts
* Plot styling
* Histogram grouping
* Cut descriptions
* Live-charge normalization
* Invariant-mass integration windows
* Crystal Ball fitting for Møller candidates

These parameters can be modified to customize the generated reports.

## Error Handling

Before processing begins, the script verifies that:

* The input ROOT file exists.
* The ROOT file can be opened successfully.
* Optional live-charge JSON files are valid before normalization is applied.

Histograms are detached from the input ROOT file before it is closed, ensuring that all plotting operations are performed safely in memory.

## Summary

This utility provides an automated method for producing publication-quality diagnostic reports from the X17 analysis output. By organizing hundreds of histograms into structured PDF documents, it greatly simplifies validation of event selection, detector performance, and reconstructed physics quantities while providing a consistent visualization format for every analyzed run.

## Author Information

This plotting utility was developed for the X17 analysis framework to automate the production of standardized PDF summaries from the analysis ROOT output.

The software recursively extracts histograms from the analysis file, organizes them by candidate type, computes live-charge-normalized invariant-mass distributions, and generates comprehensive cut-by-cut reports suitable for detector validation, analysis development, and physics studies.
