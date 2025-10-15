# EdgeDAS: A Python-Based Edge Framework for Continuous Active Seismic DAS Monitoring

This repository contains the source code for **EdgeDAS**, an automated, onsite processing framework for Distributed Acoustic Sensing (DAS) data.

EdgeDAS is designed to run on an edge computing device located onsite, enabling near-real-time analysis of seismic data. It provides a robust, configurable, and automated workflow for processing continuous active-source DAS recordings.

---

## Overview

The framework provides a comprehensive solution for the automated processing of continuous active-source DAS data. It leverages parallel processing to handle large data volumes efficiently and includes a graphical user interface (GUI) for easy configuration and operation. The primary workflow involves cross-correlating continuous data against a reference source signal to generate virtual shot gathers, which are then stacked to improve the signal-to-noise ratio.

---

## Key Features

- **GUI for Configuration:** A simple interface built with Tkinter for setting processing parameters.  
- **Flexible Source Signal Handling:**
  - Use a pre-prepared source signal file (`source_pass.mseed`).
  - Automatically extract the source signal from a specified DAS channel during processing.
- **Automated Workflow:** Automatically discovers and processes new, unprocessed days of data by comparing available data against a log of completed days.
- **Parallel Processing:** Uses Python’s multiprocessing library to distribute the cross-correlation workload across multiple CPU cores, with each core processing one hour of data.
- **Advanced Correlation Options:** Supports both standard and binary (one-bit) cross-correlation methods.
- **Data Stacking:** Stacks hourly shot gathers over user-defined periods (e.g., 24 hours) to enhance signal quality.
- **Automated Quality Control:** Generates and saves seismogram plots of the stacked data for visual inspection.
- **Process-Safe Logging:** Uses Python’s standard logging module with a queue-based listener to safely handle log messages from multiple concurrent processes.
- **Performance Monitoring:** Logs detailed performance metrics (duration, CPU usage, memory consumption) for each processing stage into a structured `performance_log.csv` file.
- **Daily Scheduling:** Can be scheduled to run automatically at a specific time each day.
- **Persistent Settings:** Remembers the last used configuration by saving and loading settings from a file.
- **utilizes.py (optional):** A supplementary module containing utility functions that can be integrated into the main workflow if needed. These utilities are not required for the basic operation of the EdgeDAS framework but can assist with additional data handling, diagnostics, or visualization tasks.
---

## Repository Structure

```
EdgeDAS/
├── Edgedas_framework.py        # Main Python application script
├── requirements.txt            # Required Python packages
├── README.md                   # This documentation file
├── utilizes.py                 # Optional utility functions (for custom extensions)
├── LICENSE                     # Project license
├── _log/
│   ├── source_pass.mseed       # (Optional) Pre-prepared source signal
│   ├── donedays.txt            # Log of processed days
│   ├── processing_log.txt      # Main application log file
│   └── processing_settings.txt # Saved GUI settings
└── _Processed_Data/
    ├── performance_log.csv         # Performance metrics log
    ├── Cross_Correlation_folder/   # Hourly cross-correlation files
    ├── Stacked_Cross_Correlation/  # Final stacked data files
    ├── Seismogram_Plots/           # Final seismogram plots
    └── Source_Signals/             # Extracted source signals (if using extract mode)
```

---

## Installation

### Clone the repository:
```bash
git clone https://github.com/your-username/EdgeDAS.git
cd EdgeDAS
```

### Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
```

### Install the required dependencies:
```bash
pip install -r requirements.txt
```

---

## Usage

### Prepare the Source Signal (Choose one option):

**Option A:** Pre-prepared File — Place your reference source signal file named `source_pass.mseed` inside the `_log` directory.  
**Option B:** Extract from Data — Specify which DAS channel contains the source signal.

### Run the Application:
```bash
python Edgedas_framework.py
```

### Configure Parameters via the GUI:
- **Target Data Folder:** Root directory containing raw DAS data (e.g., YYYYMMDD folders).  
- **Shift Seconds:** Time shift in seconds for aligning data windows.  
- **Shooting Interval (s):** Time interval in seconds between consecutive virtual shots.  
- **Hours to Stack:** Number of hourly cross-correlation files to stack.  
- **Source Function:**
  - Select “Use Pre-prepared File” to use `_log/source_pass.mseed`.
  - Select “Extract from DAS Data” and enter a channel number.  
- **Advanced Settings:**
  - **Daily Run Time:** The time (HH:MM) for automatic daily processing.
  - **Cores to Use:** Number of CPU cores for multiprocessing.
  - **Bandpass Filter:** Minimum and maximum frequencies (Hz).
  - **Use Binary Cross-Correlation:** Optional one-bit correlation mode.

### Start Processing:
Click **“Run Processing & Schedule Daily”** to process all unprocessed data and schedule automatic runs.

---

## Input Data Format

The script reads custom binary-format DAS data organized in daily/hourly subfolders.  
Each channel block includes:
- 8-byte header (two 4-byte integers): `num_sections`, `num_samples`
- Raw seismic data as 4-byte floats

Configured for 1200 channels, 1000 Hz sampling (downsampled to 500 Hz).

---

## Output

All processed results and logs are stored in `_Processed_Data` and `_log`:

- `_Processed_Data/performance_log.csv` – Performance metrics  
- `_Processed_Data/Cross_Correlation_folder/` – Hourly MSEED cross-correlations  
- `_Processed_Data/Stacked_Cross_Correlation/` – Final stacked gathers  
- `_Processed_Data/Seismogram_Plots/` – Seismogram PNG plots  
- `_Processed_Data/Source_Signals/` – Extracted source signals (if applicable)  
- `_log/processing_log.txt` – Full process log  
- `_log/donedays.txt` – List of processed days  
- `_log/processing_settings.txt` – Saved GUI configuration

---

## How to Cite
please cite:

> *[paper's citation information here once published]*

---

## License

© 2025 EdgeDAS Project. Licensed under the MIT License.
See the [LICENSE](LICENSE) file for details.

---

## Author

**Ahmad Bahaa**  
Project Researcher, University of Tokyo  
📧 [ahmadbahaa@g.ecc.u-tokyo.ac.jp  ORCID:0000-0002-5374-9379]  
🌐 [https://github.com/AhmadBahaa1995]((https://github.com/AhmadBahaa1995))

---

*EdgeDAS – Bridging edge computing and continuous seismic monitoring for next-generation DAS analysis.*
