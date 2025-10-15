[README.md](https://github.com/user-attachments/files/22923786/README.md)
# EdgeDAS: A Python-Based Edge Framework for Continuous Active Seismic DAS Monitoring

This repository contains the source code for **EdgeDAS**, an automated, onsite processing framework for Distributed Acoustic Sensing (DAS) data, as presented in the paper:

> *"EdgeDAS: A Python-Based Edge Framework for Continuous Active Seismic DAS Monitoring"*

EdgeDAS is designed to run on an **edge computing device** located onsite, enabling **near-real-time analysis** of continuous seismic data.

---

## 🧭 Overview

The framework provides a comprehensive solution for automated processing of continuous active-source DAS data.  
It leverages **parallel processing** to efficiently handle large data volumes and includes a **graphical user interface (GUI)** for easy configuration and operation.

The core workflow involves **cross-correlating continuous data** against a reference source signal to generate **virtual shot gathers**, which are then **stacked** to improve the signal-to-noise ratio (SNR).

---

## ⚙️ Key Features

- **GUI for Configuration** – Simple Tkinter-based interface for setting all processing parameters.  
- **Automated Workflow** – Automatically detects and processes new, unprocessed days of data.  
- **Parallel Processing** – Utilizes Python’s `multiprocessing` to process one hour of data per core.  
- **Ambient Noise Interferometry** – Cross-correlates continuous data against a known source signal.  
- **Advanced Correlation Options** – Supports both standard and one-bit (binary) correlation.  
- **Data Stacking** – Combines hourly shot gathers over user-defined periods (e.g., 24 hours).  
- **Automated Quality Control** – Generates and saves seismogram plots for visual inspection.  
- **Daily Scheduling** – Can be configured to run automatically at a specified daily time.  
- **Persistent Settings** – Remembers the last configuration by saving settings to a file.

---

## 📁 Repository Structure

```
EdgeDAS/
├── edgedas_framework.py        # Main Python application script
├── requirements.txt            # Required Python packages
├── .gitignore                  # Git ignore file
├── README.md                   # This documentation file
├── LICENSE                     # Project license
└── _log/
    └── source_pass.mseed       # Example/required reference source signal file
```

---

## 🧩 Installation

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/EdgeDAS.git
cd EdgeDAS
```

### 2. Create a Virtual Environment (Recommended)
```bash
python -m venv venv
source venv/bin/activate        # On Windows: venv\Scripts\activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

---

## 🚀 Usage

### 1. Prepare the Source Signal
Place your **reference source signal file** (`source_pass.mseed`) inside the `_log/` directory.  
This signal will be used for cross-correlation.

### 2. Run the Application
```bash
python edgedas_framework.py
```

### 3. Configure Parameters via the GUI
| Parameter | Description |
|------------|-------------|
| **Shift Seconds** | Time shift in seconds for aligning data windows. |
| **Shotting Interval (s)** | Time interval between consecutive virtual shots. |
| **Target Data Folder** | Directory containing raw DAS data organized by date (e.g., `YYYYMMDD`). |
| **Hours to Stack** | Number of hourly cross-correlation files to stack. |
| **Daily Run Time** | Time (`HH:MM`) for automatic daily runs. |
| **Cores to Use** | Number of CPU cores for parallel processing. |
| **Use Binary Cross-Correlation** | Enables one-bit cross-correlation. |

Click **“Run Processing & Schedule Daily”** to begin.  
The program will process all unprocessed data and schedule future runs automatically.

---

## 📡 Input Data Format

- Data should be stored in **daily subfolders** (e.g., `YYYYMMDD`), each containing **hourly subfolders**.  
- Each file is named with its precise start time (e.g., `YYYYMMDDHHMMSSffffff`).  
- The binary file format is as follows:
  - **8-byte header per channel block:** two 4-byte integers (`num_sections`, `num_samples`)
  - **Data block:** 4-byte floating-point seismic data  
- Default configuration:
  - 1200 channels  
  - Original sampling rate: 1000 Hz (downsampled to 250 Hz)

---

## 📊 Output

All processed data, plots, and logs are saved in the `_log/` directory.

| Path | Description |
|------|--------------|
| `_log/processing_log.txt` | Log of all processing activities (newest entries on top). |
| `_log/donedays.txt` | List of processed days. |
| `_log/processing_settings.txt` | Saved GUI configuration. |
| `_log/Saved_data/Cross_Correlation_folder/` | Hourly cross-correlation files (MSEED format). |
| `_log/Saved_data/Stacked_Cross_Correlation/` | Final stacked shot gathers (MSEED format). |
| `_log/Saved_data/Seismogram_Plots/` | Final stacked seismogram plots (.png). |

---

## 🧾 How to Cite

If you use this code in your research, please cite:

> *[paper's citation information here once published]*

---

## 🪪 License

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

---

## 👨‍💻 Author

**Ahmad Bahaa**  
Project Researcher, University of Tokyo  
📧 [ahmadbahaa@g.ecc.u-tokyo.ac.jp  ORCID:0000-0002-5374-9379]  
🌐 [https://github.com/AhmadBahaa1995]((https://github.com/AhmadBahaa1995))

---

*EdgeDAS – Bridging edge computing and continuous seismic monitoring for next-generation DAS analysis.*
