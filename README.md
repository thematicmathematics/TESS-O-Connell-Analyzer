
# TESS O'Connell Analyzer 

A fully automated Python pipeline for analyzing starspot evolution (the O'Connell effect) using TESS light curves.

## Scientific Background

**What is the O'Connell Effect?**
The O'Connell effect (O'Connell, 1951) refers to the unequal brightness of out-of-eclipse maxima (Max I and Max II) in the light curves of eclipsing binary star systems. This asymmetry is typically caused by magnetic spot activity (starspots) on the stellar surface or circumstellar mass transfer.

**Metrics Used:**
* **Global Diff:** The difference between Max II and Max I in the phase-folded light curve (incorporating all sector data).
* **Local Diff:** The average of individual measurements taken across each orbital cycle.
* **SNR > 5σ:** Confirms significant spot activity (Bevington & Robinson, 2003).
* **|Asym| < 0.5%:** Statistically considered as no significant spot activity (Ricker et al., 2015).

##  Table of Contents
* [Features](#-features)
* [Installation](#-installation)
* [Usage](#-usage)
* [Outputs](#-outputs)
* [Example Console Output](#-example-console-output)
* [References](#-references)

---

##  Features
* **Automated Data Acquisition:** Downloads data from the MAST archive via Lightkurve using the TIC ID.
* **Fully Automated Analysis:** Performs light curve cleaning, normalization, and phase folding.
* **Dual Metric Calculation:** Computes both global and local O'Connell effect metrics.
* **Publication-Ready Visuals:** Generates detailed analysis panels (.png and .eps) for each sector.
* **Evolution Tracking:** Creates a time-series graph displaying asymmetry trends.
* **Interactive Terminal:** User-friendly prompts for missing parameters.
* **Multi-Pipeline Support:** Automatically searches across SPOC, TESS-SPOC, and QLP pipelines.

---

##  Installation

1. Clone the repository to your local machine:
```bash
git clone [https://github.com/thematicmathematics/TESS-O-Connell-Analyzer.git](https://github.com/thematicmathematics/TESS-O-Connell-Analyzer.git)
cd TESS-O-Connell-Analyzer
```

2. Install the required packages:
```bash
pip install -r requirements.txt
```

---

##  Usage

You can run the program in two different ways:

**1. Command Line Arguments (Fast Execution)**
```bash
python main.py \
  --tic 12345678 \
  --period 3.1294323 \
  --epoch 2452964.6091 \
  --author SPOC \
  --outdir my_analysis
```

**2. Interactive Mode**
If you run the script without arguments, the program will guide you:
```bash
python main.py
```

*Example Terminal Output:*
```text
====================================
 TESS O'CONNELL ASYMMETRY ANALYSIS 
====================================
TIC Number: 12345678
Period (days): 3.16857
T0 Reference Time (BJD): 2452964.24896
Data Source / Author [AUTO]: spoc or tess-spoc or qlp 
```

---

##  Outputs

**1. Sector Analysis Panels (`S##_Evidence.png` & `.eps`)**
Each panel contains:
* **Top:** Time series data with Max I/II markers.
* **Bottom Left:** Phase-folded light curve and peak points.
* **Bottom Right:** Detailed metrics report (Global/Local O'Connell differences, SNR, minima shifts, and significance status).

**2. CSV Report (`TIC_XXXXXXXX_OConnell_Report.csv`)**
A table containing the numerical data for all sectors:
```csv
Sector,Global_Diff,Global_Err,Local_Diff,Error_Pct,Min1_Shift...
01,0.4523,0.0821,0.4201,0.1024,-0.0012...
02,0.5102,0.0756,0.4987,0.0892,-0.0015...
```

**3. Asymmetry Evolution Graph (`TIC_XXXXXXXX_Asymmetry_Evolution.png` & `.eps`)**
Displays the change in the O'Connell difference over time (or across sectors) for all observed sectors.

---

##  Example Console Output

```text
 Searching for data from TIC 393617775 on MAST servers...
  Trying author: TESS-SPOC...
 Data found in TESS-SPOC!
 2 sector(s) found. Downloading (This may take some time)...

 Data is being processed and panels are being created...
### Sector 8 completed: 
(Global: -0.308% | Shifts: +0.000, +0.000)
### Sector 35 completed: 
(Global: 0.124% | Shifts: -0.000, +0.000)

 The results were saved as a CSV file: TIC_9540226_Analysis\TIC_9540226_OConnell_Report.csv
 Operation completed successfully! All evidence panels have been saved to the 'TIC_9540226_Analysis' folder.
```

---

##  References
1. O'Connell, D. J. K. (1951). Publications of the Riverview College Observatory, 2(6), 85.
2. Ricker, G. R., et al. (2015). Journal of Astronomical Telescopes, Instruments, and Systems, 1(1), 014003.
3. Bevington, P. R., & Robinson, D. K. (2003). Data Reduction and Error Analysis for the Physical Sciences.
```

If you use this software in your research, please cite it as follows:

S. Ceren Çalışkan. (2024). TESS-O-Connell-Analyzer (v1.0.0). GitHub. https://github.com/thematicmathematics/TESS-O-Connell-Analyzer

