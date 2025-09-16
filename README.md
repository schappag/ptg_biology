# Parathyroid Gland Systems Biology (Python Implementation)

## Overview
This repository provides a Python implementation of the quantitative systems physiology model of parathyroid gland (PTG) biology published in its basic form in:

> Schappacher-Tilp G, Bieglmayer C, Binder C, Eberle C, Rudas M, Tilp M.  
> *A mathematical model of parathyroid gland biology.*  
> Physiological Reports. 2019; 7(11): e14089. https://doi.org/10.14814/phy2.14089
> 
The model describes the regulation of parathyroid hormone (PTH) secretion, intracellular degradation, synthesis, clearance, and proliferation incorporating effects of extracellular ionozed calcium, phosphate, and calcitriol.  
 It is particularly focused on the altered PTG biology found in patients with chronic kidney disease (CKD) on hemodialysis and is designed to serve as a tool for studying secondary hyperparathyroidism and exploring treatment strategies

## Key Features
• Comprehensive Physiology: Captures the complex network regulating PTH, including the effects of ionized calcium (Ca2+), phosphate, and 1,25-dihydroxyvitamin D3 (1,25D).
• CaSR-Centric Model: The model's core is the Calcium-Sensing Receptor (CaSR), the crucial regulator of PTG function and biology.
• Multi-Timescale Dynamics: Accurately simulates the different timescales of PTG adaptation, from the release of stored PTH within seconds to cellular proliferation over days and weeks.
• Adaptive Mechanisms: Explicitly models key adaptations such as changes in intracellular PTH degradation, PTH production rate, and PTG cellular proliferation.
• Validated Predictions: The model's predictions have been validated against published experimental data for various scenarios, including acute hypocalcemia, hysteresis effects, and the development of secondary hyperparathyroidism [1,2].
• Extensible Framework: Designed to be combined with models of medications (e.g., calcimimetics, [3]) or other physiological systems, such as bone remodelling, to study chronic kidney disease-mineral and bone disorder (CKD-MBD) in greater detail

## Community Activity
We welcome feedback and collaboration:  
- Use **Issues** to suggest new applications (e.g., drug effects).  
- Report bugs or model inconsistencies.  
- Propose extensions such as parameter estimation routines, PK/PD model integration, or new regulatory pathways.  

## Implementation
- Implemented in **Python**
- Includes parameter sets and initial conditions matching published simulations.  
- Example notebooks (`/notebooks`) show:  
  - Baseline homeostasis  
  - Acute calcium clamp experiments  
  - Chronic kidney disease simulation  
  - Pharmacological interventions (e.g., calcimimetics, vitamin D analogues)  

## Installation
```bash
git clone https://github.com/schappag/ptg-biology.git
cd ptg-biology
pip install -r requirements.txt



