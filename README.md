# Parathyroid Gland Systems Biology (Python Implementation)

## Aim
This repository provides a Python implementation of the quantitative systems physiology model of parathyroid gland (PTG) biology published in its basic form in:

> Schappacher-Tilp G, Bieglmayer C, Binder C, Eberle C, Rudas M, Tilp M.  
> *A mathematical model of parathyroid gland biology.*  
> Physiological Reports. 2019; 7(11): e14089. https://doi.org/10.14814/phy2.14089
> 
The model describes the regulation of parathyroid hormone (PTH) secretion, intracellular degradation, synthesis, clearance, and proliferation incorporating effects of extracellular ionozed calcium, phosphate, and calcitriol.  
 It is particularly focused on the altered PTG biology found in patients with chronic kidney disease (CKD) on hemodialysis and is designed to serve as a tool for studying secondary hyperparathyroidism and exploring treatment strategies
 
## Community Activity
We welcome feedback and collaboration:  
- Use **Issues** to suggest new applications (e.g., drug effects, CKD progression, gland hypertrophy).  
- Report bugs or model inconsistencies.  
- Propose extensions such as parameter estimation routines, PK/PD model integration, or new regulatory pathways.  

## Background and Motivation
The parathyroid gland plays a central role in calcium–phosphate homeostasis.  
The model captures key mechanisms:  
- Nonlinear Ca–PTH secretion dynamics (sigmoidal set-point response).
- Direct effect of phosphate on the Calcium sensing receptor (CaSR)
- PTH synthesis and intracellular pool regulation.  
- PTH clearance.  
- Long-term gland growth and turnover.  

Originally designed to reproduce short- and long-term PTH responses, the framework enables virtual clinical trial simulations and translational studies of parathyroid-related disorders.

## Implementation
- Implemented in **Python**, using `scipy.integrate` for ODE solving.  
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



