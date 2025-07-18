# Stochastic Modeling for COVID-19: Self-Reporting and Contact Tracing Effectiveness in Heterogeneous Risk Groups

A stochastic modeling study of COVID-19 outbreak in Daegu, Korea using Modified Gillespie Algorithm. This model analyzes the effectiveness of self-reporting and contact tracing strategies across heterogeneous risk groups.

---

## Model Overview

### Population Stratification
- High-risk group: Religious community members (Shincheonji-related cases)  
  ◦ Characteristics: High transmissibility, low self-reporting rate  
- Low-risk group: General community population  
  ◦ Characteristics: Lower transmissibility, higher reporting compliance

### Compartmental Structure
- S: Susceptible population  
- C: Contacts under surveillance  
- E1/E2: Exposed population (distinguished by self-reporting behavior)  
- I1/I2: Infectious population (distinguished by self-reporting behavior)  
- Q: Quarantined population  
- R: Recovered population



## File Structure

### Main Execution Files
- Main_ModiGillespie.m: Main simulation execution script  
- Parameters.m: Model parameter settings  
- ModiGillespie_algorithm.m: Core Gillespie algorithm implementation  
- ContactTracing_Dice.m: Contact tracing module  

### Data Files
- Daegu_COVID19_data.xlsx: Real COVID-19 data from Daegu region



## Key Parameters

### Disease Parameters
- Latent period: Gamma distribution (mean: 5.48 days, SD: 2.72 days)  
- Infectious period: Uniform distribution (7.39-17.39 days)  
- Self-reporting delay: Log-normal distribution  

### Contact and Transmission Parameters
- Contact rates: Differentiated before/after social distancing  
- Transmission probability: High-risk 0.721, Low-risk 0.118  

### Contact Tracing Parameters
- Self-reporting rate: High-risk 0.4, Low-risk 0.8  
- Contact tracing delay: Uniform distribution (4-7 days)  
- Isolation period: 14 days  

---

## Quick Start

### Prerequisites
- MATLAB R2018b or later  
- Parallel Computing Toolbox (optional, for parallel processing)  

### Basic Execution
```matlab
% Run main simulation
run('Main_ModiGillespie.m')

% Modify simulation settings in Main_ModiGillespie.m
num_of_run = 1000;  % Change to desired number of simulations
```

### Parameter Adjustment
Key variables in Parameters.m:
```matlab
input.rhoH = 0.4;           % High-risk self-reporting rate
input.rhoL = 0.8;           % Low-risk self-reporting rate
input.t1 = 4; input.t2 = 7; % Contact tracing delay range
input.SDtime = 11;          % Social distancing start time (day 11)
```



## Results and Analysis

### Output Variables
- Hcase: Cumulative confirmed cases in high-risk group  
- Lcase: Cumulative confirmed cases in low-risk group  
- 95% Credible intervals for both groups  
- Distribution histograms of final case counts  

### Key Findings
- 22% increase in total infections when high-risk self-reporting drops to 0.1  
- 164% increase in unreported cases under low self-reporting scenarios  
- 85% increase in unreported cases with 4-7 day contact tracing delays  
- Contact tracing effectiveness diminishes significantly beyond 4-day delays  



## Scenario Analysis

### Scenario 1: Varying High-risk Group Self-reporting
- Fix low-risk group at 0.8 self-reporting rate  
- Vary high-risk group from 0.1 to 1.0  
- Assess impact on both groups  

### Scenario 2: Varying Low-risk Group Self-reporting
- Fix high-risk group at 0.4 self-reporting rate  
- Vary low-risk group from 0.1 to 1.0  
- Evaluate cross-group transmission effects  

### Contact Tracing Delay Analysis
- Five delay intervals: 1-4, 4-7, 7-10, 10-13, 13-16 days  
- Uniform distribution within each interval  
- 10,000 simulations per scenario  



## Key Results Summary

```
Scenario              High-risk Cases      Low-risk Cases       Unreported Cases
-------------------------------------------------------------------------------------
Baseline              6,191 (5,774-6,580)  2,733 (2,305-3,188)  2,168 (628-3,693)
Low reporting (0.1)   +22%                 +46%                 +164%
High reporting (0.8)  -21%                 -36%                 -86%
Delayed tracing (4-7d)+15%                 +21%                 +85%
```

---

## Publication

This code accompanies the research paper:  
**"Quantitative analysis of self-reporting and contact tracing in heterogeneous risk groups: a stochastic modeling study of the COVID-19 outbreak in Daegu, Korea"**  

Authors: Jiwon Han, Eunok Jung  
Department of Mathematics, Konkuk University, Seoul, Republic of Korea



## Data Source

Epidemiological data sourced from Korea Disease Control and Prevention Agency (KDCA) public press releases covering the period February 18 - March 23, 2020.  
Data availability: KDCA Press Release Archive



## Citation

If you use this code in your research, please cite:  
Han, J., Jung, E. (2024). *Quantitative analysis of self-reporting and contact tracing in heterogeneous risk groups: a stochastic modeling study of the COVID-19 outbreak in Daegu, Korea.* [Journal TBD]



## Funding

This research is supported by:
- Korea National Research Foundation (NRF) grant (NRF-2021M3E5E308120711)  
- Government-wide R&D to Advance Infectious Disease Prevention and Control (HG23C1629)



## Contact

For questions about model implementation or usage, please open an issue in this repository.

---

**Keywords**: COVID-19, stochastic modeling, heterogeneous risk groups, self-reporting, contact-tracing, non-pharmaceutical interventions
