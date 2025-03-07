# assr-ord
Auditory Steady State Response simulation and objective response detection in time and frequency domains.

Current status (Feb 12th, 2025): 
I'm recycling old code into interoperable models.

## SinalEEG_Analysis 
This unified MATLAB application: 
- Loads or simulates EEG data 
- Preprocesses the signal (DC removal, bandpass filtering) 
- Computes FFT and magnitude-squared coherence (MSC) 
- Computes beta-distribution based sequential test thresholds 
- Runs a multi-stage sequential test with early stopping 
- Evaluates performance metrics and visualizes the results 

> Developed based on modern software engineering practices. 

> Author: [Alexandre Gomes Caldeira] 

### Update log
- Version: [Date]
    > Commit message or Comment 

- v0.0: [14 02 2025 10:40] 
    > Project created, exposing incompatibilities...

- v0.1: [07 03 2025 10:45] 
    > DataLoader with FFT functional for exp and sim! :)

- v0.2: [07 03 2025 16:42]
    > Preprocessing functional for exp and sim, signal recompute added

- v0.3: [] 
    > MSC calculator now accepts variable parameters

- v0.4: [] 
    > ORD calculator functional for exp and sim

- v0.5: [] 
    > SHT (single hyp. test [Decisions, Time]) functional for exp and sim 

- v0.6: [] 
    > add Metrics: [Pareto, ConfusionMatrix, Specificity, Repeatability]

- v0.7: [] 
    > GST (group sequential test) functional for exp and sim

- v0.8: [] 
    > GST metrics verified 

- v0.9: []
    >
- v1: [] 
    >