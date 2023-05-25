# Inference of Process Variations in Silicon Photonics from Characterization Measurements
This is the code and data repo for the paper of same name in Optics Express. _(TODO: add reference and link after published)_

## Data folder
`data/` folder contains original measurement data, simulation results, and some intermediate results from the extraction.
### Measurement data
`data/ring_data.mat` contains the original measurement (TE) data of all 48 RRs, including:
- `data`: matrix with each column corresponding to IL measurement (dB) of one device, and each row corresponding to one wavelength point
- `lambda`: vector specifying the wavelength points (nm) in the rows of `data`
- `location`: matrix with 2 columns (x, y), with each row specifying the location of the device in the column of the same index in `data`
- `radius`: vector specifying the ring radius (um) at location x
- `gap`: vector specifying the gap width (um) at location y

`data/TM_data.mat` is the measurement in TM, with the same format as `data/ring_data.mat`.

### Simulation results
`data/bend.mat` contains raw simulation data from the bent waveguise simulation from MODE.
- `neff`, `ng`: effective index and group index. first dimension corresponds to wavelength, second dimension corresponds to radius, thir dimension corresponds to variations
- `lambda`: list of wavelengths (m) corresponds to the first dimension of data
- `radius`: list of radius (um) corresponds to thr second dimension of data
- variation as the third dimension, in the order: nominal, dw=65nm, dw=-65nm, dT=10nm, dT=-10nm, dh=3.5nm, dh=-3.5nm

`data/coupling.mat` contains raw simulation data for the coupling coefficients
- `T`: (complex) coupling coefficient. first dimension corresponds to wavelength, second dimension corresponds to radius and gap width, third dimension corresponds to variations
- `lambda`: same as in `data/bend.mat`
- `DoE`: list of radius and gap width (um) corresponds to the second dimension of data. radius as the first column, gap wisth as the second.
- `ratio`: list of variations (as scaling ratio to the nominal value) corresponds to the third dimension of data. columns in order: w/w0, T/T0, h/h0

`data/coupling_base.mat` contains raw simulation data of a reference waveguide for coupling coefficient calculation.
- `Tb`: (complex) reference through port coefficient. first dimension corresponds to wavelength, second dimension corresponds to variations.
- `lambda`: same as in `data/bend.mat`
- ratio`: same as in `data/coupling.mat`

### Intermediate results
- `data/intGaussGrid.mat`: from `intGauss.m`
- `data/extract_raw.mat`: from `feature_extract.m`
- `data/extract_TM.mat`: from `TMinfo.m`
- `data/extract_coupling_loss.mat`: from `coupling_loss_assign.m`
- `data/extract_final.mat`: from `coef_save.m`
- `data/inference_input.mat`: from `sim_analysis.m`
- `data/baseline.mat`: this is the baseline response from smoothing the response from RRs with no observed resonance.

## Simulation scripts
- bend waveguide simulation in MODE: `sim/wg_draw.lsf`, `sim/wg_2D.lsf`, `sim/wg_sweep_bend.lsf`
- coupling simulation in FDTD: `sim/rr_start.lsf`, `sim/dc_draw.lsf`, `sim/dc_setup.lsf`, `sim/dc_run.lsf`, `sim/rr_sweep_ridge.lsf`, `sim/wg_ref.lsf`
## Code files
The running order of the code files and their functions:

0. `intGauss.m`: _(optional)_ produces the pre-computed grid data for fast discrete Gaussian distribution calculation
1. `feature_extract.m`: peak parameters fitting and error estimation, also include some initial parameter conversion
2. `TMinfo.m`: extracts the peak information from the TM data. Even though the whole system is set in TE mode, we observed some conversion to TM mode at certain resonant peaks, possibly from some variations. These TM mode leakage can affect background noise level of TE mode and cause inaccurate extinction ratio and thus coupling & loss. Therefore, we extract the peak information from the TM data, and exclude the coupling & loss extraction data where the TM peak is significant (but effective & group index extraction are still included).
3. `coupling_loss_assign.m`: an automatic procedure (similar to cluster problem) to assign correct values for coupling and loss based on the wavelength response and similarity across devices. For detail, please refer to my PhD thesis. _(TODO: add link when available)_
4. `coef_save.m`: calculates effective and group index from peak information.
5. `sim_analysis.m`: generates sensitivity matrices and nominal fitting design matrices from simulation data. Also performs adjustment scaling of error estimation.
6. `EM_analysis.m`: EM algorithm of the Bayesian inference.

### Variable names in `EM_analysis.m`
For numerical computation purpose, we are using a different system for the variables used in the EM algorithm in `EM_analysis.m`:
