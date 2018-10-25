-------------------------------------------------
          r-wave acoustic solver
-------------------------------------------------
This repository serves as supplementary material
for the article "Hamilton-Green solver for the 
forward and adjoint problems in Photoacoustic
Tomography".

In order to obtain the figures presented in the
paper please follow these steps:

1. Download the latest Matlab k-wave package from 
      
     http://www.k-wave.org/
 
2. Add k-wave and r-wave to your matlab path folder.

3. Change path to the folder Examples/

4. Run Example_kWave.m
This script generates the k-wave forward and 
adjoint data.

5. Run Example_RT_forward.m
This script generates the HG forward data for the 
smoothed phantom. It produces the figures
corresponding to the forward problem.

6. Run Example_RT_adjoint.m
This script generates the HG adjoint data and
the k-wave/HG adjoint data for mixed forward data.

7. Run Example_figures.
This script produces the figures for the PAT data 
and the adjoint problem.