# Radial-basis-TMPC
Learning-based robust tube-based MPC of nonlinear systems via difference of convex radial basis functions approximations. 


<br />
<p align="center">
   <img src="https://github.com/martindoff/Radial-basis-TMPC/blob/main/RBF.png" alt="Logo" width="400" height="300">
  <p align="center">
   Difference-of-convex-functions (DC) decomposition of system dynamics via radial basis functions (RBF) approximations. 
    <br />  
  </p>
</p>

<!-- ABOUT THE PROJECT -->
## About The Project

Learning-based robust tube-based MPC of dynamic systems approximated by difference-of-convex (DC) radial basis functions (RBF) models. Successive linearisations of the learned dynamics in DC form are performed to express the MPC scheme as a sequence of convex programs.  Convexity in the learned dynamics is exploited to bound successive linearisation errors tightly and treat them as bounded disturbances in a robust MPC scheme. Application to the coupled tank problem. This novel computationally tractabe tube-based MPC algorithm is presented in the paper "Safe Learning in Nonlinear Model Predictive Control" by Johannes Buerger, Mark Cannon and Martin Doff-Sotta. It is an extension of our previous work [here] (https://github.com/martindoff/DC-TMPC) and [here](https://ora.ox.ac.uk/objects/uuid:a3a0130b-5387-44b3-97ae-1c9795b91a42/download_file?safe_filename=Doff-Sotta_and_Cannon_2022_Difference_of_convex.pdf&file_format=application%2Fpdf&type_of_work=Conference+item)

### Built With

* MATLAB
* CVX
* MOSEK

## Running the code

1. Clone the repository
   ```sh
   git clone https://github.com/martindoff/Radial-basis-TMPC.git
   ```
2. In the MATLAB command window, run

   ```matlab
   convex_anmpc_main
   ```
