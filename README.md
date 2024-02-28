# Wind Turbine Operational Optimization
Evaluation/optimization framework for fatigue damage and revenue accumulation using adaptive control based on electricity prices and wind conditions. It includes the work done in the PhD thesis "Wind turbine operational optimization considering revenue and fatigue damage objectives" at the University of Stuttgart and related publications by Vasilis Pettas.

The main script is **_Evaluation_Optimization_Framework.m_**. Here the user can define the input data and choose the different options allowing for simulation with different operational approaches. The various operational modes and optimization approaches can be selected and adjusted in the script's **INPUTS** and **Control mode options** sections. The out-of-the-box setup allows users to test all modes with the predefined configurations and familiarize themselves with the available settings. 

Additional user-defined cases can be implemented by adding a new case in the list and creating a new optimization function that should be included in the folder **Objective_functions**. The settings for the different optimizers considered can be defined within each objective function. 

The folder **Wind_Price_data** includes the datasets used in the PhD thesis. New data to be evaluated can be added here following the file format of the existing files. Scripts are provided to create the statistics files by post-processing result files from the aeroelastic software FAST v8 (also compatible with openFAST).  

The folder **Functions** includes the functions required by the main script to run.

The folder **Standalone_tools** includes scripts and functions relevant for post-processing and visualising result files generated with the framework. Additionally, the tools to create the simulation database by post-processing the FAST simulation outputs are provided. 

The framework is developed and tested with Matlab 2021b on a Windows machine. 

Details about the different options and definitions can be found in the comments included in the source code.


# Related repositories

The surrogate models including the Gaussian process regression model (not included here) and relevant code for creating and generating the prediction of the surrogates can be found at https://doi.org/10.5281/zenodo.10092271

Result files from the application of the framework presented in the PhD thesis and other publications using the datasets provided here can be found at https://doi.org/10.5281/zenodo.10580236


# Related publications

Pettas V., Wind turbine operational optimization considering revenue and fatigue damage objectives, Ph.D. dissertation Wind Energy, University of Stuttgart, 2024 http://dx.doi.org/10.18419/opus-13959

Pettas V., Cheng P. W.: Surrogate modeling and aeroelastic response of a wind turbine with down-regulation, power boosting and individual blade pitch control capabilities, submitted, Energies 2024

Pettas V., Cheng P. W.: Operational optimization of wind turbines for revenue and fatigue objectives considering wind conditions and electricity prices, in preparation, Wind Energy 2024

Kölle K., Göçmen T., T., Eguinoa, I., Alcayaga Roman, L. A., Aparicio-Sanchez, M., Feng, J., Meyers, J., Pettas, V., and Sood, I.: FarmConners market showcase results: wind farm flow control considering electricity prices, Wind Energ. Sci., 7, 2181–2200, https://doi.org/10.5194/wes-7-2181-2022,2022.
