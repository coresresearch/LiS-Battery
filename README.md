# LiS-Battery Model using Cantera

The model contained in this repository allows the user to run simulations of a lithium-sulfur cell. The following list shows the main sections of the README. You can click the hyperlink to skip directly to a section of interest.

1. [Repository contents](#repository-contents)
2. [Getting started](#getting-started)
3. [Using the model](#using-the-model)
4. [Further work on model](#further-work-on-model)


# Repository contents

The Li-S model is composed of five primary python files and at least one yaml file in order to be used. The yaml file is the Cantera input file that provides all of the relevant species, phase, thermodynamic, and kinetic parameters. 

The five python files are:

- `li_s_battery_inputs.py`: this file is where parameters are provided to set geometry, initial conditions, Cantera input file name, phase names (used to import information from the Cantera input), etc. It is also used to set simulation conditions such as C-rate, sulfur loading, and the nucleation density of the active cathode phases.
- `li_s_battery_init.py`: this file takes the inputs from `li_s_battery_inputs` and organizes them into class structures for convenient carrying of parameters in the model. It also completes any pre-processing that is needed, such as calculating relevant values such as total initial sulfur mass, external current based on C-rate input, and calculates any additional values needed based on input parameters. Some of these can also be overridden (such as initial volume fractions).
- `li_s_battery_model.py`: this is the primary model file that contains the numerical solver initialization, the residual function passed to the solver (this contains the governing equations of the model), and can write the data produced to a `.csv` file if desired.
- `li_s_battery_functions.py`: this file contains any helper functions used repeatedly in the `li_s_battery_model.py` file, such as transport calculations.
- `li_s_battery_post.py`: this file performs post-processing to produce summary plots after each run of `li_s_battery_model.py` and to process the data in Pandas dataframes for easy inspection and further processing.

As mentioned previously, at least one `.yml` file will need to be provided, as this is the input file used by Cantera to create relevant phases, species, and reactions. The repository contains multiple examples, the primary one being `Bessler_Dennis_lithiated.yml`. For more information on Cantera, please visit the Cantera website [here](https://cantera.org/tutorials/index.html)


# Getting started

In order to use the model, several python packages are necessary. The key packages with version numbers used are listed below.

- assimulo (version 2.9)
- numpy (version 1.15.2)
- pandas (version 0.23.4)
- python (version 3.5.6)
- cantera (version 2.5.0a4)

Information on installing Assimulo, the DAE solver used here, can be found [here](https://jmodelica.org/assimulo/download.html#)


# Using the model

Once a working environment has been created with the necessary packages, the next step is to provide a `.yml` input file that Cantera can use to import the required phases, species, interfaces, and reactions. Several input files are included in this repository and can be used as templates. The structure and phase names provided in the input file must match the strings assigned to phase names in the `li_s_battery_inputs.py` file. 

General steps to run the model:

1. Open `li_s_battery_inputs.py` and set desired input parameters. Primary parameters that can be altered are:
  - `C_rate`: the C-rate is used to set the external current in conjunction with the sulfur loading
  - `H_cat`: the cathode thickness
  - `A_C_0`: the pristine interface area between the carbon and electrolyte phases
  - `A_cat`: planar area of the cell
  - `m_S_0`: sulfur loading (mass per area) or bulk (mass) depending on `sulfur_method` string value
  - Other geometric and material parameters can be set for the cell components in this file as well
2. The model can be run from `li_s_battery_inputs.py` directly, or by opening and running `li_s_battery_model.py` where solver parameters can be modified if desired.
3. After the model is finished running, it will print summary plots of voltage, solid phase volume fractions, and species concentrations. Additionally, the data will be saved in the current folder as a `.csv` file.

# Further work on model

Work is being done to further generalize the code and move all parameters users may wish to adjust directly into the `li_s_battery_inputs.py` file, including solver settings, etc. Additionally, further updates to the model will likely alter the numerical solver package used to allow for more flexibility.
