# BedrockCoreModel
This code accompanies the manuscript Balter-Kennedy et al., (submitted), 
"Orbital- and centennial-scale erosion rates beneath the Greenland Ice Sheet", submitted
to Journal of Geophysical Research: Earth Surface.

A general workflow for using these scripts is as follows:

- If you would like to initialize the model with non-zero Be-10 concentrations, you can 
use the script make_steady_state_profiles.m to calculate Be-10 concentrations with depth 
under user-prescribed steady erosion rates. Pre-calculated profiles for subaerial erosion
rates of 5, 10, and 50 m/Myr are provided here as .mat files in the SteadyStateProfiles 
folder.

- Set up and run the model using the wrapper script JAK_CR1_CoreModelRun.m. It should be 
relatively straightforward to adapt this for a different field site / dataset. 

- CoreModel.m is the actual model. All inputs are defined in JAK_CR1_CoreModelRun.m. Outputs
are the final Be-10 concentrations and depths. Output for relevant samples depths will be saved 
in .mat files.

- Figures can be generated using JAK_CR1_CoreModelFigureGenerator.m. 
