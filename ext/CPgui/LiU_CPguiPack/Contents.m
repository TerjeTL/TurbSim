% LiU Compressor Parameterization Graphical User Interface package 
% ---- LiU CPgui -------------------------------------------
% ==========================================================
%       Version 1.2: 2017-12-20
%       Copyright (C) 2017, Xavier Llamas
%
%       LiU CPgui is free software: you can redistribute it and/or modify
%       it under the terms of the GNU Lesser General Public License as
%       published by the Free Software Foundation, version 3 of the License.
%
%       This package is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%       GNU Lesser General Public License for more details.
%
%       You should have received a copy of the GNU Lesser General Public License
%       along with LiU CPgui.  If not, see <http://www.gnu.org/licenses/>.
% ==========================================================
% Package Contents
% ==========================================================
% init_LiU_CPgui            	- Initializes the compressor parameterization package and sets the required path.
% example_LiU_CPgui            	- Example call of the CP GUI with a example compressor map, created to resemble a measured one, the user is supposed to load the compressor map in the provided fields.
% ExampleMap                    - mat file which contains the example compressor map (not measured from any real compressor).                 
% LiU_CPgui                   	- Compressor parameterization Graphical User Interface.

% Model Implementations
% ==========================================================
% CompressorModel              	- m-function implementation of the compressor model in the forward implementation.
% CompModel                    	- Simulink implementation of the compressor model in usual forward mode.
% CompModel_Surge               - Simulink implementation of the compressor model for surge simulation in the More-Greitzer framework.
% BaseFunc_Blocks             	- Simulink library of the different model base functions, implemented to be used in case the user selects different base functions in the LiU CPgui tool.

% Subfunctions used inside the GUI
% ==========================================================
% lsoptim                    	- Solves an uncostrained non-linear least squares problem using a Levenberg-Marquardt like method.
% initialize_algorithm          - Creates all initial signals required for the parameterization process depending on the user inputs
% get_signals                   - Creates all initial signals required for the given compressor map, in thegiven map struct
% get_vector                    - Returns a vector of the mean of each column in the Matrix M
% get_mapStruct                 - Creates the map structure required by the parameterization GUI as input
% Ellipse_Initialization        - Returns the Ellipse parameters given the Initial guess, the Constraints and solver options
% Ellipse_Parameterization      - Solves a Total Least Squares problem to optimize the Ellipse parameters given the Initial guess, the Constraintsand solver options
% Efficiency_Initialization    	-  Provides an initial guess of the efficiency parameters by fitting each SpL independently and then parameterizing the base functions.
% Complete_Parameterization   	- Solves a Total Least Squares problem to optimize the complete compressor model parameters given the Initial guess, the Constraintsand solver options

% Plotting Subfunctions
% ==========================================================
% plot_modelMassFlow            - Plots the model pressure ratio vs mass flow in the given axes handle
% plot_modelEfficiency          - Plots the model efficiency vs mass flow in the given axes handle
% plot_modelContour             - Plots the model pressure ratio vs mass flow in the given axes handle with the efficiency extrapolation as contour levels.
% plot_model3D                  - Plots a number of speed lines in the 3D space formed by mass flow, pressure ratio and efficiency
% plot_lambda                   - Plots the work coefficient vs flow coefficient or the enthalpy vs mass flow for the given map and parameter sets
% plot_baseFunctions            - Plots each model base function vs corrected speed
% compute_grid              	- Creates the compressor mass flow grid area given the model parameters and the options
% compute_PiC                   - Computes the pressure ratio given the gridded compressor mass flow area and the model parameters
% compute_etaC                  - Computes the efficiency given the gridded compressor mass flow and pressure ratio areas and the model parameters
% compute_basefunctions      	- Computes the base functions in a dense grid given the model parameters and the gridding options

% Base functions
% ==========================================================
% F_WZSL                        - Computes the Zero Slope Line Mass Flow given the parameters and normalized corrected speed vector
% F_WChL_switch                 - Computes the choke line mass flow given the parameters and normalized corrected speed vector
% F_WChL_atan                   - Computes the choke line mass flow given the parameters and normalized corrected speed vector
% F_PIZSL                       - Computes the Zero Slope Line pressure ratio given the parameters and normalized corrected speed vector
% F_PIChL                       - Computes the choking pressure ratio given the parameters and normalized corrected speed vector
% F_PI0                         - Computes the pressure ratio at zero flow given the parameters and normalized corrected speed vector
% F_kloss_none                  - Computes k_loss values given the parameters and the corrected speed and mass flow values
% F_kloss                       - Computes k_loss values given the parameters and the corrected speed and mass flow values
% F_CUR                         - Computes CUR values given the parameters and the normalized corrected speed vector
% F_b_SAE                       - Computes B value given the parameters and the normalized corrected speed vector
% F_b_ASME                      - Computes B value given the parameters and the normalized corrected speed vector
% F_a_SAE                       - Computes A value given the parameters and the normalized corrected speed vector
% F_a_ASME                      - Computes A value given the parameters and the normalized corrected speed vector


% If you found this package useful for your research, please cite our work using one or more of the references below.
%
% More information:
% ==========================================================
% Llamas, X. and Eriksson, L., Control-Oriented Compressor Model with Adiabatic Efficiency Extrapolation, 
% SAE Int. J. Engines 10(4):2017, doi:10.4271/2017-01-1032.
% 
% Llamas X, Eriksson L. Parameterizing Compact and Extensible Compressor Models Using Orthogonal Distance Minimization. 
% ASME. J. Eng. Gas Turbines Power. 2016;139(1):012601-012601-10. doi:10.1115/1.4034152. 
%
% Eriksson, L. , and Nielsen, L. , 2014, Modeling and Control of Engines and Drivelines, Wiley, Hoboken, NJ.

% To do in future releases:
% ==========================================================
% - Show model absolute relative errors for the selected parameters in the GUI.
%
% - Plot the current feasible space for the model parameterization based on the current parameter bound constraints
%
% - More solver options like maxiter, tolerances, etc.. in the GUI panel.
