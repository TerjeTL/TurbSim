function [Init]  = initialize_algorithm(map,BaseFunc,k_q,Constraints,LSOptim)
%[Init]  = initialize_algorithm(map,BaseFunc,k_q,constraints,LSOptim)
%
% Creates all initial signals required for the parameterization process
% depending on the user inputs
%
%   map             - Struct containing the map signals
%   BaseFunc        - Struct containing the base function choice
%   k_q             - Heat transfer correction parameter
%   Constraints   	- Struct containing the parameter bound constraints
%   LSOptim       	- Boolean, true if LSOptim solver is selected
%
% Version 1.2: 2017-12-20
% Copyright (C) 2017, Xavier Llamas
%
% This file is part of LiU CPgui.
%
% LiU CPgui is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as
% published by the Free Software Foundation, version 3 of the License.
%
% This package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the74
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with LiU CPgui.  If not, see <http://www.gnu.org/licenses/>.
%
%% Create all necessary map signals
map             = get_signals(map,k_q);
%% Options
if isfield(map,'MODEL')
    map.MODEL(ismember(map.MODEL,' ,.:;!')) = '_'; % remove blank Spaces
    options.Name_save       = map.MODEL;
else
    options.Name_save       = 'Unknown_Compressor';
end
% Option to plot reverse mass flow model, requires a parameter: param.K0
options.SurgeReverse    = 0;    % Plot reverse flow area.
if options.SurgeReverse
    map.K0 = 0.3;
end
%% Choice of Base functions from user input
options.BaseFunc        = BaseFunc;  
%% Choice of Solver, LSOptim or lsqnonlin, Solver options are editable.
options.LSOptim         = LSOptim;        % 1 = LSOptim, 0 = lsqnonlin
if options.LSOptim
    % Solver Options
    options.n_eval_init     = 40;       % Maximum number of evaluations in LSOptim Initialization of base functions
    options.n_eval_2D       = 60;       % Maximum number of evaluations in LSOptim Massflow
    options.n_eval_3D       = 100;      % Maximum number of evaluations in LSOptim Complete model
else
    % Options for lsqnonlin solver
    options.BaseSolverOpt   = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','TolFun',1e-7, 'display', 'none','MaxIter',1000,'MaxFunEvals',50000); 
    options.SolverOptions   = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','TolFun',1e-7, 'display', 'iter','MaxIter',1000,'MaxFunEvals',50000); 
    options.SolverOptions3D = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','TolFun',1e-6, 'display', 'iter','MaxIter',1000,'MaxFunEvals',50000);

end
%% Gather all the remaining parameters into the Init Structure
options                     = get_options(options);
InitGuess                   = get_initialguess(map,options);
Constraints                 = get_constraints(map,options,Constraints);
Init.map                    = map;
Init.Options                = options;
Init.InitGuess              = InitGuess;
Init.Constraints            = Constraints;
end
%% Subfunctions

function [options] = get_options(options)   
% December 2016 Xavier Llamas
% Version 1.0
%  + 
%% Function handles
options.F_Wch           = str2func(options.BaseFunc.Wch);
options.F_PIch          = str2func(options.BaseFunc.PIch);
options.F_WZS           = str2func(options.BaseFunc.WZS);
options.F_PIZS          = str2func(options.BaseFunc.PIZS);
options.F_CUR           = str2func(options.BaseFunc.CUR);
options.F_PI0           = str2func(options.BaseFunc.PI0);
options.F_a             = str2func(options.BaseFunc.A);
options.F_b             = str2func(options.BaseFunc.B);
options.F_kloss         = str2func(options.BaseFunc.K_loss);
%% Size of the parameter vector for the Choke mass flow base function
if strcmp(options.BaseFunc.Wch,'F_WChL_atan')
    options.n_Wchparam 	= 4;
elseif strcmp(options.BaseFunc.Wch,'F_WChL_switch')
    options.n_Wchparam 	= 5;
else
    %Error
end

if strcmp(options.BaseFunc.A,'F_a_SAE')
    options.n_Aparam  	= 3;
elseif strcmp(options.BaseFunc.A,'F_a_ASME')
    options.n_Aparam   	= 2;
else
    %Error
end

if strcmp(options.BaseFunc.B,'F_b_SAE')
    options.n_Bparam  	= 2;
elseif strcmp(options.BaseFunc.B,'F_b_ASME')
    options.n_Bparam 	= 2;
else
    %Error
end

if strcmp(options.BaseFunc.K_loss,'F_kloss')
    options.n_Klossparam	= 1;
elseif strcmp(options.BaseFunc.K_loss,'F_kloss_none')
    options.n_Klossparam	= 0;
else
    %Error
end
options.n_PIchparam    	= 3;
options.n_WZSparam    	= 2;
options.n_PIZSparam    	= 2;
options.n_CURparam    	= 3;
options.n_PI0param  	= 1;
% Compute length of the parameter vector
options.n_param_ellipse = options.n_Wchparam+options.n_PIchparam+options.n_WZSparam+options.n_PIZSparam+options.n_CURparam+options.n_PI0param;
options.n_param_eff     = options.n_Bparam+options.n_Aparam+options.n_Klossparam;
options.n_param        	= options.n_param_ellipse+options.n_param_eff;

%% Declare Interval vector, used when parsing the parameter vector and the constraints.
options.ParvecI       	= zeros(9,1);
options.ParvecI(1)    	= options.n_Wchparam;
options.ParvecI(2)     	= options.ParvecI(1)+options.n_PIchparam;
options.ParvecI(3)    	= options.ParvecI(2)+options.n_WZSparam;
options.ParvecI(4)    	= options.ParvecI(3)+options.n_PIZSparam;
options.ParvecI(5)   	= options.ParvecI(4)+options.n_CURparam;
options.ParvecI(6)    	= options.ParvecI(5)+options.n_PI0param;
options.ParvecI(7)     	= options.ParvecI(6)+options.n_Bparam;
options.ParvecI(8)     	= options.ParvecI(7)+options.n_Aparam;
options.ParvecI(9)     	= options.ParvecI(8)+options.n_Klossparam;
%% Penalty factor for the verticality of the choking line during the TLS algorithm.
options.Pen            	= 0.01; % the smaller the number, the more vertical is the choke line during the parameterization
end

function [InitGuess] = get_initialguess(map,Options)
%% Initial Guess for the Ellipse Mass Flow model
if map.Wc_max_map<2
    % Vaneless Compressors Initial Guess
    if strcmp(Options.BaseFunc.Wch,'F_WChL_atan')
        InitGuess.param_Wch 	= [0.795,0.278,2.491,1.441]';
    elseif strcmp(Options.BaseFunc.Wch,'F_WChL_switch')
        InitGuess.param_Wch 	= [0.491,0.617,1.274,0.314,0.804]';
    else
        %Error
    end
    InitGuess.param_PIch   	= [0.109,0.387,3.493]';
    InitGuess.param_WZSL   	= [0.671,2.155]';
    InitGuess.param_PIZSL  	= [0.995,2.274]';
    InitGuess.param_CUR   	= [2.092,0.984,5.001]';
    InitGuess.C_s          	= 1;
    % Initial Guess and Constraints for the Efficiency Model
    if strcmp(Options.BaseFunc.A,'F_a_SAE')
        InitGuess.param_a  	= [0.403,0.0177,2.568]'; 
    elseif strcmp(Options.BaseFunc.A,'F_a_ASME')
        InitGuess.param_a  	= [0.757,0.216]'; 
    else
        %Error
    end

    if strcmp(Options.BaseFunc.B,'F_b_SAE')
        InitGuess.param_b  	= [1.022,0.0979]';
    elseif strcmp(Options.BaseFunc.B,'F_b_ASME')
        InitGuess.param_b	= [0.091,1.249]';
    else
        %Error
    end

    if strcmp(Options.BaseFunc.K_loss,'F_kloss')
        InitGuess.C_loss   	= 0.0114;
    elseif strcmp(Options.BaseFunc.K_loss,'F_kloss_none')
        InitGuess.C_loss    = [];
    else
        %Error
    end
else
    % Vaned Compressors Initial Guess (To be updated with more Vaned compressor data)
    if strcmp(Options.BaseFunc.Wch,'F_WChL_atan')
        InitGuess.param_Wch 	= [0.769,0.354,3.726,2.929]';
    elseif strcmp(Options.BaseFunc.Wch,'F_WChL_switch')
        InitGuess.param_Wch 	= [0.336,0.758,2.371,0.856,0.896]';
    else
        %Error
    end
    InitGuess.param_PIch 	= [0.066,0.641,2.931]';
    InitGuess.param_WZSL  	= [0.874,2.149]';
    InitGuess.param_PIZSL 	= [1.020,2.620]';
    InitGuess.param_CUR  	= [2.453,1.887,2.421]';
    InitGuess.C_s         	= 1;
    % Initial Guess and Constraints for the Efficiency Model ( From Vaneless Compressors, to be updated with more Vaned compressor data)
    if strcmp(Options.BaseFunc.A,'F_a_SAE')
        InitGuess.param_a  	= [0.311,0.071,5.209]'; 
    elseif strcmp(Options.BaseFunc.A,'F_a_ASME')
        InitGuess.param_a  	= [0.635,0.216]'; 
    else
        %Error
    end

    if strcmp(Options.BaseFunc.B,'F_b_SAE')
        InitGuess.param_b  	= [0.988,0.086]';
    elseif strcmp(Options.BaseFunc.B,'F_b_ASME')
        InitGuess.param_b	= [0,1.382]';
    else
        %Error
    end

    if strcmp(Options.BaseFunc.K_loss,'F_kloss')
        InitGuess.C_loss   	= 0.0161;
    elseif strcmp(Options.BaseFunc.K_loss,'F_kloss_none')
        InitGuess.C_loss   	= [];
    else
        %Error
    end    
end
%% Generate Choking mass flow and Choking pressure ratio for each SpL and zero
InitGuess.N_Ch                	= [0;map.Nc_vec];
InitGuess.N_Ch_n              	= InitGuess.N_Ch./map.Nc_max_map;
InitGuess.W_chk              	= Options.F_Wch(InitGuess.param_Wch,InitGuess.N_Ch_n,map);
InitGuess.PI_chk              	= Options.F_PIch(InitGuess.param_PIch,InitGuess.N_Ch_n,map);
%% Generate Zero Slope mass flow and Zero Slope pressure ratio for each SpL and zero
InitGuess.N_ZSL                	= [0;map.Nc_vec];
InitGuess.N_ZSL_n             	= InitGuess.N_ZSL./map.Nc_max_map;
InitGuess.W_ZSL                	= Options.F_WZS(InitGuess.param_WZSL,InitGuess.N_ZSL_n,map);
InitGuess.PI_ZSL               	= Options.F_PIZS(InitGuess.param_PIZSL,InitGuess.N_ZSL_n,map);  
%% Gererate Initial guess Curvature values for each SpL
InitGuess.CUR_i                 = Options.F_CUR(InitGuess.param_CUR,InitGuess.N_Ch_n);
end

function [Constraints] = get_constraints(map,Options,Constraints)
%% Ellipse Mass Flow Model Hard Bound Constraints
if map.Wc_max_map<2
    % Vaneless Compressors Hard constraints
    if strcmp(Options.BaseFunc.Wch,'F_WChL_atan')
        Constraints.C_Wch.lb  	= [0,0,0,0]';
        Constraints.C_Wch.ub  	= [5,5,5,5]';
    elseif strcmp(Options.BaseFunc.Wch,'F_WChL_switch')
        Constraints.C_Wch.lb  	= [0.35, 0, 1, 0, 0.5]';
        Constraints.C_Wch.ub  	= [0.7, inf, 2.5, inf, 1]';
    else
        %Error
    end
    Constraints.C_PIch.lb 	= [0.25/map.PI_max_map,0.7/map.PI_max_map,2]';
    Constraints.C_PIch.ub	= [0.8/map.PI_max_map,2/map.PI_max_map,8]';
    Constraints.C_Wzs.lb   	= [0.55 2]';
    Constraints.C_Wzs.ub  	= [1 3]';
    Constraints.C_PIzs.lb 	= [0 2]';
    Constraints.C_PIzs.ub 	= []';
    Constraints.C_cur.lb  	= [2,1e-4,1]';
    Constraints.C_cur.ub   	= [4,6,9]';
    Constraints.C_s.lb  	= 0.85; 
    Constraints.C_s.ub    	= 1.5;
else
    % Vaned Compressors Hard constraints
    if strcmp(Options.BaseFunc.Wch,'F_WChL_atan')
        Constraints.C_Wch.lb  	= [0,0,0,0]';
        Constraints.C_Wch.ub  	= [5,5,10,10]';
    elseif strcmp(Options.BaseFunc.Wch,'F_WChL_switch')
        Constraints.C_Wch.lb  	= [0.2, 0, 0, 0, 0.1]';
        Constraints.C_Wch.ub  	= [0.6,inf,4,inf,0.99]';
    else
        %Error
    end
    Constraints.C_PIch.lb 	= []';
    Constraints.C_PIch.ub 	= []';
    Constraints.C_Wzs.lb  	= [0.55 2]';
    Constraints.C_Wzs.ub  	= []';
    Constraints.C_PIzs.lb 	= [0 1]';
    Constraints.C_PIzs.ub 	= []';
    Constraints.C_cur.lb  	= [2,1e-4,1]';
    Constraints.C_cur.ub   	= [3.5,6,9]';
    Constraints.C_s.lb    	= 0.85; 
    Constraints.C_s.ub    	= 1.5;
end
%% Enthalpy Based Efficiency Hard Bound Constraints
if strcmp(Options.BaseFunc.A,'F_a_SAE')
    Constraints.C_a.lb 	= [0,0,1]';
    Constraints.C_a.ub 	= [inf,inf,25]';
elseif strcmp(Options.BaseFunc.A,'F_a_ASME')
    Constraints.C_a.lb 	= [-inf,-inf]';
    Constraints.C_a.ub	= []';
else
    %Error
end

if strcmp(Options.BaseFunc.B,'F_b_SAE')
    Constraints.C_b.lb  	= [0 -inf]';
    Constraints.C_b.ub  	= []';
elseif strcmp(Options.BaseFunc.B,'F_b_ASME')
 	Constraints.C_b.lb  	= [0,0]';
    Constraints.C_b.ub  	= []';
else
    %Error
end

if strcmp(Options.BaseFunc.K_loss,'F_kloss')
    Constraints.C_loss.lb  	= 0;
    Constraints.C_loss.ub 	= []';
elseif strcmp(Options.BaseFunc.K_loss,'F_kloss_none')
  	Constraints.C_loss.lb  	= 0;
    Constraints.C_loss.ub 	= []';
else
    %Error
end
end