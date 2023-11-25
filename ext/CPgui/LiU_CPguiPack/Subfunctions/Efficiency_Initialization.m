function [paramT_i] = Efficiency_Initialization(map, Param_mf, InitialGuess, Constraints, Options)
% [paramT_i] = Efficiency_Initialization(Data, Param_mf, InitialGuess, Constraints, Options)
%
% Provides an initial guess of the efficiency parameters by fitting each SpL 
% independently and then parameterizing the base functions.
%
%   map             - Struct containing the map signals
%   Param_mf    	- Vector containing the Ellipse model parameters
%   InitialGuess   	- Struct containing the initial guess of the base
%                     functions
%   Constraints   	- Struct containing the parameter bound constraints
%   Options       	- Struct containing the solver options
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with LiU CPgui.  If not, see <http://www.gnu.org/licenses/>.
%
%% Overall Optimization
%Gather Initial parameters
Param_init      =[InitialGuess.param_b;InitialGuess.param_a;InitialGuess.C_loss];
% Create upper and lower bounds
param_lb        = -inf.*ones(Options.n_Aparam+Options.n_Bparam+Options.n_Klossparam,1);
param_ub        = inf.*ones(Options.n_Aparam+Options.n_Bparam+Options.n_Klossparam,1);
% B bound constaints
if isempty(Constraints.C_b.lb)   
else
    param_lb(1:1+Options.n_Bparam-1) = Constraints.C_b.lb;
end
if isempty(Constraints.C_b.ub)   
else
    param_ub(1:1+Options.n_Bparam-1) = Constraints.C_b.ub;
end
%% A bound constaints
if isempty(Constraints.C_a.lb)   
else
    param_lb(1+Options.n_Bparam:Options.n_Bparam+Options.n_Aparam) = Constraints.C_a.lb;
end
if isempty(Constraints.C_a.ub)   
else
    param_ub(1+Options.n_Bparam:Options.n_Bparam+Options.n_Aparam) = Constraints.C_a.ub;
end

%% K_loss bound constaints
if isempty(Constraints.C_loss.lb)   
else
    param_lb(Options.n_Bparam+Options.n_Aparam+1:Options.n_Bparam+Options.n_Aparam+Options.n_Klossparam) = Constraints.C_loss.lb;
end
if isempty(Constraints.C_loss.ub)   
else
    param_ub(Options.n_Bparam+Options.n_Aparam+1:Options.n_Bparam+Options.n_Aparam+Options.n_Kfricparam) = Constraints.k_fric.ub;
end
%%
Residuals   = @(param)(residuals_efficiency_TLS(param,map,Options));
% Solve nonlinear Least squares
% Choose Solver
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(param_lb,param_ub,0);
    Param_eff   = lsoptim(Residuals,Param_init,Options.n_eval_init,A,b);
else
    Param_eff  	= lsqnonlin(Residuals,Param_init,param_lb,param_ub,Options.SolverOptions);
end
% Return Parameter Struct
Param_eff       = [Param_mf.vector;Param_eff];
In              = Options.ParvecI;
if Options.n_Klossparam == 0
    paramT_i.vector   	= Param_eff;
    paramT_i.delta_Wc 	= Param_mf.delta_Wc;
    paramT_i.C_Wch    	= Param_eff(1:In(1));
    paramT_i.C_PIch    	= Param_eff(In(1)+1:In(2));
    paramT_i.C_Wzs     	= Param_eff(In(2)+1:In(3));
    paramT_i.C_PIzs   	= Param_eff(In(3)+1:In(4));
    paramT_i.C_cur    	= Param_eff(In(4)+1:In(5));
    paramT_i.Gamma_PIcs	= 0.5;
    paramT_i.C_s        = Param_eff(In(5)+1:In(6));
    paramT_i.C_b      	= Param_eff(In(6)+1:In(7));
    paramT_i.C_a      	= Param_eff(In(7)+1:In(8));
    paramT_i.C_loss     = 0;
else
    paramT_i.vector   	= Param_eff;
    paramT_i.delta_Wc 	= Param_mf.delta_Wc;
    paramT_i.C_Wch   	= Param_eff(1:In(1));
    paramT_i.C_PIch    	= Param_eff(In(1)+1:In(2));
    paramT_i.C_Wzs     	= Param_eff(In(2)+1:In(3));
    paramT_i.C_PIzs   	= Param_eff(In(3)+1:In(4));
    paramT_i.C_cur   	= Param_eff(In(4)+1:In(5));
    paramT_i.Gamma_PIcs	= 0.5;
    paramT_i.C_s        = Param_eff(In(5)+1:In(6));
    paramT_i.C_b      	= Param_eff(In(6)+1:In(7));
    paramT_i.C_a      	= Param_eff(In(7)+1:In(8));
    paramT_i.C_loss     = Param_eff(In(8)+1:In(9));
end
paramT_i.F_Wch          = Options.F_Wch;
paramT_i.F_PIch         = Options.F_PIch;
paramT_i.F_WZS          = Options.F_WZS;
paramT_i.F_PIZS         = Options.F_PIZS;
paramT_i.F_CUR          = Options.F_CUR;
paramT_i.F_PI0          = Options.F_PI0;
paramT_i.F_a            = Options.F_a;
paramT_i.F_b            = Options.F_b;
paramT_i.F_kloss      	= Options.F_kloss;
paramT_i.T_ref          = map.T_ref;
paramT_i.p_ref          = map.p_ref;
paramT_i.Nc_max_map     = map.Nc_max_map;
paramT_i.Wc_max_map     = map.Wc_max_map;
paramT_i.PI_max_map     = map.PI_max_map;
paramT_i.Delta_h_max    = map.Delta_h_max;
paramT_i.gamma_air      = map.gamma_air;
paramT_i.rho1           = map.rho1;
paramT_i.D_2            = map.D_2;
paramT_i.Cp_air         = map.Cp_air;
paramT_i.alpha_kPiW0    = 0.15*(map.PI_max_map/map.Wc_max_map);
end

function [r] = residuals_efficiency_TLS(param,map,Options)
%[r] = residuals_efficiency_TLS(param,map,Options)
% March 2016 Xavier Llamas
% Version 1.0
%  + 
%% Compute modeled efficiency
eta_c_model         = F_eta_c_enthalpy_TLS(param,map,Options);
eta_c_vec               = [];
eta_c_model_vec         = [];
    for i = 1:length(map.Nc_vec)
            eta_c_curr          = map.etaC_M(:,i);
            eta_c_curr          = eta_c_curr(isnan(eta_c_curr)==false);
            eta_c_model_curr    = eta_c_model(:,i);
            eta_c_model_curr    = eta_c_model_curr(isnan(eta_c_model_curr)==false);
            eta_c_vec           = [eta_c_vec; eta_c_curr];
            eta_c_model_vec     = [eta_c_model_vec; eta_c_model_curr];
    end
Nc_vec      = map.NcCorr_M(isnan(map.NcCorr_M(:))<1);
weight      = ones(length(Nc_vec),1);
Nc_vec_n    = Nc_vec./map.Nc_max_map;
%weight(Nc_vec_n>map.weight_switch) = map.weight;
r = (eta_c_vec-eta_c_model_vec).*weight;
end

function [eta_c_model] = F_eta_c_enthalpy_TLS(param,map,Options)
%[eta_c_model] = F_eta_c_enthalpy_TLS(param,map,Options)
% July 2016 Xavier Llamas
% Version 1.0
%  + 
%% Compute Variables
B_model         = Options.F_b(param(1:Options.n_Bparam),map.NcCorr_M./map.Nc_max_map,map);
A_model         = -Options.F_a(param(Options.n_Bparam+1:Options.n_Bparam+Options.n_Aparam),map.NcCorr_M./map.Nc_max_map,map);
K_fric_model    = Options.F_kloss(param(end),map.NcCorr_M, map.WcCorr_M,map);
D_H_model       = K_fric_model.*(B_model+A_model.* map.WcCorr_M);
Delta_h_is      = map.Cp_air.*map.T01_M.*(map.PiC_M.^((map.gamma_air - 1)./map.gamma_air) - 1);
eta_c_model     = Delta_h_is./D_H_model;
end

function [A,b]= CreateLinearConstraints(param_lb,param_ub,m)
% [A,b]= CreateLinearConstraints(param_lb,param_ub,m)
%
% Creates linear system of constrains given the upper and lower constraint
% vectors and the size of the Deltas in the TLS problem. Infs in the bounds
% are removed.
%
%   param_lb                    - Lower bound on the parameters, length n
%   param_ub                    - Upper bound on the parameters, length n
%   m                           - Length of the Deltas in the TLS problem
%   A                           - Matrix of the Linear system
%   b                           - B vector of the linear system
%
% March 2016 Xavier Llamas
% Version 1.0
%  + 
%% Get the output
Lower_bound         = isempty(param_lb);
Upper_bound         = isempty(param_ub);
Deltas              = m>1;
Both_Constraints    = Lower_bound*Upper_bound;
if isempty(param_lb)
    A_lb    = [];
    b_lb    = [];
    nred_lb = 0;
else
    n           = length(param_lb);
    pos         = find(not(isinf(param_lb)));  
    param_lb    = param_lb(isinf(param_lb)<1);
    nred_lb     = length(param_lb);
    A_lb        = zeros(nred_lb,n);
    for i = 1:nred_lb
        A_lb(i,pos(i))        = 1;
    end
    b_lb        = param_lb;
end
if isempty(param_ub)
    A_ub    = [];
    b_ub    = [];
    nred_ub = 0;
else
    n           = length(param_ub);
    pos         = find(not(isinf(param_ub)));  
    param_ub    = param_ub(isinf(param_ub)<1);
    nred_ub     = length(param_ub);
    A_ub        = zeros(nred_ub,n);
    for i = 1:nred_ub
        A_ub(i,pos(i))        = -1;
    end
    b_ub        = -param_ub;
end

if Both_Constraints<1
    A           = [A_lb;A_ub];
    b           = [b_lb;b_ub];
elseif (Both_Constraints>0)&&(Lower_bound<1)
    A           = A_lb;
    b           = b_lb;
elseif (Both_Constraints>0)&&(Upper_bound<1)
    A           = A_ub;
    b           = b_ub;
else
    A           = [];
    b           = [];
end

if (Deltas>0)&&(Both_Constraints<1)
    A_deltas    = zeros(nred_lb+nred_ub,m); 
    A           = [A,A_deltas];
elseif (Deltas>0)&&(Both_Constraints>0)
    A_deltas    = zeros(nred_lb+nred_ub,m); 
    A           = [A,A_deltas];
else
end
end