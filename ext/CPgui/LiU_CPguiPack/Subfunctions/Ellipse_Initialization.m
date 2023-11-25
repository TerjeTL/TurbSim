function [paramEll_i] = Ellipse_Initialization(map, InitialGuess, Constraints, Options)
% [param_I] = Ellipse_Initialization(map, InitialGuess, Constraints, Options)
%
% Returns the Ellipse parameters given the Initial guess, the Constraints
% and solver options
%
%   map             - Struct containing the map signals
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
%% Fit Base functions
if map.n_SpL>2
    %% Fit Base functions to the Initial guess vectors if there are more than 2 SpL available in the map
    [param_Wch,	 map.Wch_model]    	= Fit_F_WChL(map, InitialGuess, Constraints, Options);
    [param_Pich, map.PIch_model]   	= Fit_F_PIChL(map, InitialGuess, Constraints, Options);
    [param_Wzsl, map.WZSL_model]   	= Fit_F_WZSL(map, InitialGuess, Constraints, Options);
    [param_PIzsl,map.PIZSL_model]   = Fit_F_PIZSL(map, InitialGuess, Constraints, Options);
  
else
    param_Wch       = InitialGuess.param_Wch;
    param_Pich      = InitialGuess.param_PIch;
    param_Wzsl      = InitialGuess.param_WZSL;
    param_PIzsl     = InitialGuess.param_PIZSL;
    map.Wch_model   = Options.F_Wch(param_Wch,map.Nc_vec./map.Nc_max_map,map);
    map.PIch_model  = Options.F_PIch(param_Pich,map.Nc_vec./map.Nc_max_map,map);
    map.WZSL_model  = Options.F_WZS(param_Wzsl,map.Nc_vec./map.Nc_max_map,map);
    map.PIZSL_model = Options.F_PIZS(param_PIzsl,map.Nc_vec./map.Nc_max_map,map);
end
%% Get Pi0 with current parameters and the initial guess
map.PI0_model   = Options.F_PI0(param_PIzsl,map.Nc_vec./map.Nc_max_map,map);
%% Nonlinear least squares to find each speed line curvature  and surge parameter // CUR C_s
[C_s_vec,map.CUR_vec,Delta_Wc_init] = Fit_SpL_CUR(map, InitialGuess, Constraints, Options);
%% F_CUR Least squares fitting
param_CUR   = Fit_F_CUR(map, InitialGuess, Constraints, Options);
%% Output 
C_s0    = max(Constraints.C_s.lb ,mean(C_s_vec(end-1:end)));
C_s0    = min(C_s0,Constraints.C_s.ub );
Param_init  = [param_Wch; param_Pich; param_Wzsl; param_PIzsl;param_CUR;C_s0];
In       	= Options.ParvecI;
paramEll_i.vector   	= Param_init;
paramEll_i.delta_Wc 	= Delta_Wc_init;
paramEll_i.C_Wch      	= Param_init(1:In(1));
paramEll_i.C_PIch   	= Param_init(In(1)+1:In(2));
paramEll_i.C_Wzs     	= Param_init(In(2)+1:In(3));
paramEll_i.C_PIzs     	= Param_init(In(3)+1:In(4));
paramEll_i.C_cur      	= Param_init(In(4)+1:In(5));
paramEll_i.Gamma_PIcs  	= 0.5;
paramEll_i.C_s          = Param_init(In(5)+1:In(6));
paramEll_i.F_Wch        = Options.F_Wch;
paramEll_i.F_PIch      	= Options.F_PIch;
paramEll_i.F_WZS        = Options.F_WZS;
paramEll_i.F_PIZS     	= Options.F_PIZS;
paramEll_i.F_CUR        = Options.F_CUR;
paramEll_i.F_PI0        = Options.F_PI0;
paramEll_i.C_cur_vec  	= map.CUR_vec;
paramEll_i.T_ref        = map.T_ref;
paramEll_i.p_ref        = map.p_ref;
paramEll_i.Nc_max_map   = map.Nc_max_map;
paramEll_i.Wc_max_map  	= map.Wc_max_map;
paramEll_i.PI_max_map  	= map.PI_max_map;
paramEll_i.Delta_h_max 	= map.Delta_h_max;
paramEll_i.gamma_air  	= map.gamma_air;
paramEll_i.rho1         = map.rho1;
paramEll_i.D_2          = map.D_2;
paramEll_i.Cp_air       = map.Cp_air;
paramEll_i.alpha_kPiW0  = 0.15*(map.PI_max_map/map.Wc_max_map);
end
%% Subfunctions
function [param_Wch,Wch_model] = Fit_F_WChL(map, InitialGuess, Constraints, Options)
res_W_chk   = @(param)(InitialGuess.W_chk-Options.F_Wch(param,InitialGuess.N_Ch_n,map));
% Choose Solver
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(Constraints.C_Wch.lb ,Constraints.C_Wch.ub,0);
    param_Wch   = lsoptim(res_W_chk,InitialGuess.param_Wch,Options.n_eval_init,A,b);
else
    param_Wch   = lsqnonlin(res_W_chk,InitialGuess.param_Wch,Constraints.C_Wch.lb,Constraints.C_Wch.ub,Options.BaseSolverOpt);
end
Wch_model   = Options.F_Wch(param_Wch,map.Nc_vec./map.Nc_max_map,map);
end

function [param_Pich,PIch_model] = Fit_F_PIChL(map, InitialGuess, Constraints, Options)
res_PI_chk  = @(param)(InitialGuess.PI_chk-Options.F_PIch(param,InitialGuess.N_Ch_n,map));
% Choose Solver
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(Constraints.C_PIch.lb,Constraints.C_PIch.ub,0);
    param_Pich  = lsoptim(res_PI_chk,InitialGuess.param_PIch,Options.n_eval_init,A,b);
else
    param_Pich  = lsqnonlin(res_PI_chk,InitialGuess.param_PIch,Constraints.C_PIch.lb,Constraints.C_PIch.ub,Options.BaseSolverOpt);
end
PIch_model  = Options.F_PIch(param_Pich,map.Nc_vec./map.Nc_max_map,map);
end

function [param_Wzsl,WZSL_model] = Fit_F_WZSL(map, InitialGuess, Constraints, Options)
res_W_ZSL   = @(param)(InitialGuess.W_ZSL-Options.F_WZS(param,InitialGuess.N_ZSL./map.Nc_max_map,map));
% Choose Solver
if Options.LSOptim
    [A,b]           = CreateLinearConstraints(Constraints.C_Wzs.lb ,Constraints.C_Wzs.ub,0);
    param_Wzsl	= lsoptim(res_W_ZSL,InitialGuess.param_WZSL,Options.n_eval_init,A,b);
else
    param_Wzsl	= lsqnonlin(res_W_ZSL,InitialGuess.param_WZSL,Constraints.C_Wzs.lb,Constraints.C_Wzs.ub,Options.BaseSolverOpt);
end
WZSL_model  = Options.F_WZS(param_Wzsl,map.Nc_vec./map.Nc_max_map,map);
end

function [param_PIzsl,PIZSL_model] = Fit_F_PIZSL(map, InitialGuess, Constraints, Options)
res_PI_ZSL  = @(param)(InitialGuess.PI_ZSL-Options.F_PIZS(param,InitialGuess.N_ZSL./map.Nc_max_map,map));
% Choose Solver
if Options.LSOptim
    [A,b]               = CreateLinearConstraints(Constraints.C_PIzs.lb,Constraints.C_PIzs.ub,0);
    param_PIzsl    = lsoptim(res_PI_ZSL,InitialGuess.param_PIZSL,Options.n_eval_init,A,b);
else
    param_PIzsl    = lsqnonlin(res_PI_ZSL,InitialGuess.param_PIZSL,Constraints.C_PIzs.lb,Constraints.C_PIzs.ub,Options.BaseSolverOpt);
end
PIZSL_model = Options.F_PIZS(param_PIzsl,map.Nc_vec./map.Nc_max_map,map);
end

function [C_s_vec,CUR_vec,Delta_Wc_init] = Fit_SpL_CUR(map, InitialGuess, Constraints, Options)

%% Separated Total Least Squares to find each speed line curvature  and surge parameter // CUR C_s
CUR_vec           	= zeros(length(map.Nc_vec),1);
C_s_vec         = zeros(length(map.Nc_vec),1);
Delta_Wc_init       = zeros(map.n_measured,1);
index_Delta_curr    = 1;
for i = 1:length(map.Nc_vec)
    % Filter NaNs from the map data into the current SpL vector
    Wc_corr_curr    = map.WcCorr_M(:,i);
    Wc_corr_curr    = Wc_corr_curr(isnan(Wc_corr_curr)==false);
    PiC_curr        = map.PiC_M(:,i);
    PiC_curr        = PiC_curr(isnan(PiC_curr)==false);  
    n_SpL_curr      = length(Wc_corr_curr);
    % Current Choke and ZSL values
    Wch_curr        = map.Wch_model(i);
    PIch_curr       = map.PIch_model(i);
    PIZSL_curr      = map.PIZSL_model(i);
    WZSL_curr       = map.WZSL_model(i);
    PI0_curr        = map.PI0_model(i);
    % Initial guess and bounds from the guessed data
    deltaWc_0       = -1e-3.*ones(length(Wc_corr_curr),1);
    deltaWc_0(Wc_corr_curr>Wch_curr) = (Wch_curr-Wc_corr_curr(Wc_corr_curr>Wch_curr))*1.1;
    deltaWc_0(PiC_curr<PIch_curr) = (Wch_curr-Wc_corr_curr(PiC_curr<PIch_curr))*0.8;
    param0_C_s  = 1;
    param0_CUR      = [InitialGuess.CUR_i(i);param0_C_s;deltaWc_0];  
    Residuals_CUR   = @(param)(Residuals_CUR_TLS(param,Wc_corr_curr,PiC_curr,WZSL_curr,Wch_curr,PIch_curr,PIZSL_curr,PI0_curr,Options));
    if Options.LSOptim
        param_CUR_curr = lsoptim(Residuals_CUR,param0_CUR,Options.n_eval_init,[],[]);
    else
        param_CUR_curr = lsqnonlin(Residuals_CUR,param0_CUR,[],[],Options.BaseSolverOpt);
    end
    % Store the results
    CUR_vec(i)      = real(param_CUR_curr(1)); 
    C_s_vec(i)  = real(param_CUR_curr(2)); 
    C_s_vec(i)  = max(0.7,C_s_vec(i));
    Delta_Wc_init(index_Delta_curr:index_Delta_curr+n_SpL_curr-1)   = param_CUR_curr(3:end);
    index_Delta_curr = index_Delta_curr+n_SpL_curr;
end
end

function [r] = Residuals_CUR_TLS(param,Wc_corr_curr,PiC_curr,WZSL_curr,Wch_curr,PIch_curr,PIZSL_curr,Pi0_curr,Options)
% [r] = Residuals_CUR_TLS(param,Wc_corr_curr,PiC_curr,WZSL_curr,Wch_curr,PIch_curr,PIZSL_curr)
%
% July 2016 Xavier Llamas
% Version 1.0
%  + 
%% Compute intersection points.
Pic = f_TLS_CUR(param,Wc_corr_curr,WZSL_curr,Wch_curr,PIch_curr,PIZSL_curr,Pi0_curr,Options);
%% Compute squared distance  normalized with the maximun in the current SpL
% compute the maximums
n_points    = length(PiC_curr);
r           = (Pic-PiC_curr);
r           = [r;param(3:end)];
W_norm      = max(Wc_corr_curr).*ones(n_points,1);
Pi_norm     = max(PiC_curr).*ones(n_points,1);
Norm        = [Pi_norm;W_norm];
r           = r./Norm;
end

function [Pic] = f_TLS_CUR(param,Wc_corr_curr,WZSL_curr,Wch_curr,PIch_curr,PIZSL_curr,PI0_curr,Options)
%
% July 2016 Xavier Llamas
% Version 1.0
%  + 
%% Current parameters for the given SpL
CUR_SpL       	= param(1);
C_s       	= param(2); 
%% Initialize errors in variables as parameters
Wc_inner       	= Wc_corr_curr + param(3:end);
%% Check wether a map point is outside the range and compute the intersection in another manner
in_limit_Wch  	= Wc_inner<=Wch_curr;
in_limit_Wzs   	= Wc_inner>WZSL_curr;
in_limit_surge 	= Wc_inner<=WZSL_curr;
in_zone        	= in_limit_Wch.*in_limit_Wzs;
%% Return ellipse evaluation
PiC_i_curr_e  	= PIch_curr+(PIZSL_curr-PIch_curr).*(1-((Wc_inner(in_zone>0)-WZSL_curr)./(Wch_curr-WZSL_curr)).^CUR_SpL).^(1./CUR_SpL);
%% Return surge polynomial evaluation
W_tilde        	= WZSL_curr.*(1-(1-(Wc_inner(in_limit_surge>0)./WZSL_curr)).^C_s).^(1./C_s);
PiC_i_curr_p  	= PI0_curr +3.*((PIZSL_curr-PI0_curr)./WZSL_curr.^2).*W_tilde.^2-2.*((PIZSL_curr-PI0_curr)./WZSL_curr.^3).*W_tilde.^3;
%% Concatenate vectors
PiC_i_curr     	= zeros(length(Wc_corr_curr),1); 
PiC_i_curr(in_zone>0)        = PiC_i_curr_e;
PiC_i_curr(in_limit_surge>0) = PiC_i_curr_p;
%% Linear model when Wc_inner > Wch_curr, as mush steep as possible to resemble a vertical choking line
W_tild       	= Wch_curr+Wch_curr.*Options.Pen;
A               = [1, Wch_curr; 1,W_tild];
b               = [PIch_curr;0];
Linear_param    = (A^-1)*b;
PiC_i_curr(in_limit_Wch==0) = Linear_param(1)+Linear_param(2).*Wc_inner(in_limit_Wch==0);
Pic             = PiC_i_curr;
end

function [param_CUR] = Fit_F_CUR(map, InitialGuess, Constraints, Options)
res_CUR     = @(param)(map.CUR_vec-Options.F_CUR(param,map.Nc_vec./map.Nc_max_map));
% Choose Solver
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(Constraints.C_cur.lb,Constraints.C_cur.ub,0);
    param_CUR   = lsoptim(res_CUR,InitialGuess.param_CUR,Options.n_eval_init,A,b);
else
    param_CUR   = lsqnonlin(res_CUR,InitialGuess.param_CUR,Constraints.C_cur.lb,Constraints.C_cur.ub,Options.BaseSolverOpt);
end
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