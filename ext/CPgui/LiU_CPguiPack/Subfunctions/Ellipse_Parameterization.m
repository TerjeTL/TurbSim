function [paramEll,errorsEll] = Ellipse_Parameterization(map, InitialGuess, Constraints, Options)
%[paramEll,errorsEll] = Ellipse_Parameterization(map, InitialGuess, Constraints, Options)
%
% Solves a Total Least Squares problem to optimize the Ellipse parameters given the Initial guess, the Constraints
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
%% Compute Constraints given thresholds and current parameters
Constraints   	= cal_constraints(InitialGuess.Param_i,Constraints,Options);
%% Complete optimization of the Ellipse mass flow model
curr_param_num  = length(InitialGuess.Param_i);
Param_init      = [InitialGuess.Param_i;InitialGuess.Delta_i];
m               = length(InitialGuess.Delta_i);
% Nested residuals function
Residuals       = @(param)(Residuals_2D_TLS(param,map,Options));
param_ub        = Constraints.Vec.param_ub;
param_lb        = Constraints.Vec.param_lb;
%% Choose Solver and optimize
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(param_lb,param_ub,m);
    % Solve nonlinear Least squares with LSOPTIM
    Param_mf    = lsoptim(Residuals,Param_init,Options.n_eval_2D,A,b);
else
    param_lb    = [param_lb;-inf.*ones(m,1)];
    param_ub    = [param_ub;inf.*ones(m,1)];
    Param_mf    = lsqnonlin(Residuals,Param_init,param_lb,param_ub,Options.SolverOptions);
end
% Extract Parameters
Delta_Wc    = Param_mf(curr_param_num+1:end);
Param_mf    = Param_mf(1:curr_param_num);
In       	= Options.ParvecI;
paramEll.vector   	= Param_mf;
paramEll.delta_Wc 	= Delta_Wc;
paramEll.C_Wch   	= Param_mf(1:In(1));
paramEll.C_PIch   	= Param_mf(In(1)+1:In(2));
paramEll.C_Wzs     	= Param_mf(In(2)+1:In(3));
paramEll.C_PIzs   	= Param_mf(In(3)+1:In(4));
paramEll.C_cur   	= Param_mf(In(4)+1:In(5));
paramEll.Gamma_PIcs	= 0.5;
paramEll.C_s        = Param_mf(In(5)+1:In(6));
paramEll.F_Wch   	= Options.F_Wch;
paramEll.F_PIch  	= Options.F_PIch;
paramEll.F_WZS    	= Options.F_WZS;
paramEll.F_PIZS 	= Options.F_PIZS;
paramEll.F_CUR    	= Options.F_CUR;
paramEll.F_PI0    	= Options.F_PI0;
paramEll.T_ref     	= map.T_ref;
paramEll.p_ref    	= map.p_ref;
paramEll.Nc_max_map	= map.Nc_max_map;
paramEll.Wc_max_map	= map.Wc_max_map;
paramEll.PI_max_map	= map.PI_max_map;
paramEll.Delta_h_max= map.Delta_h_max;
paramEll.gamma_air 	= map.gamma_air;
paramEll.rho1      	= map.rho1;
paramEll.D_2     	= map.D_2;
paramEll.Cp_air  	= map.Cp_air;
paramEll.alpha_kPiW0= 0.15*(map.PI_max_map/map.Wc_max_map);
% Compute Model Errors
errorsEll 	= compute_modelErrors(map,paramEll,Options);
end

%% Subfunctions
function [r] = Residuals_2D_TLS(param,map,Options)
%% Compute intersection points.
Pic         = f_TLS(param,map,Options);
%% Compute squared distance  normalized with the maximun in the current SpL
% compute the maximums
In       	= Options.ParvecI;
Pic_vec    	= map.PiC_M(isnan(map.WcCorr_M(:))<1);
r           = (Pic-Pic_vec);
r           = [r;param(In(6)+1:end)];

param_Wch   = param(1:In(1));
param_PIZSL = param(In(3)+1:In(4));
Wch_SpL     = Options.F_Wch(param_Wch,map.Nc_vec./map.Nc_max_map,map);
PIZSL_SpL   = Options.F_PIZS(param_PIZSL,map.Nc_vec./map.Nc_max_map,map);

size_M      = size(map.WcCorr_M);
W_norm      = repmat(Wch_SpL',size_M(1),1);
W_norm      = W_norm(isnan(map.WcCorr_M(:))<1);
Pi_norm     = repmat(PIZSL_SpL',size_M(1),1);
Pi_norm     = Pi_norm(isnan(map.WcCorr_M(:))<1);

Norm        = [Pi_norm;W_norm];
r           = r./Norm;
end

function [PiC] = f_TLS(param,map,Options)

% July 2016 Xavier Llamas
% Version 1.0
%  + 
%% Extract Data
Wc_corr         = map.WcCorr_M;
%PiC             = map.PiC_M;
In              = Options.ParvecI;
n_parameters    = In(6);
%% Preprocess data
Nc_corr_n       = map.Nc_vec_n;
PiC            	= zeros(map.n_measured,1);
n_SpL_curr      = 0;
index_ini_curr  = n_parameters+1;
index_PiC_curr  = 1;
for i = 1: length(Nc_corr_n)
    %% Current speed line
    Wc_corr_curr    = Wc_corr(:,i);
    Wc_corr_curr    = Wc_corr_curr(isnan(Wc_corr_curr)==false);
    index_ini_curr  = index_ini_curr +n_SpL_curr;
    n_SpL_curr      = length(Wc_corr_curr);
    %% Compute -> Wch_SpL PIch_SpL WZSL_SpL PIZSL_SpL CUR_SpL given Nc_corr
    param_Wch   = param(1:In(1));
    param_PIch  = param(In(1)+1:In(2));
    param_WZSL  = param(In(2)+1:In(3));
    param_PIZSL = param(In(3)+1:In(4));
    param_CUR   = param(In(4)+1:In(5));
    C_s     = param(In(5)+1:In(6));
    % With this the speed line is fully determined.
    Wch_SpL     = Options.F_Wch(param_Wch,Nc_corr_n(i),map);
    PIch_SpL    = Options.F_PIch(param_PIch,Nc_corr_n(i),map);
    WZSL_SpL    = Options.F_WZS(param_WZSL,Nc_corr_n(i),map);
    PIZSL_SpL   = Options.F_PIZS(param_PIZSL,Nc_corr_n(i),map);
    CUR_SpL     = Options.F_CUR(param_CUR,Nc_corr_n(i));
    PI0_SpL     = Options.F_PI0(param_PIZSL,Nc_corr_n(i),map);
    %% Initialize errors in variables as parameters
    Wc_inner   = Wc_corr_curr + param(index_ini_curr:index_ini_curr+n_SpL_curr-1);
    %% Check wether a map point is outside the range and compute the intersection in another manner
    in_limit_Wch    = Wc_inner<=Wch_SpL;
    in_limit_Wzs    = Wc_inner>WZSL_SpL;
    in_limit_surge  = Wc_inner<=WZSL_SpL;
    in_zone         = in_limit_Wch.*in_limit_Wzs;
    %% Return ellipse evaluation
    PiC_i_curr_e     = PIch_SpL+(PIZSL_SpL-PIch_SpL).*(1-((Wc_inner(in_zone>0)-WZSL_SpL)./(Wch_SpL-WZSL_SpL)).^CUR_SpL).^(1./CUR_SpL);
    %% Return surge polynomial evaluation
    W_tilde          = WZSL_SpL.*(1-(1-(Wc_inner(in_limit_surge>0)./WZSL_SpL)).^C_s).^(1./C_s);
    PiC_i_curr_p     = PI0_SpL +3.*((PIZSL_SpL-PI0_SpL)./WZSL_SpL.^2).*W_tilde.^2-2.*((PIZSL_SpL-PI0_SpL)./WZSL_SpL.^3).*W_tilde.^3;
    %% Concatenate vectors
    PiC_i_curr       = zeros(length(Wc_corr_curr),1); 
    PiC_i_curr(in_zone>0)        = PiC_i_curr_e;
    PiC_i_curr(in_limit_surge>0) = PiC_i_curr_p;
    %% Linear model when Wc_inner > Wch_curr, as much steep as possible to resemble a vertical choking line
    W_tild = Wch_SpL+Wch_SpL.*Options.Pen;
    A = [1, Wch_SpL; 1,W_tild];
    b = [PIch_SpL;0];
    Linear_param = (A^-1)*b;
    PiC_i_curr(in_limit_Wch==0)	= Linear_param(1)+Linear_param(2).*Wc_inner(in_limit_Wch==0);
    %Pic                          = [Pic;PiC_i_curr];
    PiC(index_PiC_curr:index_PiC_curr+n_SpL_curr-1)     	= PiC_i_curr;
    index_PiC_curr = index_PiC_curr+n_SpL_curr;
end
end

function [errors] = compute_modelErrors(map,param,Options)
% get Intersection Points    
inter.PiC           = f_TLS([param.vector;param.delta_Wc],map,Options);
inter.Wc            = map.WcCorr_M(isnan(map.WcCorr_M(:))<1)+param.delta_Wc;
% Compute absolute relative errors
PiC_vec             = map.PiC_M(isnan(map.PiC_M(:))<1);
errors.PiC_rel      = abs((PiC_vec-inter.PiC)./(1/length(PiC_vec)*sum(PiC_vec)))*100;
errors.PiC_mean_rel = 1/length(PiC_vec)*sum(errors.PiC_rel);

Wc_vec              = map.WcCorr_M(isnan(map.WcCorr_M(:))<1);
errors.Wc_rel       = abs((Wc_vec-inter.Wc)./(1/length(Wc_vec)*sum(Wc_vec)))*100;
errors.Wc_mean_rel  = 1/length(Wc_vec)*sum(errors.Wc_rel);
end

function [Constraints] = cal_constraints(Param_init,Constraints,Options)
n_param                         = length(Param_init);
Constraints.Vec.param_lb        = -inf.*ones(n_param,1);
Constraints.Vec.param_ub        = inf.*ones(n_param,1);
In                              = Options.ParvecI;
%% W_ch bound constraints, between the hard limits and the calculated surrounding
% Lower Bound
if isempty(Constraints.C_Wch.lb)
    Constraints.LBcurr.C_Wch      = Param_init(1:In(1))-abs(Param_init(1:In(1)).*Constraints.C_Wch.Th_lb);
else
    param_lb_2                  = Param_init(1:In(1))-abs(Param_init(1:In(1)).*Constraints.C_Wch.Th_lb);
    param_lb_1                  = Constraints.C_Wch.lb;
    Constraints.LBcurr.C_Wch      = max(param_lb_2,param_lb_1);   
end
% Upper Bound
if isempty(Constraints.C_Wch.ub)
    Constraints.UBcurr.C_Wch      = Param_init(1:In(1))+abs(Param_init(1:In(1)).*Constraints.C_Wch.Th_ub);
else
    param_ub_2                  = Param_init(1:In(1))+abs(Param_init(1:In(1)).*Constraints.C_Wch.Th_ub);
    param_ub_1                  = Constraints.C_Wch.ub;
    Constraints.UBcurr.C_Wch      = min(param_ub_2,param_ub_1);        
end
%% Pi_ch bound Constraints 
% Lower Bound
if isempty(Constraints.C_PIch.lb) 
    Constraints.LBcurr.C_PIch     = Param_init(In(1)+1:In(2))-abs(Param_init(In(1)+1:In(2)).*Constraints.C_PIch.Th_lb);
else
    param_lb_2                  = Param_init(In(1)+1:In(2))-abs(Param_init(In(1)+1:In(2)).*Constraints.C_PIch.Th_lb);
    param_lb_1                  = Constraints.C_PIch.lb;
    Constraints.LBcurr.C_PIch     = max(param_lb_2,param_lb_1);      
end  
% Upper Bound
if isempty(Constraints.C_PIch.ub)  
    Constraints.UBcurr.C_PIch     = Param_init(In(1)+1:In(2))+abs(Param_init(In(1)+1:In(2)).*Constraints.C_PIch.Th_ub);
else
    param_ub_2                  = Param_init(In(1)+1:In(2))+abs(Param_init(In(1)+1:In(2)).*Constraints.C_PIch.Th_ub);
    param_ub_1                  = Constraints.C_PIch.ub;
    Constraints.UBcurr.C_PIch     = min(param_ub_2,param_ub_1);    
end
%% W_ZSL bound Constraints
% Lower Bound
if isempty(Constraints.C_Wzs.lb) 
    Constraints.LBcurr.C_Wzs     = Param_init(In(2)+1:In(3))-abs(Param_init(In(2)+1:In(3)).*Constraints.C_Wzs.Th_lb);
else
    param_lb_2                  = Param_init(In(2)+1:In(3))-abs(Param_init(In(2)+1:In(3)).*Constraints.C_Wzs.Th_lb);
    param_lb_1                  = Constraints.C_Wzs.lb;
    Constraints.LBcurr.C_Wzs     = max(param_lb_2,param_lb_1);
end
% Upper Bound
if isempty(Constraints.C_Wzs.ub)
    Constraints.UBcurr.C_Wzs     = Param_init(In(2)+1:In(3))+abs(Param_init(In(2)+1:In(3)).*Constraints.C_Wzs.Th_ub);
else
    param_ub_2                  = Param_init(In(2)+1:In(3))+abs(Param_init(In(2)+1:In(3)).*Constraints.C_Wzs.Th_ub);
    param_ub_1                  = Constraints.C_Wzs.ub;
    Constraints.UBcurr.C_Wzs     = min(param_ub_2,param_ub_1);
   
end
%% Pi_ZSL bound Constraints
% Lower Bound
if isempty(Constraints.C_PIzs.lb) 
    Constraints.LBcurr.C_PIzs    = Param_init(In(3)+1:In(4))-abs(Param_init(In(3)+1:In(4)).*Constraints.C_PIzs.Th_lb);
else
    param_lb_2                  = Param_init(In(3)+1:In(4))-abs(Param_init(In(3)+1:In(4)).*Constraints.C_PIzs.Th_lb);
    param_lb_1                  = Constraints.C_PIzs.lb;
    Constraints.LBcurr.C_PIzs    = max(param_lb_2,param_lb_1);   
end
% Upper Bound
if isempty(Constraints.C_PIzs.ub)
    Constraints.UBcurr.C_PIzs    = Param_init(In(3)+1:In(4))+abs(Param_init(In(3)+1:In(4)).*Constraints.C_PIzs.Th_ub);
else
    param_ub_2                  = Param_init(In(3)+1:In(4))+abs(Param_init(In(3)+1:In(4)).*Constraints.C_PIzs.Th_ub);
    param_ub_1                  = Constraints.C_PIzs.ub;
    Constraints.UBcurr.C_PIzs    = min(param_ub_2,param_ub_1);
end
%% C_cur bound constaints
% Lower Bound
if isempty(Constraints.C_cur.lb)  
    Constraints.LBcurr.C_cur      = [0,0,0];
else
    Constraints.LBcurr.C_cur      = Constraints.C_cur.lb;
end
% Upper Bound
if isempty(Constraints.C_cur.ub)   
    Constraints.UBcurr.C_cur      = [2,2,8];
else
    Constraints.UBcurr.C_cur      = Constraints.C_cur.ub;
end

%% C_s bound constaints
% Lower Bound
if isempty(Constraints.C_s.lb) 
    Constraints.LBcurr.C_s  = 0.75;
else
    Constraints.LBcurr.C_s  = Constraints.C_s.lb;
end
% Upper Bound
if isempty(Constraints.C_s.ub) 
    Constraints.UBcurr.C_s  = 2.5;
else
    Constraints.UBcurr.C_s  = Constraints.C_s.ub;
end

%% Concatenate Vectors and return
Constraints.Vec.param_lb(1:In(1))           = Constraints.LBcurr.C_Wch;
Constraints.Vec.param_lb(In(1)+1:In(2))     = Constraints.LBcurr.C_PIch;
Constraints.Vec.param_lb(In(2)+1:In(3))     = Constraints.LBcurr.C_Wzs;
Constraints.Vec.param_lb(In(3)+1:In(4))     = Constraints.LBcurr.C_PIzs;
Constraints.Vec.param_lb(In(4)+1:In(5))     = Constraints.LBcurr.C_cur; 
Constraints.Vec.param_lb(In(5)+1:In(6))     = Constraints.LBcurr.C_s;

Constraints.Vec.param_ub(1:In(1))           = Constraints.UBcurr.C_Wch;
Constraints.Vec.param_ub(In(1)+1:In(2))     = Constraints.UBcurr.C_PIch;
Constraints.Vec.param_ub(In(2)+1:In(3))     = Constraints.UBcurr.C_Wzs;
Constraints.Vec.param_ub(In(3)+1:In(4))     = Constraints.UBcurr.C_PIzs;
Constraints.Vec.param_ub(In(4)+1:In(5))  	= Constraints.UBcurr.C_cur;
Constraints.Vec.param_ub(In(5)+1:In(6)) 	= Constraints.UBcurr.C_s;

%% Include Efficiency constraints if it is the case
if n_param > In(6)
  %% B bound constaints
    if isempty(Constraints.B.lb)  
        Constraints.LBcurr.B 	= -inf.*ones(length(In(6)+1:In(7)),1);
    else
        Constraints.LBcurr.B 	= Constraints.B.lb;
    end
    if isempty(Constraints.B.ub) 
        Constraints.UBcurr.B 	= inf.*ones(length(In(6)+1:In(7)),1);
    else
        Constraints.UBcurr.B    = Constraints.B.ub;
    end
    %% A bound constaints
    if isempty(Constraints.A.lb)  
        Constraints.LBcurr.A 	= -inf.*ones(length(In(7)+1:In(7)),1);
    else
        Constraints.LBcurr.A    = Constraints.A.lb;
    end
    if isempty(Constraints.A.ub) 
        Constraints.UBcurr.A 	= inf.*ones(length(In(7)+1:In(8)),1);
    else
        Constraints.UBcurr.A    = Constraints.A.ub;
    end

    %% K_fric bound constaints
    if isempty(Constraints.k_loss.lb)   
        Constraints.LBcurr.k_loss 	= 0;
    else
        Constraints.LBcurr.k_loss   = Constraints.k_loss.lb;
    end
    if isempty(Constraints.k_loss.ub)   
        Constraints.UBcurr.k_loss 	= inf;
    else
        Constraints.UBcurr.k_loss   = Constraints.k_loss.ub;
    end  
 %% Concatenate Vectors 
 Constraints.Vec.param_lb(In(6)+1:In(7))  = Constraints.LBcurr.B;
 Constraints.Vec.param_lb(In(7)+1:In(8))  = Constraints.LBcurr.A;
 Constraints.Vec.param_lb(In(8)+1:In(9))  = Constraints.LBcurr.k_loss;
  
 Constraints.Vec.param_ub(In(6)+1:In(7))  = Constraints.UBcurr.B;
 Constraints.Vec.param_ub(In(7)+1:In(8))  = Constraints.UBcurr.A;
 Constraints.Vec.param_ub(In(8)+1:In(9))  = Constraints.UBcurr.k_loss;
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