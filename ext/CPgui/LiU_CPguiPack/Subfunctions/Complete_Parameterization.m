function [paramT,errorsT] = Complete_Parameterization(map, InitialGuess, Constraints, Options)
%[paramT,errorsT] = Complete_Parameterization(map, InitialGuess, Constraints, Options)
%
% Solves a Total Least Squares problem to optimize the complete compressor model parameters given the Initial guess, the Constraints
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
Constraints         = cal_constraints(InitialGuess.Param_i,Constraints,Options);
%%  3D algorithm
Param_init          = [InitialGuess.Param_i;InitialGuess.Delta_i];
m                 	= length(InitialGuess.Delta_i);
%% Nested anonymous function
Residuals          	= @(param)(Residuals_3D_TLS(param,map,Options));
param_ub         	= Constraints.Vec.param_ub;
param_lb          	= Constraints.Vec.param_lb;
%% Choose Solver and solve
if Options.LSOptim
    [A,b]       = CreateLinearConstraints(param_lb,param_ub,m);
    Param_F     = lsoptim(Residuals,Param_init,Options.n_eval_3D,A,b);
else
    param_lb    = [param_lb;-inf.*ones(m,1)];
    param_ub    = [param_ub;inf.*ones(m,1)];
    Param_F     = lsqnonlin(Residuals,Param_init,param_lb,param_ub,Options.SolverOptions3D);
end
Delta_Wc     	= Param_F(Options.n_param+1:end);
Parameters      = Param_F(1:Options.n_param);
%% Return Parameter Struct
In              = Options.ParvecI;
if Options.n_Klossparam == 0
    paramT.vector   	= Parameters;
    paramT.delta_Wc 	= Delta_Wc;
    paramT.C_Wch      	= Parameters(1:In(1));
    paramT.C_PIch       = Parameters(In(1)+1:In(2));
    paramT.C_Wzs     	= Parameters(In(2)+1:In(3));
    paramT.C_PIzs       = Parameters(In(3)+1:In(4));
    paramT.C_cur      	= Parameters(In(4)+1:In(5));
    paramT.Gamma_PIcs  	= 0.5;
    paramT.C_s          = Parameters(In(5)+1:In(6));
    paramT.C_b        	= Parameters(In(6)+1:In(7));
    paramT.C_a        	= Parameters(In(7)+1:In(8));
    paramT.C_loss       = 0;
else
    paramT.vector   	= Parameters;
    paramT.delta_Wc 	= Delta_Wc;
    paramT.C_Wch      	= Parameters(1:In(1));
    paramT.C_PIch       = Parameters(In(1)+1:In(2));
    paramT.C_Wzs     	= Parameters(In(2)+1:In(3));
    paramT.C_PIzs       = Parameters(In(3)+1:In(4));
    paramT.C_cur      	= Parameters(In(4)+1:In(5));
    paramT.Gamma_PIcs 	= 0.5;
    paramT.C_s          = Parameters(In(5)+1:In(6));
    paramT.C_b        	= Parameters(In(6)+1:In(7));
    paramT.C_a        	= Parameters(In(7)+1:In(8));
    paramT.C_loss       = Parameters(In(8)+1:In(9));
end
paramT.F_Wch      	= Options.F_Wch;
paramT.F_PIch     	= Options.F_PIch;
paramT.F_WZS      	= Options.F_WZS;
paramT.F_PIZS      	= Options.F_PIZS;
paramT.F_CUR      	= Options.F_CUR;
paramT.F_PI0     	= Options.F_PI0;
paramT.F_a       	= Options.F_a;
paramT.F_b        	= Options.F_b;
paramT.F_kloss   	= Options.F_kloss;
paramT.T_ref       	= map.T_ref;
paramT.p_ref       	= map.p_ref;
paramT.Nc_max_map 	= map.Nc_max_map;
paramT.Wc_max_map  	= map.Wc_max_map;
paramT.PI_max_map 	= map.PI_max_map;
paramT.Delta_h_max	= map.Delta_h_max;
paramT.gamma_air   	= map.gamma_air;
paramT.rho1       	= map.rho1;
paramT.D_2        	= map.D_2;
paramT.Cp_air     	= map.Cp_air;
paramT.alpha_kPiW0  = 0.15*(map.PI_max_map/map.Wc_max_map);
% Compute Model Errors
errorsT 	= compute_modelErrors(map,paramT,Options);
end

%% Subfunctions
function [etac,Pic] = g_TLS(param,map,Options)
% July 2016 Xavier Llamas
% Version 1.0
%  + 
%% Extract Data
Wc_corr         = map.WcCorr_M;
PiC             = map.PiC_M;
n_parameters    = Options.n_param;
In              = Options.ParvecI;
Cp_air          = map.Cp_air;
T_1             = map.T_ref;
gamma_air       = map.gamma_air;
%% Preprocess data
Nc_corr_n       = map.Nc_vec./map.Nc_max_map;
Pic             = [];
etac            = [];
n_SpL_curr      = 0;
index_ini_curr  = n_parameters+1;
for i = 1: length(Nc_corr_n)
    %% Current speed line
    Wc_corr_curr    = Wc_corr(:,i);
    Wc_corr_curr    = Wc_corr_curr(isnan(Wc_corr_curr)==false);
    PiC_curr        = PiC(:,i);
    PiC_curr        = PiC_curr(isnan(PiC_curr)==false);
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
    % Get the linear parameters for enthalpy as function of compressor speed    
    a_Nc_corr	= Options.F_a(param(In(7)+1:In(8)),Nc_corr_n(i),map);
    b_Nc_corr  	= Options.F_b(param(In(6)+1:In(7)),Nc_corr_n(i),map);
    
    %% Initialize errors in variables as parameters
    Wc_inner        = Wc_corr_curr + param(index_ini_curr:index_ini_curr+n_SpL_curr-1);
    %% Check wether a map point is outside the range and compute the intersection in another manner for massflow model
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
    %% Linear model when Wc_inner > Wch_curr, as mush steep as possible to resemble a vertical choking line
    W_tild = Wch_SpL+Wch_SpL.*Options.Pen;
    A = [1, Wch_SpL; 1,W_tild];
    b = [PIch_SpL;0];
    Linear_param = (A^-1)*b;
    PiC_i_curr(in_limit_Wch==0) = Linear_param(1)+Linear_param(2).*Wc_inner(in_limit_Wch==0);   
    %% Check wether a map point is outside the range and compute the intersection in another manner for the efficiency model
    in_limit_Wc     = Wc_inner<b_Nc_corr./a_Nc_corr; % otherwise divide by zero ( normally far away from the Wc_inner values ( PiC_curr is already below 1
    in_limit_Pi     = PiC_i_curr>1; % delta h is not defined
    in_zone_pi      = in_limit_Wc.*in_limit_Pi;
    %% Return function evaluation
    %Wc_uncorr       = Wc_inner.*(p01_curr./map.pCref)./sqrt(T01_curr./map.TCref);
    K_fric          = Options.F_kloss(param(n_parameters),map.Nc_vec(i),Wc_inner,map);
    etac_i_curr_z	= (Cp_air*T_1.*(PiC_i_curr(in_zone_pi>0).^((gamma_air - 1)./gamma_air) - 1))./(K_fric(in_zone_pi>0).*(b_Nc_corr-a_Nc_corr.*Wc_inner(in_zone_pi>0)));    
   
    %% Concatenate vectors
    etac_i_curr       = zeros(length(PiC_curr),1);    
    etac_i_curr(in_zone_pi>0)   = etac_i_curr_z;
    etac_i_curr(in_limit_Wc==0) = 0;
    etac_i_curr(in_limit_Pi==0) = 0;
    etac                        = [etac;etac_i_curr];
    Pic                        	= [Pic;PiC_i_curr];
end
end

function [r] = Residuals_3D_TLS(param,map,Options)
% [r] = Residuals_3D_TLS(param,Data)

%
% March 2016 Xavier Llamas
% Version 1.0
%  + 
%% Compute intersection points.
[etac,Pic] = g_TLS(param,map,Options);
%% Compute squared distance  normalized with the maximun in the current SpL
% compute the maximums
Pic_vec    	= map.PiC_M(isnan(map.WcCorr_M(:))<1);
etac_vec    = map.etaC_M(isnan(map.etaC_M(:))<1);
r1          = (Pic-Pic_vec);
r2          = (etac-etac_vec);
r           = [r1;r2;param(Options.n_param+1:end)];
In        	= Options.ParvecI;

param_Wch   = param(1:In(1));
param_PIZSL = param(In(3)+1:In(4));
Wch_SpL     = Options.F_Wch(param_Wch,map.Nc_vec./map.Nc_max_map,map);
PIZSL_SpL   = Options.F_PIZS(param_PIZSL,map.Nc_vec./map.Nc_max_map,map);
eta_max     = max(map.etaC_M)';
size_M      = size(map.PiC_M);
W_norm      = repmat(Wch_SpL',size_M(1),1);
W_norm      = W_norm(isnan(map.WcCorr_M(:))<1);
Pi_norm     = repmat(PIZSL_SpL',size_M(1),1);
Pi_norm     = Pi_norm(isnan(map.WcCorr_M(:))<1);
eta_norm    = repmat(eta_max',size_M(1),1);
eta_norm    = eta_norm(isnan(map.WcCorr_M(:))<1);
Norm = [Pi_norm;eta_norm;W_norm];
r = r./Norm;
end

function [errors] = compute_modelErrors(map,param,Options)
% get Intersection Points    
[inter.etaC,inter.PiC]	= g_TLS([param.vector;param.delta_Wc],map,Options);
inter.Wc                = map.WcCorr_M(isnan(map.WcCorr_M(:))<1)+param.delta_Wc;
% Compute absolute relative errors
PiC_vec                 = map.PiC_M(isnan(map.PiC_M(:))<1);
errors.PiC_rel          = abs((PiC_vec-inter.PiC)./(1/length(PiC_vec)*sum(PiC_vec)))*100;
errors.PiC_mean_rel     = 1/length(PiC_vec)*sum(errors.PiC_rel);

Wc_vec                  = map.WcCorr_M(isnan(map.WcCorr_M(:))<1);
errors.Wc_rel           = abs((Wc_vec-inter.Wc)./(1/length(Wc_vec)*sum(Wc_vec)))*100;
errors.Wc_mean_rel      = 1/length(Wc_vec)*sum(errors.Wc_rel);

etac_vec                = map.etaC_M(isnan(map.etaC_M(:))<1);
errors.eta_rel          = abs((etac_vec-inter.etaC)./(1/length(etac_vec)*sum(etac_vec)))*100;
errors.eta_mean_rel     = 1/length(etac_vec)*sum(errors.eta_rel);
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
%% CUR bound constaints
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
    if isempty(Constraints.C_b.lb)  
        Constraints.LBcurr.C_b 	= -inf.*ones(length(In(6)+1:In(7)),1);
    else
        Constraints.LBcurr.C_b 	= Constraints.C_b.lb;
    end
    if isempty(Constraints.C_b.ub) 
        Constraints.UBcurr.C_b 	= inf.*ones(length(In(6)+1:In(7)),1);
    else
        Constraints.UBcurr.C_b    = Constraints.C_b.ub;
    end
    %% A bound constaints
    if isempty(Constraints.C_a.lb)  
        Constraints.LBcurr.C_a 	= -inf.*ones(length(In(7)+1:In(7)),1);
    else
        Constraints.LBcurr.C_a    = Constraints.C_a.lb;
    end
    if isempty(Constraints.C_a.ub) 
        Constraints.UBcurr.C_a 	= inf.*ones(length(In(7)+1:In(8)),1);
    else
        Constraints.UBcurr.C_a    = Constraints.C_a.ub;
    end

    %% K_fric bound constaints
    if isempty(Constraints.C_loss.lb)   
        Constraints.LBcurr.C_loss 	= 0;
    else
        Constraints.LBcurr.C_loss   = Constraints.C_loss.lb;
    end
    if isempty(Constraints.C_loss.ub)   
        Constraints.UBcurr.C_loss 	= inf;
    else
        Constraints.UBcurr.C_loss   = Constraints.C_loss.ub;
    end  
 %% Concatenate Vectors 
 Constraints.Vec.param_lb(In(6)+1:In(7))  = Constraints.LBcurr.C_b;
 Constraints.Vec.param_lb(In(7)+1:In(8))  = Constraints.LBcurr.C_a;
 Constraints.Vec.param_lb(In(8)+1:In(9))  = Constraints.LBcurr.C_loss;
  
 Constraints.Vec.param_ub(In(6)+1:In(7))  = Constraints.UBcurr.C_b;
 Constraints.Vec.param_ub(In(7)+1:In(8))  = Constraints.UBcurr.C_a;
 Constraints.Vec.param_ub(In(8)+1:In(9))  = Constraints.UBcurr.C_loss;
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