function [Wc,eta_c,T_out,Pc] = CompressorModel(param,p_in,p_out,T_in,omega_c)
% [Wc,eta_c,T_out,Pc] = CompressorModel(param,p_in,p_out,T_in,omega_c)
%
% Compressor model implementation in forward mode 
%
%   param         	- Struct containing the model parameters, from LiU CPgui
%   p_in            - Compressor inlet pressure [Pa]
%   p_out           - Compressor outlet pressure [Pa]
%   T_in            - Compressor inlet temperature [K]
%   omega_c       	- Compressor rotational speed [rad/s]
%
%   Wc              - Model calculated compressor mass flow [kg/s]
%   eta_c          	- Model calculated compressor efficiency [-]
%   T_out          	- Model calculated outlet temperature [K]
%   Pc              - Model calculated compressor consumed power [W]
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

%% Extract Model Parameters and preprocess the inputs
T_ref     	= param.T_ref;
p_ref     	= param.p_ref; 
Cp_air      = param.Cp_air; 
gamma_air   = param.gamma_air;

PiC         = p_out./p_in;
Nc_corr     = (omega_c.*(30./pi))./(sqrt(T_in./T_ref));  % corrected rad/s
Nc_corr_n   = Nc_corr./param.Nc_max_map;
Wc_corr     = zeros(length(PiC),1);
%% Mass flow model
% Compute base functions
Wch_SpL     = param.F_Wch(param.C_Wch,Nc_corr_n,param);
PIch_SpL    = param.F_PIch(param.C_PIch,Nc_corr_n,param);
WZSL_SpL    = param.F_WZS(param.C_Wzs,Nc_corr_n,param);
PIZSL_SpL   = param.F_PIZS(param.C_PIzs,Nc_corr_n,param);
CUR_SpL     = param.F_CUR(param.C_cur,Nc_corr_n);
% Calculate the ellipse value
Wc_corr_el 	= WZSL_SpL +(Wch_SpL-WZSL_SpL).*(1-((PiC-PIch_SpL)./(PIZSL_SpL-PIch_SpL)).^CUR_SpL).^(1./CUR_SpL);
% Calculate the linear surge part
Wc_corr_su 	= WZSL_SpL-(PiC-PIZSL_SpL)./(param.alpha_kPiW0);
% Determine the output part
Wc_corr(PiC<PIZSL_SpL)  = Wc_corr_el(PiC<PIZSL_SpL);
Wc_corr(PiC>=PIZSL_SpL) = Wc_corr_su(PiC>=PIZSL_SpL);
Wc_corr(PiC<PIch_SpL)   = Wch_SpL(PiC<PIch_SpL);
% Uncorrect the massflow
Wc          = Wc_corr.*(p_in./p_ref)./(sqrt(T_in./T_ref));
%% Efficiency model 
% Compute base functions
B_model    	= param.F_b(param.C_b,Nc_corr_n,param);
A_model    	= -param.F_a(param.C_a,Nc_corr_n,param);
K_fric   	= param.F_kloss(param.C_loss,Nc_corr,Wc_corr,param);
% Compute actual Enthalphy
D_H_act   = K_fric.*(B_model+A_model.*Wc_corr);
D_H_act(D_H_act<1e-4) = 1e-4;
% Compute inentropic Enthalphy
Delta_h_is  = Cp_air.*T_ref.*(PiC.^((gamma_air - 1)./gamma_air) - 1);
Delta_h_is(Delta_h_is<1e-4) = 1e-4;
% Compute efficiency with limitations
eta_c       = Delta_h_is./D_H_act;
eta_c(eta_c<1e-3)   = 1e-3;
eta_c(eta_c>1)      = 1;
%% Outlet temperature
T_out = T_in+(T_in./eta_c).*(PiC.^((gamma_air - 1)/gamma_air) - 1);
%% Consumed power
Pc  = Cp_air.*Wc.*(T_out-T_in);
