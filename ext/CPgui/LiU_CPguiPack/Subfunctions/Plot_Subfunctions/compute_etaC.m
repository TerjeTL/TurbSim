function [grid] = compute_etaC(map,param,grid,options)
%[grid] = compute_etaC(map,param,grid,options)
%
% Computes the efficiency given the gridded compressor mass flow and
% pressure ratio areas and the model parameters
%
%   map         	- Struct containing the map signals
%   param       	- Struct with the parameters
%   grid            - Struct with the gridded extrapolation of the
%                     compressor performance variables
%   options         - Struct containing the grid options
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
%%
B           = param.F_b(param.C_b,grid.Nc./map.Nc_max_map,map);
A           = -param.F_a(param.C_a,grid.Nc./map.Nc_max_map,map);
%Wc_uncorr 	= grid.Wc.*(grid.p01./map.pCref)./sqrt(grid.T01./map.TCref);
K_fric   	= param.F_kloss(param.C_loss,grid.Nc,grid.Wc,map);
grid.Dh     = max(K_fric.*(B+A.*grid.Wc),0);
grid.Dhis   = max(map.Cp_air.*grid.T01.*(grid.PiC.^((map.gamma_air - 1)./map.gamma_air) - 1),0);
grid.etaC   = max(grid.Dhis./grid.Dh ,0);

if isfield(options,'lambda')
    rho             = grid.p01./map.R_air./grid.T01;
    grid.Uc         = grid.Nc.*sqrt(grid.T01./map.TCref).*pi./30.*map.D_2./2;
    grid.Mach       = grid.Uc./sqrt(map.gamma_air.*map.R_air.*grid.T01);  
    grid.lambda     = grid.Dh./grid.Uc.^2;  
    grid.Phi        = grid.Wc./rho./map.D_2^2./grid.Uc;
    grid.T02        = grid.T01.*(1+(grid.PiC.^((map.gamma_air-1)./map.gamma_air)-1)./grid.etaC);
    grid.eta_poly   = ((map.gamma_air-1)./map.gamma_air).*log(grid.PiC)./log(grid.T02./grid.T01);
    if isfield(options,'Ns')
        % Compute specific speed and eta_poly
        grid.omegaC     = (grid.Nc.*sqrt(grid.T01./map.TCref)).*pi./30;
        grid.rhoe       = grid.PiC.*grid.p01./map.R_air./grid.T02;
        %grid.Ns         = grid.omegaC.*sqrt(Wc_uncorr./grid.rhoe).*grid.Dhis.^(-3/4);
    end
end
end