function [PI0_SpL] = F_PI0(param,Nc_corr_n,map)
% [PI0_SpL] = F_PI0(param,Nc_corr_n,map)
%
% Computes the pressure ratio at zero flow given the parameters and normalized corrected speed vector
%
%   param           - Parameters of the subfunction, requires PIZSL parameters
%   Nc_corr_n     	- Scalar or vector containing Speed values (normalized)
%   map             - Struct containing the map signals
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
%% Get Tau from the selected value in the map or give it a default value
if isfield(map,'Gamma_PIcs')
    Gamma_PIcs    = param.Gamma_PIcs;
else
    Gamma_PIcs    = 0.5;
end
PIZSL_SpL   = F_PIZSL(param,Nc_corr_n,map);
% Compute the zero flow pressure ratio
PI0_SpL     = PIZSL_SpL-Gamma_PIcs.*(PIZSL_SpL-1);
end