function [WZSL_SpL] = F_WZSL(param,Nc_corr_n,map)
% [WZSL_SpL] = F_WZSL_SpL(param,Nc_corr_n)
%
% Computes the Zero Slope Line Mass Flow given the parameters and normalized corrected speed vector
%
%   param           - Parameters of the subfunction [c1,c2]
%   Nc_corr_n       - Scalar or vector containing Speed values
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
%% Extract parameters
c1          = param(1:2);
% Compute the zero slope line mass flow
WZSL_SpL=map.Wc_max_map.*(c1(1).*Nc_corr_n.^c1(2));
end