function [a_p] = F_a_ASME(param,Nc_corr_n,map)
% [a_p] = F_a_ASME(param,Nc_corr_n,map)
%
% Computes A value given the parameters and the normalized corrected speed vector
%
%   param        	- Parameters of the subfunction [c1,c2]
%   Nc_corr_n       - Scalar or vector containing Speed values (normalized)
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
c1  = param(1:2);
% Compute the a_p for given N_c_corr (normalized)
a_p = map.Delta_h_max.*(c1(1).*Nc_corr_n-c1(2).*Nc_corr_n.^2)./map.Wc_max_map;
end