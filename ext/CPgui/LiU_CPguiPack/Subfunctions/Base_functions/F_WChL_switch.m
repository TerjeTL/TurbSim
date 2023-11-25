function [Wch_SpL] = F_WChL_switch(param,Nc_corr_n,map)
% [Wch_SpL] = F_WChL_switch(param,Nc_corr_n,map)
%
% Computes the choke line mass flow given the parameters and normalized corrected speed vector
%
%   param        	- Parameters of the subfunction [c1,c2,c3,c4,Nc_shift]
%   Nc_corr_n    	- Scalar or vector containing Speed values (normalized)
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
c1          = param(1:3);
c2          = param(4);
Nc_shift    = param(5);
param_fixed = c1(1)+c1(2).*Nc_shift.^c1(3)-c2(1).*Nc_shift;
% Compute choke mass flow shifting the equation when necessary.
Wch_SpL=map.Wc_max_map.*((c1(1)+c1(2).*Nc_corr_n.^c1(3)).*(Nc_corr_n<= Nc_shift)+(param_fixed+c2(1).*Nc_corr_n).*(Nc_corr_n> Nc_shift));
end