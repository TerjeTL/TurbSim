function [base] = compute_basefunctions(map,param,options)
%[base] = compute_basefunctions(map,param,options)
%
% Computes the base functions in a dense grid given the model parameters and
% the gridding options
%
%   map         	- Struct containing the map signals
%   param       	- Struct with the parameters
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
base.Nc_plot    = linspace(options.Nc_vec(1),options.Nc_vec(end),options.points_base_plot);
base.Wch_plot  	= param.F_Wch(param.C_Wch,base.Nc_plot./map.Nc_max_map,map);
base.PIch_plot 	= param.F_PIch(param.C_PIch,base.Nc_plot./map.Nc_max_map,map);
base.WZSL_plot 	= param.F_WZS(param.C_Wzs,base.Nc_plot./map.Nc_max_map,map);
base.PIZSL_plot = param.F_PIZS(param.C_PIzs,base.Nc_plot./map.Nc_max_map,map);
base.CUR_plot  	= param.F_CUR(param.C_cur,base.Nc_plot./map.Nc_max_map);
base.PI0_plot   = param.F_PI0(param.C_PIzs,base.Nc_plot./map.Nc_max_map,map);

base.Nc_vec     = options.Nc_vec;
base.Wch_vec    = param.F_Wch(param.C_Wch,base.Nc_vec./map.Nc_max_map,map);
base.PIch_vec   = param.F_PIch(param.C_PIch,base.Nc_vec./map.Nc_max_map,map);
base.WZSL_vec   = param.F_WZS(param.C_Wzs,base.Nc_vec./map.Nc_max_map,map);
base.PIZSL_vec  = param.F_PIZS(param.C_PIzs,base.Nc_vec./map.Nc_max_map,map);
base.CUR_vec    = param.F_CUR(param.C_cur,base.Nc_vec./map.Nc_max_map);
base.PI0_vec    = param.F_PI0(param.C_PIzs,base.Nc_vec./map.Nc_max_map,map);

if isfield(param,'C_a')
    base.B_plot  	= param.F_b(param.C_b,base.Nc_plot./map.Nc_max_map,map);
    base.A_plot  	= param.F_a(param.C_a,base.Nc_plot./map.Nc_max_map,map);
    base.B_vec      = param.F_b(param.C_b,base.Nc_vec./map.Nc_max_map,map);
    base.A_vec      = param.F_a(param.C_a,base.Nc_vec./map.Nc_max_map,map);
end
end