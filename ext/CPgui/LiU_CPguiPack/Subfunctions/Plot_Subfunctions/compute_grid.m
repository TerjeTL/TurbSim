function [grid] = compute_grid(map,param,options)
%[grid] =  compute_grid(map,param,options)
%
% Creates the compressor mass flow grid area given the model parameters and
% the options
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
if ~isfield(options,'SurgeReverse')
    options.SurgeReverse = 0;
end
if options.SurgeReverse && isfield(param,'K0')
    options.points_in_SpL  = [100,options.points_in_SpL];
 	Nc_vec_n    = options.Nc_vec./map.Nc_max_map;
    Wch         = param.F_Wch(param.C_Wch,Nc_vec_n,map);
    WZSL        = param.F_WZS(param.C_Wzs,Nc_vec_n,map);
    grid.Wc     = zeros(sum(options.points_in_SpL),length(options.Nc_vec));
    grid.Nc     = repmat(options.Nc_vec',sum(options.points_in_SpL),1);
    grid.T01    = mean(map.T01);
    grid.p01    = mean(map.p01);
    for i = 1:length(Nc_vec_n)
        Wc_surgeRev   	= linspace(-param.K0.*map.Wc_max_map,0,options.points_in_SpL(1))';
        Wc_surge        = linspace(0,WZSL(i),options.points_in_SpL(2))';
        Wc_ellipse_1    = linspace(WZSL(i),0.98.*Wch(i),options.points_in_SpL(3)-ceil(0.2*options.points_in_SpL(3)))';
        Wc_ellipse_2    = linspace(0.98.*Wch(i),Wch(i),ceil(0.2*options.points_in_SpL(3)))';
        Wc_ellipse      = [Wc_ellipse_1;Wc_ellipse_2];
        Wc_choke        = Wch(i).*ones(options.points_in_SpL(4),1);
        grid.Wc(:,i)    = [Wc_surgeRev;Wc_surge;Wc_ellipse;Wc_choke];
    end
else
    Nc_vec_n    = options.Nc_vec./map.Nc_max_map;
    Wch         = param.F_Wch(param.C_Wch,Nc_vec_n,map);
    WZSL        = param.F_WZS(param.C_Wzs,Nc_vec_n,map);
    grid.Wc     = zeros(sum(options.points_in_SpL),length(options.Nc_vec));
    grid.Nc     = repmat(options.Nc_vec',sum(options.points_in_SpL),1);
    grid.T01    = mean(map.T01);
    grid.p01    = mean(map.p01);
    for i = 1:length(Nc_vec_n)
        W_surgStart     = (WZSL(i).*options.surgePercent).*Nc_vec_n(i);
        Wc_surge        = linspace(W_surgStart,WZSL(i),options.points_in_SpL(1))';
        Wc_ellipse_1    = linspace(WZSL(i),0.98.*Wch(i),options.points_in_SpL(2)-ceil(0.2*options.points_in_SpL(2)))';
        Wc_ellipse_2    = linspace(0.98.*Wch(i),Wch(i),ceil(0.2*options.points_in_SpL(2)))';
        Wc_ellipse      = [Wc_ellipse_1;Wc_ellipse_2];
        Wc_choke        = Wch(i).*ones(options.points_in_SpL(3),1);
        grid.Wc(:,i)    = [Wc_surge;Wc_ellipse;Wc_choke];
    end
end
end