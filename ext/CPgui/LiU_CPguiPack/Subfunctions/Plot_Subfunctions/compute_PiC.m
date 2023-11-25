function [grid] = compute_PiC(map,param,grid,options)
%[grid] = compute_PiC(map,param,grid,options)
%
% Computes the pressure ratio given the gridded compressor mass flow area
% and the model parameters
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
if ~isfield(options,'SurgeReverse')
    options.SurgeReverse = 0;
elseif ~isfield(options,'linearSurge')
    options.linearSurge  = 0;
end
if options.SurgeReverse && isfield(param,'K0')
 	Wch         = param.F_Wch(param.C_Wch,grid.Nc./map.Nc_max_map,map);
    PIch        = param.F_PIch(param.C_PIch,grid.Nc./map.Nc_max_map,map);
    WZSL        = param.F_WZS(param.C_Wzs,grid.Nc./map.Nc_max_map,map);
    PIZSL       = param.F_PIZS(param.C_PIzs,grid.Nc./map.Nc_max_map,map);
    CUR         = param.F_CUR(param.C_cur,grid.Nc./map.Nc_max_map);
    PI0         = param.F_PI0(param.C_PIzs,grid.Nc./map.Nc_max_map,map);
    C_s     = param.C_s;
    Kt          = param.Kt;
    K0          = param.K0;
    % Initialize PiC grid
    grid.PiC    = zeros(size(grid.Wc));
    % In Reverse flow Region
	index_Rev   = grid.Wc<=0;
    grid.PiC(index_Rev>0) = -1+PI0(index_Rev>0) +(1-(-grid.Wc(index_Rev>0)./(K0.*map.Wc_max_map)).^2).^(-1/Kt);
    % In surge Region
    index_ZSL   = (0<grid.Wc).*(grid.Wc<WZSL);
    W_tilde     = WZSL(index_ZSL>0).*(1-(1-(grid.Wc(index_ZSL>0)./WZSL(index_ZSL>0))).^C_s).^(1./C_s);
    grid.PiC(index_ZSL>0) = PI0(index_ZSL>0) +3.*((PIZSL(index_ZSL>0)-PI0(index_ZSL>0))./WZSL(index_ZSL>0).^2).*W_tilde.^2-2.*((PIZSL(index_ZSL>0)-PI0(index_ZSL>0))./WZSL(index_ZSL>0).^3).*W_tilde.^3;
    % In Ellipse Region
    index_ellipse = (grid.Wc>=WZSL).*(grid.Wc<Wch);
    grid.PiC(index_ellipse>0) = (((1-(((grid.Wc(index_ellipse>0)-WZSL(index_ellipse>0))./(Wch(index_ellipse>0)-WZSL(index_ellipse>0))).^CUR(index_ellipse>0)))).^(1./CUR(index_ellipse>0))).*(PIZSL(index_ellipse>0)-PIch(index_ellipse>0))+PIch(index_ellipse>0);
    % In Choke Region
    %index_choke = grid.Wc==Wch;
    n_choke = options.points_in_SpL(3);
    for i = 1:length(options.Nc_vec)
        grid.PiC(end-n_choke:end,i) = linspace(PIch(end,i),options.lower_PiC,n_choke+1)';
    end
elseif options.linearSurge && isfield(param,'alpha_kPiW0')
	Wch         = param.F_Wch(param.C_Wch,grid.Nc./map.Nc_max_map,map);
    PIch        = param.F_PIch(param.C_PIch,grid.Nc./map.Nc_max_map,map);
    WZSL        = param.F_WZS(param.C_Wzs,grid.Nc./map.Nc_max_map,map);
    PIZSL       = param.F_PIZS(param.C_PIzs,grid.Nc./map.Nc_max_map,map);
    CUR         = param.F_CUR(param.C_cur,grid.Nc./map.Nc_max_map);
    
    C_s     = param.C_s;
    % Initialize PiC grid
    grid.PiC    = zeros(size(grid.Wc));
    % In surge Region
    index_ZSL   = grid.Wc<WZSL;
    % Get an equal slope for all SpLs
    %alpha = -0.7*map.PI_max_map*(param.kPiW0-1)/(0.5*map.Wc_max_map);
    alpha = -param.alpha_kPiW0;
    grid.PiC(index_ZSL>0)         = PIZSL(index_ZSL>0)+alpha.*(grid.Wc(index_ZSL>0)-WZSL(index_ZSL>0));
    % In Ellipse Region
    index_ellipse = (grid.Wc>=WZSL).*(grid.Wc<Wch);
    grid.PiC(index_ellipse>0) = (((1-(((grid.Wc(index_ellipse>0)-WZSL(index_ellipse>0))./(Wch(index_ellipse>0)-WZSL(index_ellipse>0))).^CUR(index_ellipse>0)))).^(1./CUR(index_ellipse>0))).*(PIZSL(index_ellipse>0)-PIch(index_ellipse>0))+PIch(index_ellipse>0);
    % In Choke Region
    %index_choke = grid.Wc==Wch;
    n_choke = options.points_in_SpL(3);
    for i = 1:length(options.Nc_vec)
        grid.PiC(end-n_choke:end,i) = linspace(PIch(end,i),options.lower_PiC,n_choke+1)';
    end
else
    Wch         = param.F_Wch(param.C_Wch,grid.Nc./map.Nc_max_map,map);
    PIch        = param.F_PIch(param.C_PIch,grid.Nc./map.Nc_max_map,map);
    WZSL        = param.F_WZS(param.C_Wzs,grid.Nc./map.Nc_max_map,map);
    PIZSL       = param.F_PIZS(param.C_PIzs,grid.Nc./map.Nc_max_map,map);
    CUR         = param.F_CUR(param.C_cur,grid.Nc./map.Nc_max_map);
    PI0         = param.F_PI0(param.C_PIzs,grid.Nc./map.Nc_max_map,map);
    C_s     = param.C_s;
    % Initialize PiC grid
    grid.PiC    = zeros(size(grid.Wc));
    % In surge Region
    index_ZSL   = grid.Wc<WZSL;
    W_tilde     = WZSL(index_ZSL>0).*(1-(1-(grid.Wc(index_ZSL>0)./WZSL(index_ZSL>0))).^C_s).^(1./C_s);
    grid.PiC(index_ZSL>0) = PI0(index_ZSL>0) +3.*((PIZSL(index_ZSL>0)-PI0(index_ZSL>0))./WZSL(index_ZSL>0).^2).*W_tilde.^2-2.*((PIZSL(index_ZSL>0)-PI0(index_ZSL>0))./WZSL(index_ZSL>0).^3).*W_tilde.^3;
    % In Ellipse Region
    index_ellipse = (grid.Wc>=WZSL).*(grid.Wc<Wch);
    grid.PiC(index_ellipse>0) = (((1-(((grid.Wc(index_ellipse>0)-WZSL(index_ellipse>0))./(Wch(index_ellipse>0)-WZSL(index_ellipse>0))).^CUR(index_ellipse>0)))).^(1./CUR(index_ellipse>0))).*(PIZSL(index_ellipse>0)-PIch(index_ellipse>0))+PIch(index_ellipse>0);
    % In Choke Region
    %index_choke = grid.Wc==Wch;
    n_choke = options.points_in_SpL(3);
    for i = 1:length(options.Nc_vec)
        grid.PiC(end-n_choke:end,i) = linspace(PIch(end,i),options.lower_PiC,n_choke+1)';
    end
end
end