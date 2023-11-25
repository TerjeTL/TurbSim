function [] = plot_modelContour(map,Options,Param,axesHandle)
%[] = plot_modelContour(map,Options,Param,axesHandle)
%
% Plots the model pressure ratio vs mass flow in the given axes handle
% with the efficiency extrapolation as contour levels.
%
%   map         	- Struct containing the map signals
%   Options         - Struct containing the plot options
%   Param       	- Struct with the parameters
%   axesHandle  	- axes handle where the plot has to be drawn
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
switch nargin
    case 0
    case 1
        Subplot = 0;
        Lwidth  = 2;
        axesHandle = [];
    case 2
        Subplot = 0;
        Lwidth  = 2;
        axesHandle = [];
    case 3
        Subplot = 0;
        Lwidth  = 2;
        axesHandle = [];
    case 4
        Subplot = 1;
        Lwidth  = 1;
    otherwise
end

%% Compute model signals
% Options
if ~isfield(Options,'SurgeReverse')
    Options.SurgeReverse = 0;
end
options_gen                     = Options;
options_gen.surgePercent        = 0;
options_gen.Nc_vec              = [0;map.Nc_vec];
options_gen.points_in_SpL       = [100,150,30];
options_gen.lower_PiC           = 0.3;
options_base                  	= Options;
options_base.Nc_vec             = options_gen.Nc_vec;
options_base.points_base_plot   = 150;
% Create Grid SpL
grid_SpL   = compute_grid(map,Param,options_gen);
grid_SpL   = compute_PiC(map,Param,grid_SpL,options_gen);
base_SpL   = compute_basefunctions(map,Param,options_base);
% Create thin grid for the contour
options_extense                	= Options;
options_extense.SurgeReverse    = 0;
options_extense.surgePercent    = 0.2;
options_extense.Nc_vec          = linspace(0,map.Nc_vec(end),100)';
options_extense.points_in_SpL   = [70,200,50];
options_extense.lower_PiC       = 0.25;

grid_extense = compute_grid(map,Param,options_extense);
grid_extense = compute_PiC(map,Param,grid_extense,options_extense);
grid_extense = compute_etaC(map,Param,grid_extense,options_extense);

eta_cMAX                        = max(map.etaC_M(:))*1.02;
if ~Subplot
    options_extense.zlevs           = [0:0.2:0.4,0.3:0.05:eta_cMAX-0.1,eta_cMAX-0.1:0.01:eta_cMAX-0.05,eta_cMAX-0.05:0.005:eta_cMAX-0.02,eta_cMAX-0.02:0.002:eta_cMAX];
else
    options_extense.zlevs           = [0:0.2:0.5,0.5:0.1:eta_cMAX-0.1,eta_cMAX-0.1:0.05:eta_cMAX-0.05,eta_cMAX-0.05:0.01:eta_cMAX];
end
%% Do the plot
if isempty(axesHandle)
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name','Compressor Efficiency Contour Map') 
    set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
    title('Compressor Efficiency Contour Map');
    axesHandle = gca;
end
set(axesHandle,'defaulttextinterpreter','latex');
if verLessThan('matlab','8.4')
else
    set(axesHandle,'DefaultAxesColorOrder',colormap(parula(length(map.Nc_vec))));
end 
plot(axesHandle,map.WcCorr_M,map.PiC_M,'-x','LineWidth',Lwidth+0.5);
hold(axesHandle,'on');
if verLessThan('matlab','8.4')
else
    set(axesHandle,'DefaultAxesColorOrder',colormap(jet(length(options_extense.zlevs))));
end 
contour(axesHandle,grid_extense.Wc,grid_extense.PiC,grid_extense.etaC,options_extense.zlevs,'ShowText','on');
plot(axesHandle,grid_SpL.Wc,grid_SpL.PiC,'-b','LineWidth',Lwidth);
plot(axesHandle,base_SpL.Wch_vec,base_SpL.PIch_vec,'ok','LineWidth',Lwidth);
plot(axesHandle,base_SpL.Wch_plot,base_SpL.PIch_plot,'-k','LineWidth',Lwidth);
plot(axesHandle,base_SpL.WZSL_vec,base_SpL.PIZSL_vec,'ok','LineWidth',Lwidth);
plot(axesHandle,base_SpL.WZSL_plot,base_SpL.PIZSL_plot,'-k','LineWidth',Lwidth);
xlabel(axesHandle,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(axesHandle,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);
grid(axesHandle,'on');
if Options.SurgeReverse
    axis(axesHandle,[-1.05.*Param.K0.*map.Wc_max_map,max(grid_SpL.Wc(:))*1.05,0,max(map.PiC_M(:))*1.05])
else
    axis(axesHandle,[0,max(grid_SpL.Wc(:))*1.05,0,max(map.PiC_M(:))*1.05])
end