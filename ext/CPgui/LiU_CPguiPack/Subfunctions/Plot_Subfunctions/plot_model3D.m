function [] = plot_model3D(map,Options,Param,axesHandle)
%[] = plot_model3D(map,Options,Param,axesHandle)
%
% Plots a number of speed lines in the 3D space formed by mass flow,
% pressure ratio and efficiency
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
        axesHandle = [];
    case 2
        axesHandle = [];
    case 3
        axesHandle = [];
    case 4
    otherwise
end

%% Compute model signals
% Options
Options.SurgeReverse            = 0;
options_gen                     = Options;
options_gen.surgePercent        = 0.6;
options_gen.Nc_vec              = map.Nc_vec;
options_gen.points_in_SpL       = [50,150,30];
options_gen.lower_PiC           = 0.3;
options_base                    = Options;
options_base.Nc_vec             = options_gen.Nc_vec;
options_base.points_base_plot   = 150;
% Create Grid SpL
grid_SpL   = compute_grid(map,Param,options_gen);
grid_SpL   = compute_PiC(map,Param,grid_SpL,options_gen);
grid_SpL   = compute_etaC(map,Param,grid_SpL,options_gen);
% Create thin grid for the contour
options_extense                	= Options;
options_extense.surgePercent    = 0.6;
options_extense.Nc_vec          = linspace(0,map.Nc_vec(end),100)';
options_extense.points_in_SpL   = [70,200,50];
options_extense.lower_PiC       = 0.25;
eta_cMAX                        = max(map.etaC_M(:))*1.02;
options_extense.zlevs           = [0:0.2:0.4,0.3:0.05:eta_cMAX-0.1,eta_cMAX-0.1:0.01:eta_cMAX-0.05,eta_cMAX-0.05:0.005:eta_cMAX-0.02,eta_cMAX-0.02:0.001:eta_cMAX];

grid_extense = compute_grid(map,Param,options_extense);
grid_extense = compute_PiC(map,Param,grid_extense,options_extense);
grid_extense = compute_etaC(map,Param,grid_extense,options_extense);
%% Do the plot
if isempty(axesHandle)
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name','Compressor model 3D plot vs Measurements') 
    set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
    title('Compressor model 3D plot vs Measurements');
    axesHandle = gca;
end
set(axesHandle,'defaulttextinterpreter','latex');
plot3(axesHandle,map.WcCorr_M,map.PiC_M,map.etaC_M,'--x','LineWidth',3);
hold(axesHandle,'on');
plot3(axesHandle,grid_SpL.Wc,grid_SpL.PiC,grid_SpL.etaC,'-','LineWidth',3);
plot3(axesHandle,grid_extense.Wc,grid_extense.PiC,grid_extense.etaC,'-k');
xlabel(axesHandle,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(axesHandle,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);
zlabel(axesHandle,'$\eta_c \quad [-]$','interpreter', 'latex','FontSize',13);
grid(axesHandle,'on');
axis(axesHandle,[0,max(map.WcCorr_M(:))*1.05,0,max(map.PiC_M(:))*1.05,0,max(map.etaC_M(:))*1.05])
