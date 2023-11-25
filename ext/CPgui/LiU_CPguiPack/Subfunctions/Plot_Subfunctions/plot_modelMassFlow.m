function [] = plot_modelMassFlow(map,Options,Param,axesHandle,Plot_Color)
%[] = plot_modelMassFlow(map,Options,Param,axesHandle)
%
% Plots the model pressure ratio vs mass flow in the given axes handle
%
%   map         	- Struct containing the map signals
%   Options         - Struct containing the plot options
%   Param       	- Compressor model parameters
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
        Plot_Color = '-b';
    case 2
        axesHandle = [];
        Plot_Color = '-b';
    case 3
        axesHandle = [];
        Plot_Color = '-b';
    case 4
        Plot_Color = '-b';
    case 5
  otherwise
end
%% Compute model signals
% Options
options_gen                     = Options;
options_gen.surgePercent        = 0;
options_gen.Nc_vec              = [0;map.Nc_vec];
options_gen.points_in_SpL       = [100,150,30];
options_gen.lower_PiC           = 0.3;
options_base                   	= Options;
options_base.Nc_vec             = options_gen.Nc_vec;
options_base.points_base_plot   = 150;
% Create Grids
gridF   = compute_grid(map,Param,options_gen);
gridF   = compute_PiC(map,Param,gridF,options_gen);
baseF   = compute_basefunctions(map,Param,options_base);

%% Do the plot
if isempty(axesHandle)
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name','Fitted and Optimized Compressor Map') 
    set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
    title('Fitted and Optimized Compressor Map');
    axesHandle = gca;
end

set(axesHandle,'defaulttextinterpreter','latex');
if verLessThan('matlab','8.4')
else
    set(axesHandle,'DefaultAxesColorOrder',colormap(parula(length(map.Nc_vec))));
end  
plot(axesHandle,map.WcCorr_M,map.PiC_M,'-x','LineWidth',2);
hold(axesHandle,'on');
% Plot Model
plot(axesHandle,gridF.Wc,gridF.PiC,Plot_Color);
plot(axesHandle,baseF.Wch_vec,baseF.PIch_vec,'ok');
plot(axesHandle,baseF.Wch_plot,baseF.PIch_plot,'-k');
plot(axesHandle,baseF.WZSL_vec,baseF.PIZSL_vec,'ok');
plot(axesHandle,baseF.WZSL_plot,baseF.PIZSL_plot,'-k');

xlabel(axesHandle,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(axesHandle,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);
grid(axesHandle,'on');
if options_gen.SurgeReverse
    axis(axesHandle,[-1.05.*Param.K0.*map.Wc_max_map,max(map.WcCorr_M(:))*1.05,0,max(map.PiC_M(:))*1.05])
else
    axis(axesHandle,[0,max(map.WcCorr_M(:))*1.05,0,max(map.PiC_M(:))*1.05])
end