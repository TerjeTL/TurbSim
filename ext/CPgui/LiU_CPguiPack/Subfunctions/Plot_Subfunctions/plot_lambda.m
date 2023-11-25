function [] = plot_lambda(map,Options,Param,axesHandle,plotType)
%[] = plot_lambda(map,Options,Param,axesHandle,plotType)
%
% Plots the work coefficient vs flow coefficient or the enthalpy vs mass
% flow for the given map and parameter sets
%
%   map         	- Struct containing the map signals
%   Options         - Struct containing the plot options
%   Param       	- Struct with the parameters
%   axesHandle  	- axes handle where the plot has to be drawn
%   axesHandle  	- Type of plot to be done, either work coefficient vs
%                     flow coefficient or enthalpy vs mass flow
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
if isempty(Param)
else
    % Options
    extra_SpL                       = 0;
    Options.SurgeReverse            = 0;
    options_gen                     = Options;
    options_gen.lambda              = 1;
    options_gen.Ns                  = 0;
    options_gen.surgePercent        = 0.6;
    options_gen.Nc_vec              = [linspace(0,map.Nc_vec(1),extra_SpL+1)';map.Nc_vec(2:end)];
    options_gen.points_in_SpL       = [50,150,30];
    options_gen.lower_PiC           = 0.3;
    % Create Grid SpL
    grid_SpL        = compute_grid(map,Param,options_gen);
    grid_SpL        = compute_PiC(map,Param,grid_SpL,options_gen);
    grid_SpL        = compute_etaC(map,Param,grid_SpL,options_gen);
end
%% Do the plot

if isempty(axesHandle)
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name', ['Lambda plots for Compressor ',   map.MODEL]) 
    set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
    title(['Lambda plots for Compressor ',  map.MODEL]);
    axesHandle = gca;
end
set(axesHandle,'defaulttextinterpreter','latex');
if verLessThan('matlab','8.4')
else
    set(axesHandle,'DefaultAxesColorOrder',colormap(parula(length(map.Nc_vec))));
end 
if plotType==1
    plot(axesHandle,map.PhiC_M,map.lambda_M,'-x','LineWidth',2);
    h = legend(axesHandle,num2str(map.MachC_vec),'location','best');
    v = get(h,'title');
    set(v,'string','Mach Number');
    hold(axesHandle,'on');
    if ~isempty(Param)
        plot(axesHandle,grid_SpL.Phi(:,extra_SpL+1:end),grid_SpL.lambda(:,extra_SpL+1:end),'--','LineWidth',1.5);
        plot(axesHandle,grid_SpL.Phi(:,1:extra_SpL),grid_SpL.lambda(:,1:extra_SpL),'--r','LineWidth',1.5);
    end
    ylabel(axesHandle,'$\lambda \quad [-]$','interpreter', 'latex','FontSize',13);
    xlabel(axesHandle,'$\Phi_{01} \quad [-]$','interpreter', 'latex','FontSize',13);
    grid(axesHandle,'on');
    axis(axesHandle,[min(map.PhiC_M(:))*0.7,max(map.PhiC_M(:))*1.3,0,max(map.lambda_M(:).*1.05)])

else
    plot(axesHandle,map.WcCorr_M,map.hoC_M,'-x','LineWidth',2);
  	h = legend(axesHandle,num2str(map.MachC_vec),'location','best');
    v = get(h,'title');
    set(v,'string','Mach Number');
    hold(axesHandle,'on');
    if ~isempty(Param)
        plot(axesHandle,grid_SpL.Wc(:,extra_SpL+1:end),grid_SpL.Dh(:,extra_SpL+1:end),'--','LineWidth',1.5);
        plot(axesHandle,grid_SpL.Wc(:,1:extra_SpL),grid_SpL.Dh(:,1:extra_SpL),'--r','LineWidth',1.5);
    end
	ylabel(axesHandle,'$\Delta h_0 \quad [kJ/kg]$','interpreter', 'latex','FontSize',13);
    xlabel(axesHandle,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
    grid(axesHandle,'on');
    axis(axesHandle,[0,max(map.WcCorr_M(:))*1.05,0,max(map.hoC_M(:))*1.05])
end