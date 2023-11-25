function [] = plot_baseFunctions(map,Options,Param,ParamI,InitialGuess,Constraints,axesHandle)
%[] = plot_baseFunctions(map,Options,Param,ParamI,InitialGuess,Constraints)
%
% Plots each model base function vs corrected speed
%
%   map         	- Struct containing the map signals
%   Options         - Struct containing the plot options
%   Param       	- Struct containing the model parameters
%   ParamI          - Optional second set of parameters
%   InitialGuess  	- Struct with the used initial guesses, plotted for
%                     comparison
%   Constraints  	- Struct with the used constraints, plotted for
%                     reference of feasible space during the parameterization
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
    case 2
    case 3
        ParamI          = [];
        InitialGuess    = [];
        Constraints     = [];
        axesHandle      = [];
    case 4
        InitialGuess    = [];
        Constraints     = [];
        axesHandle      = [];
	case 5
        Constraints     = [];
        axesHandle      = [];
	case 6
        axesHandle      = [];
    otherwise
end
%% Compute model signals
% Options
options_base                    = Options;
options_base.Nc_vec             = [0;map.Nc_vec];
options_base.points_base_plot   = 150;
% Create Grids Initial
if ~isempty(ParamI)   
    baseI   = compute_basefunctions(map,ParamI,options_base);
%elseif isfield(ParamF,'CUR_vec')
    options_gen                     = Options;
    options_gen.surgePercent        = 0.6;
    options_gen.Nc_vec              = [0;map.Nc_vec];
    options_gen.points_in_SpL       = [50,150,30];
    options_gen.lower_PiC           = 0.3;
    gridF   = compute_grid(map,Param,options_gen);
    gridF   = compute_PiC(map,Param,gridF,options_gen);
end
if isfield(Param,'C_a')
    plot_row = 3;
    plot_col = 3;
else
    plot_row = 2;
    plot_col = 3;
end 
% Create Grids Final
base   = compute_basefunctions(map,Param,options_base);
if ~isempty(Constraints)
    PIZS_max    = options_base.F_PIZS([Constraints.UBcurr.PIZSL(1);Constraints.LBcurr.PIZSL(2)],base.Nc_plot./map.Nc_max_map,map);
    PIZS_min    = options_base.F_PIZS([Constraints.LBcurr.PIZSL(1);Constraints.UBcurr.PIZSL(2)],base.Nc_plot./map.Nc_max_map,map);
    WZS_max     = options_base.F_WZS([Constraints.UBcurr.WZSL(1);Constraints.LBcurr.WZSL(2)],base.Nc_plot./map.Nc_max_map,map);
    WZS_min     = options_base.F_WZS([Constraints.LBcurr.WZSL(1);Constraints.UBcurr.WZSL(2)],base.Nc_plot./map.Nc_max_map,map);
end
% base.WZSL_plot 	= options.F_WZS(param.WZSL,base.Nc_plot./map.Nc_max_map,map);
% base.PIZSL_plot = options.F_PIZS(param.PIZSL,base.Nc_plot./map.Nc_max_map,map);
%% Do the plot
if isempty(axesHandle)
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name', 'Base Functions') 
    set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
    axesHandle = gca;
end
set(axesHandle,'defaulttextinterpreter','latex');
hold(axesHandle,'on');
subplot(axesHandle);
subHandle = subplot(plot_row,plot_col,1);
p1 = plot(subHandle,base.Nc_plot,base.Wch_plot,'-b','LineWidth',1.5); hold on;
plot(subHandle,base.Nc_vec,base.Wch_vec,'ob','LineWidth',1.5);
if ~isempty(ParamI) 
    p3 = plot(subHandle,baseI.Nc_plot,baseI.Wch_plot,'--r','LineWidth',1.5);
	plot(subHandle,baseI.Nc_vec,baseI.Wch_vec,'or','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
elseif ~isempty(InitialGuess)
    p3 = plot(subHandle,InitialGuess.N_Ch,InitialGuess.W_chk ,'--*r','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Guess','location','best')
else
    %legend(p1(1),'Final Param','location','best')
end
axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.Wch_plot(:))*0.9,max(base.Wch_plot(:))*1.1])
ylabel(subHandle,'$\bar{W}_{ch} \quad [kg/s]$','interpreter', 'latex','FontSize',15);
xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
grid(subHandle,'on');

subHandle = subplot(plot_row,plot_col,2);
p1 = plot(subHandle,base.Nc_plot,base.WZSL_plot,'-b','LineWidth',1.5);hold on;
if ~isempty(Constraints) 
    p4 = plot(subHandle,base.Nc_plot,WZS_max,'--k','LineWidth',1.5);hold on;
    p5 = plot(subHandle,base.Nc_plot,WZS_min,':k','LineWidth',1.5);hold on;
end
plot(subHandle,base.Nc_vec,base.WZSL_vec,'ob','LineWidth',1.5)
if ~isempty(ParamI) 
    p3 = plot(subHandle,baseI.Nc_plot,baseI.WZSL_plot,'--r','LineWidth',1.5);
    plot(subHandle,baseI.Nc_vec,baseI.WZSL_vec,'or','LineWidth',1.5)
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
elseif ~isempty(InitialGuess)
    p3 = plot(subHandle,InitialGuess.N_ZSL,InitialGuess.W_ZSL ,'--*r','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Guess','location','best')
else
    %legend(p1(1),'Final Param','location','best')
end
axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.WZSL_plot(:))*0.9,max(base.WZSL_plot(:))*1.1])
ylabel(subHandle,'$\bar{W}_{ZSL} \quad [kg/s]$','interpreter', 'latex','FontSize',15);
xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
grid(subHandle,'on');

subHandle = subplot(plot_row,plot_col,3);
p1 = plot(subHandle,base.Nc_plot,base.PIch_plot,'-b','LineWidth',1.5);hold on;
plot(subHandle,base.Nc_vec,base.PIch_vec,'ob','LineWidth',1.5)
if ~isempty(ParamI) 
    p3 = plot(subHandle,baseI.Nc_plot,baseI.PIch_plot,'--r','LineWidth',1.5);
    plot(subHandle,baseI.Nc_vec,baseI.PIch_vec,'or','LineWidth',1.5)
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
elseif ~isempty(InitialGuess)
    p3 = plot(subHandle,InitialGuess.N_Ch,InitialGuess.PI_chk ,'--*r','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Guess','location','best')
else
    %legend(p1(1),'Final Param','location','best')
end
axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.PIch_plot(:))*0.9,max(base.PIch_plot(:))*1.1])
ylabel(subHandle,'$\Pi_{ch} \quad [-]$','interpreter', 'latex','FontSize',15);
xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
grid(subHandle,'on');

subHandle = subplot(plot_row,plot_col,4);
p1 = plot(subHandle,base.Nc_plot,base.PIZSL_plot,'-b','LineWidth',1.5);hold on;
plot(subHandle,base.Nc_vec,base.PIZSL_vec,'ob','LineWidth',1.5)
if ~isempty(Constraints) 
    p4 = plot(subHandle,base.Nc_plot,PIZS_max,'--k','LineWidth',1.5);hold on;
    p5 = plot(subHandle,base.Nc_plot,PIZS_min,':k','LineWidth',1.5);hold on;
end
if ~isempty(ParamI) 
    p3 = plot(subHandle,baseI.Nc_plot,baseI.PIZSL_plot,'--r','LineWidth',1.5);
    plot(subHandle,baseI.Nc_vec,baseI.PIZSL_vec,'or','LineWidth',1.5)
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
elseif ~isempty(InitialGuess)
    p3 = plot(subHandle,InitialGuess.N_ZSL,InitialGuess.PI_ZSL ,'--*r','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Guess','location','best')
else
    %legend(p1(1),'Final Param','location','best')
end
axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.PIZSL_plot(:))*0.9,max(base.PIZSL_plot(:))*1.15])
ylabel(subHandle,'$\Pi_{ZSL} \quad [-]$','interpreter', 'latex','FontSize',15);
xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
grid(subHandle,'on');

subHandle = subplot(plot_row,plot_col,5);
p1 = plot(subHandle,base.Nc_plot,base.CUR_plot,'-b','LineWidth',1.5);hold on;
plot(subHandle,base.Nc_vec,base.CUR_vec,'ob','LineWidth',1.5)
if ~isempty(ParamI) 
    p3 = plot(subHandle,baseI.Nc_plot,baseI.CUR_plot,'--r','LineWidth',1.5);
    plot(subHandle,baseI.Nc_vec,baseI.CUR_vec,'or','LineWidth',1.5)
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
elseif isfield(Param,'CUR_vec')
    p3 = plot(subHandle,map.Nc_vec,Param.CUR_vec ,'--*r','LineWidth',1.5);
    legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Guess','location','best')
else
    %legend(p1(1),'Final Param','location','best')
end
axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.CUR_plot(:))*0.9,max(base.CUR_plot(:))*1.15])
ylabel(subHandle,'$CUR \quad [-]$','interpreter', 'latex','FontSize',15);
xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
grid(subHandle,'on');

if  isfield(base,'A_plot')
    subHandle = subplot(plot_row,plot_col,6);
    p1 = plot(subHandle,base.Nc_plot,base.B_plot,'-b','LineWidth',1.5); hold on;
    plot(subHandle,base.Nc_vec,base.B_vec,'ob','LineWidth',1.5)
    if ~isempty(ParamI) 
        if isfield(baseI,'B_plot') 
            p3 = plot(subHandle,baseI.Nc_plot,baseI.B_plot,'--r','LineWidth',1.5);
            plot(subHandle,baseI.Nc_vec,baseI.B_vec,'or','LineWidth',1.5)
            legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
        end
    else
        %legend(p1(1),'Final Param','location','best')   
    end
    axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.B_plot(:))*0.9,max(base.B_plot(:))*1.15])
    ylabel(subHandle,'$B$','interpreter', 'latex','FontSize',15);
    xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
    grid(subHandle,'on');
    
    subHandle = subplot(plot_row,plot_col,7);
    p1 = plot(subHandle,base.Nc_plot,base.A_plot,'-b','LineWidth',1.5); hold on;
    plot(subHandle,base.Nc_vec,base.A_vec,'ob','LineWidth',1.5)
    if ~isempty(ParamI) 
        if isfield(baseI,'A_plot') 
            p3 = plot(subHandle,baseI.Nc_plot,baseI.A_plot,'--r','LineWidth',1.5);
            plot(subHandle,baseI.Nc_vec,baseI.A_vec,'or','LineWidth',1.5)
            legend(subHandle,[p1(1),p3(1)],'Final Param','Initial Param','location','best')
        end
    else
        %legend(p1(1),'Final Param','location','best')   
    end
    axis(subHandle,[0,max(base.Nc_plot(:))*1.15,min(base.A_plot(:))*0.9,max(base.A_plot(:))*1.15])
    ylabel(subHandle,'$A$','interpreter', 'latex','FontSize',15);
    xlabel(subHandle,'$\bar{N}_{c} \quad [rpm]$','interpreter', 'latex','FontSize',15);
    grid(subHandle,'on');
end