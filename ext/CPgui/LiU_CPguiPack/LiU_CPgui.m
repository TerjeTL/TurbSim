function varargout = LiU_CPgui(varargin)
%       Version 1.2: 2017-12-20
%       Copyright (C) 2017, Xavier Llamas
%
%       This file is part of LiU_CPgui.
%
%       LiU_CPgui is free software: you can redistribute it and/or modify
%       it under the terms of the GNU Lesser General Public License as
%       published by the Free Software Foundation, version 3 of the License.
%
%       This package is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%       GNU Lesser General Public License for more details.
%
%       You should have received a copy of the GNU Lesser General Public License
%       along with LiU_CPgui.  If not, see <http://www.gnu.org/licenses/>.

% LIU_CPGUI MATLAB code for LiU_CPgui.fig
%      LIU_CPGUI, by itself, creates a new LIU_CPGUI or raises the existing
%      singleton*.
%
%      H = LIU_CPGUI returns the handle to a new LIU_CPGUI or the handle to
%      the existing singleton*.
%
%      LIU_CPGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIU_CPGUI.M with the given input arguments.
%
%      LIU_CPGUI('Property','Value',...) creates a new LIU_CPGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LiU_CPgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LiU_CPgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LiU_CPgui

% Last Modified by GUIDE v2.5 19-Dec-2017 09:36:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LiU_CPgui_OpeningFcn, ...
                   'gui_OutputFcn',  @LiU_CPgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LiU_CPgui is made visible.
function LiU_CPgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LiU_CPgui (see VARARGIN)

% Retrieve input data
handles.map         = varargin{1};

% Default options //get them from GUI
handles.BaseFunc.Wch            = 'F_WChL_atan'; 
handles.BaseFunc.PIch           = 'F_PIChL';     
handles.BaseFunc.WZS            = 'F_WZSL';     
handles.BaseFunc.PIZS           = 'F_PIZSL';   
handles.BaseFunc.CUR            = 'F_CUR';      
handles.BaseFunc.PI0            = 'F_PI0';        
handles.BaseFunc.A              = 'F_a_SAE';   
handles.BaseFunc.B              = 'F_b_SAE';      
handles.BaseFunc.K_loss         = 'F_kloss'; 
% Initialize, initialization boolean
handles.Initialized        = 0;
% Plot Measured map
handles.f = hObject;

hold(handles.ax,'on')
xlabel(handles.ax,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(handles.ax,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);

hold(handles.ax2,'on')
xlabel(handles.ax2,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(handles.ax2,'$\eta_c \quad [-]$','interpreter', 'latex','FontSize',13);
% Check that user has the Optimization Toolbox with lsqnonlin installed.
has_lsqnonlin = ~isempty(which('lsqnonlin'));
if ~has_lsqnonlin
	% User does not have the toolbox installed.
   	warningstring = 'The default solver lsqnonlin, from the Optimization Toolbox, is not available. LSOptim will be used instead.';
  	h = errordlg(warningstring,'lsqnonlin not available','modal');
    waitfor(h);
    set(handles.popupmenu6,'Value',2);
    handles.LSOptim          = 1;   
  	% Theshold Constraints tighter if LSOptim is used)
    set(handles.edit8,'String','0.3')
    set(handles.edit9,'String','0.3')
    set(handles.edit10,'String','0.3')
    set(handles.edit11,'String','0.3')
    set(handles.edit12,'String','0.3')
    set(handles.edit13,'String','0.3')
    set(handles.edit14,'String','0.3')
    set(handles.edit15,'String','0.3')
    handles.Constraints.C_Wch.Th_lb    	= str2double(get(handles.edit8,'String'));
    handles.Constraints.C_PIch.Th_lb   	= str2double(get(handles.edit9,'String'));
    handles.Constraints.C_Wzs.Th_lb   	= str2double(get(handles.edit10,'String'));
    handles.Constraints.C_PIzs.Th_lb   	= str2double(get(handles.edit11,'String'));
    handles.Constraints.C_Wch.Th_ub   	= str2double(get(handles.edit12,'String'));
    handles.Constraints.C_PIch.Th_ub   	= str2double(get(handles.edit13,'String'));
    handles.Constraints.C_Wzs.Th_ub   	= str2double(get(handles.edit14,'String'));
    handles.Constraints.C_PIzs.Th_ub   	= str2double(get(handles.edit15,'String'));
else
    set(handles.popupmenu6,'Value',1);
    handles.LSOptim          = 0;    
end
% Define output
handles.output = hObject;
% Init Base functions
% TEXT annotations need an axes as parent so create an invisible axes which
% is as big as the button panel
handles.laxis = axes('parent',handles.InitBaseFunctions,'units','normalized','position',[0 0 1 1],'visible','off');
% Find all static text UICONTROLS whose 'Tag' starts with latex_
lbls = findobj(hObject,'-regexp','tag','latex_Base*');
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex','fontsize',13);
end
% Heat correction options
% TEXT annotations need an axes as parent so create an invisible axes which
% is as big as the button panel
handles.laxis = axes('parent',handles.InitHeatParameters,'units','normalized','position',[0 0 1 1],'visible','off');
% Find all static text UICONTROLS whose 'Tag' starts with latex_
lbls = findobj(hObject,'-regexp','tag','latex_Heat*');
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex','fontsize',13);
end
% Constraint options
% TEXT annotations need an axes as parent so create an invisible axes which
% is as big as the button panel
handles.laxis = axes('parent',handles.ConstraintOptions,'units','normalized','position',[0 0 1 1],'visible','off');
% Find all static text UICONTROLS whose 'Tag' starts with latex_
lbls = findobj(hObject,'-regexp','tag','latex_Const*');
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex','fontsize',13);
end
[handles.map,requiredMapData,MapFormat] = CheckMapData(handles.map);
handles.AlgorithmOK = 1;
% Check map data
if ~requiredMapData
   	warningstring = 'Compressor Map has missing signals';
  	h = errordlg(warningstring,'Data error','modal');
    waitfor(h);
    handles.AlgorithmOK = 0;
elseif ~MapFormat
  	warningstring = 'Compressor Map signals are not column vectors';
  	h = errordlg(warningstring,'Data error','modal');
    waitfor(h);
    handles.AlgorithmOK = 0;
else
    set(handles.pushbutton16, 'Enable', 'on') 
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LiU_CPgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = LiU_CPgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Parameterizes the ellipse model and returns the parameters
disp ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp ('Ellipse Mass Flow Model Parameterization')
disp ('----------------------------------------')
% select initial guess based on User Input
InitGuess.Param_i                       = handles.InitEllipseParam.vector;
InitGuess.Delta_i                       = handles.InitEllipseParam.delta_Wc;
[handles.paramEll,handles.errorsEll] 	= Ellipse_Parameterization(handles.map,InitGuess,handles.Constraints,handles.Options);
% Display Current Relative Errors
disp ('---------------------------------------------------------')
disp ('Absolute Relative Errors of the Ellipse Mass Flow Model')
disp ('---------------------------------------------------------')
% Relative error on PiC dimension
disp ('Absolute Relative error of Pi_c:')
disp ('Mean value [%]:')
disp (handles.errorsEll.PiC_mean_rel)
disp ('Max. value [%]:')
disp (max(handles.errorsEll.PiC_rel))
% Relative error on Wc_corr dimension
disp ('Absolute Relative error of W_c:')
disp ('Mean value [%]:')
disp (handles.errorsEll.Wc_mean_rel)
disp ('Max. value [%]:')
disp (max(handles.errorsEll.Wc_rel))
% Update handles
handles.popupmenu_Plotparameters{2}         = 'paramEll';
handles.popupmenu_EllipseInitparameters{2}  = 'paramEll';
handles.popupmenu_EffInitparameters{2}      = 'paramEll';

set(handles.popupmenu7,'Value',2);
handles.Curr_param                          = handles.paramEll;
set(handles.popupmenu8,'String',handles.PlotTypes_sh);
handles.Plot_color 	= '-r';
set(handles.popupmenu7,'String',handles.popupmenu_Plotparameters);
set(handles.popupmenu10,'String',handles.popupmenu_EllipseInitparameters);
set(handles.popupmenu12,'String',handles.popupmenu_EffInitparameters);
handles.InitEffParam = handles.paramEll;
% Update handles structure
contents = cellstr(get(handles.popupmenu7,'String'));
if strcmp(contents{get(handles.popupmenu7,'Value')},'paramEll')
    handles.Curr_param = handles.paramEll;
end
% Enable GUI buttons
set(handles.popupmenu12, 'Enable', 'on')
set(handles.popupmenu12,'Value',2);
set(handles.pushbutton9, 'Enable', 'on') 
guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Parameterizes the ellipse model and returns the parameters
% Check if the edit tick was enabled, disable and update initial guess
% before proceding
if handles.DragInitGuess
% make sure no toolbar buttons are enabled
children = get( get(handles.uitoggletool1,'Parent'),'Children');
set(children,'State','off');
zoom off
pan off
rotate3d off
datacursormode off
% Change color of the plotted lines
set(handles.lH,'Color','k');
% Disable drag callbacks and update initial guess
    if isempty(getappdata(handles.ax,'InGuW_chk'))
    else   
        handles.InitGuess.W_chk    = getappdata(handles.ax,'InGuW_chk');
        handles.InitGuess.W_ZSL    = getappdata(handles.ax,'InGuW_ZSL');
        handles.InitGuess.PI_chk   = getappdata(handles.ax,'InGuPI_chk');
        handles.InitGuess.PI_ZSL   = getappdata(handles.ax,'InGuPI_ZSL');
    end
    set(handles.f,'WindowButtonDownFcn',''); % Defining what happens when clicking
    set(handles.f,'WindowButtonUpFcn','');
    handles.DragInitGuess = 0;
    set(handles.checkbox4,'Value',0);
end
% Initialize Ellipse model
disp ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp ('Initialization of the Ellipse Model')
disp ('-----------------------------------')
handles.paramEll_i              = Ellipse_Initialization(handles.map, handles.InitGuess, handles.Constraints, handles.Options);
handles.popupmenu_Plotparameters{1}         = 'paramEll_i';
handles.popupmenu_EllipseInitparameters{1}  = 'paramEll_i';
handles.popupmenu_EffInitparameters{1}      = 'paramEll_i';
set(handles.popupmenu7,'Value',1);
set(handles.popupmenu8,'Value',1);
handles.Curr_param                          = handles.paramEll_i;
set(handles.popupmenu8,'String',handles.PlotTypes_sh);
handles.InitEllipseParam                    = handles.paramEll_i;
handles.InitEffParam                        = handles.paramEll_i;
handles.Plot_color                          = '-g';

set(handles.popupmenu7,'String',handles.popupmenu_Plotparameters);
set(handles.popupmenu10,'String',handles.popupmenu_EllipseInitparameters);
set(handles.popupmenu12,'String',handles.popupmenu_EffInitparameters);
% Enable GUI buttons if Initialization went ok
set(handles.popupmenu10, 'Enable', 'on') 
set(handles.pushbutton3, 'Enable', 'on') 
set(handles.popupmenu7, 'Enable', 'on') 
set(handles.popupmenu8, 'Enable', 'on') 
set(handles.pushbutton28, 'Enable', 'on') 
set(handles.pushbutton29, 'Enable', 'on') 
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
handles.ax = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
handles.ax2 = hObject;
guidata(hObject, handles);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles.holdPlot = get(hObject,'Value');
%ax = get(handles.f,'CurrentAxes');
if handles.holdPlot
    hold(handles.ax,'on')
else
    if ~handles.DragInitGuess
    hold(handles.ax,'off')
    end
%     if verLessThan('matlab','8.4')
%         set(handles.ax,'DefaultAxesColorOrder',colormap(jet(length(map.Nc_vec))));
%     else
%         set(handles.ax,'DefaultAxesColorOrder',colormap(parula(length(map.Nc_vec))));
%     end  
end
guidata(hObject, handles);

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
handles.holdPlot2 = get(hObject,'Value');
%ax = get(handles.f,'CurrentAxes');
if handles.holdPlot2
    hold(handles.ax2,'on')
else
    hold(handles.ax2,'off')
%    	if verLessThan('matlab','8.4')
%         set(handles.ax,'DefaultAxesColorOrder',colormap(jet(length(map.Nc_vec))));
%     else
%         set(handles.ax,'DefaultAxesColorOrder',colormap(parula(length(map.Nc_vec))));
%     end  
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Find an initial guess for the efficiency parameters
disp ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp ('Efficiency Model Initialization')
disp ('-------------------------------')
handles.paramT_i  = Efficiency_Initialization(handles.map, handles.InitEffParam, handles.InitGuess, handles.Constraints, handles.Options);
handles.popupmenu_Plotparameters{3}         = 'paramT_i';
handles.popupmenu_TotalInitparameters{1}	= 'paramT_i';
set(handles.popupmenu7,'String',handles.popupmenu_Plotparameters);
set(handles.popupmenu11,'String',handles.popupmenu_TotalInitparameters);
% default initial guess for the total parameterization
handles.InitTotalParam       	= handles.paramT_i;
handles.Curr_param           	= handles.paramT_i;
set(handles.popupmenu7,'Value',3);
set(handles.popupmenu8,'String',handles.PlotTypes);
handles.Plot_color  = '-k';
set(handles.popupmenu8,'Value',2);
% Enable GUI buttons
set(handles.pushbutton12, 'Enable', 'on') 
set(handles.popupmenu11, 'Enable', 'on') 
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Init_eff_vec    = [handles.InitEffParam.vector;handles.InitGuess.param_b;handles.InitGuess.param_a;handles.InitGuess.C_loss];
%Init_eff_par	= compute_parameters(Init_eff_par,handles.InitEffParam.delta_Wc,handles.Options);
In              = handles.Options.ParvecI;
if handles.Options.n_Klossparam == 0
    Init_eff_par.vector   	= Init_eff_vec;
    Init_eff_par.delta_Wc 	= handles.InitEffParam.delta_Wc;
    Init_eff_par.C_Wch    	= Init_eff_vec(1:In(1));
    Init_eff_par.C_PIch    	= Init_eff_vec(In(1)+1:In(2));
    Init_eff_par.C_Wzs     	= Init_eff_vec(In(2)+1:In(3));
    Init_eff_par.C_PIzs   	= Init_eff_vec(In(3)+1:In(4));
    Init_eff_par.C_cur    	= Init_eff_vec(In(4)+1:In(5));
    Init_eff_par.Gamma_PIcs	= 0.5;
    Init_eff_par.C_s        = Init_eff_vec(In(5)+1:In(6));
    Init_eff_par.C_b      	= Init_eff_vec(In(6)+1:In(7));
    Init_eff_par.C_a      	= Init_eff_vec(In(7)+1:In(8));
    Init_eff_par.C_loss     = 0;
else
    Init_eff_par.vector   	= Init_eff_vec;
    Init_eff_par.delta_Wc 	= handles.InitEffParam.delta_Wc;
    Init_eff_par.C_Wch    	= Init_eff_vec(1:In(1));
    Init_eff_par.C_PIch    	= Init_eff_vec(In(1)+1:In(2));
    Init_eff_par.C_Wzs     	= Init_eff_vec(In(2)+1:In(3));
    Init_eff_par.C_PIzs   	= Init_eff_vec(In(3)+1:In(4));
    Init_eff_par.C_cur      	= Init_eff_vec(In(4)+1:In(5));
    Init_eff_par.Gamma_PIcs	= 0.5;
    Init_eff_par.C_s        = Init_eff_vec(In(5)+1:In(6));
    Init_eff_par.C_b      	= Init_eff_vec(In(6)+1:In(7));
    Init_eff_par.C_a      	= Init_eff_vec(In(7)+1:In(8));
    Init_eff_par.C_loss     = Init_eff_vec(In(8)+1:In(9));
end
Init_eff_par.F_Wch          = handles.Options.F_Wch;
Init_eff_par.F_PIch      	= handles.Options.F_PIch;
Init_eff_par.F_WZS          = handles.Options.F_WZS;
Init_eff_par.F_PIZS     	= handles.Options.F_PIZS;
Init_eff_par.F_CUR          = handles.Options.F_CUR;
Init_eff_par.F_PI0          = handles.Options.F_PI0;
Init_eff_par.F_a            = handles.Options.F_a;
Init_eff_par.F_b        	= handles.Options.F_b;
Init_eff_par.F_kloss        = handles.Options.F_kloss;
plot_modelEfficiency(handles.map,handles.Options,Init_eff_par,handles.ax2,'-r')
if handles.holdPlot2
    hold(handles.ax2,'on')
else
    hold(handles.ax2,'off')
end
% Enable GUI buttons
set(handles.pushbutton10, 'Enable', 'on') 
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp ('Total Mass flow and Efficiency Parameterization')
disp ('------------------------------------------------')
InitGuess.Param_i                   = handles.InitTotalParam.vector;
InitGuess.Delta_i                   = handles.InitTotalParam.delta_Wc;
[handles.paramT,handles.errorsT]  	= Complete_Parameterization(handles.map, InitGuess, handles.Constraints, handles.Options);
% Display Current Relative Errors
disp ('---------------------------------------------------------')
disp ('Absolute Relative Errors of the Complete Compressor Model')
disp ('---------------------------------------------------------')
% Relative error on PiC dimension
disp ('Absolute Relative error of Pi_c:')
disp ('Mean value [%]:')
disp (handles.errorsT.PiC_mean_rel)
disp ('Max. value [%]:')
disp (max(handles.errorsT.PiC_rel))
% Relative error on Wc_corr dimension
disp ('Absolute Relative error of W_c:')
disp ('Mean value [%]:')
disp (handles.errorsT.Wc_mean_rel)
disp ('Max. value [%]:')
disp (max(handles.errorsT.Wc_rel))
% Relative error on Wc_corr dimension
disp ('Absolute Relative error of eta_c:')
disp ('Mean value [%]:')
disp (handles.errorsT.eta_mean_rel)
disp ('Max. value [%]:')
disp (max(handles.errorsT.eta_rel))
% Update handles
handles.popupmenu_Plotparameters{4} = 'paramT';
set(handles.popupmenu7,'Value',4);
handles.Curr_param                  = handles.paramT;
set(handles.popupmenu8,'String',handles.PlotTypes);
handles.Plot_color 	= '-b';
handles.popupmenu_TotalInitparameters{2}   	= 'paramT';
set(handles.popupmenu7,'String',handles.popupmenu_Plotparameters);
set(handles.popupmenu11,'String',handles.popupmenu_TotalInitparameters);
% Update handles structure
contents = cellstr(get(handles.popupmenu7,'String'));
if strcmp(contents{get(handles.popupmenu7,'Value')},'paramT')
    handles.Curr_param = handles.paramT;
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
init  = initialize_algorithm(handles.map,handles.BaseFunc,handles.k_q,handles.Constraints,handles.LSOptim);
handles.Initialized        = 1;
handles.map                = init.map;
handles.Options            = init.Options;
handles.InitGuess          = init.InitGuess;
handles.Constraints        = init.Constraints;
if verLessThan('matlab','8.4')
else
    set(groot,'DefaultAxesColorOrder',colormap(parula(length(handles.map.Nc_vec))));
end  
% Enable buttons after the execution and no errors found
set(handles.popupmenu13, 'Enable', 'on')
set(handles.pushbutton30, 'Enable', 'on') 
set(handles.pushbutton20, 'Enable', 'on')
% Disable buttons in case the algorithm is being reinitialized
set(handles.popupmenu7, 'Enable', 'off')
set(handles.popupmenu8, 'Enable', 'off')
set(handles.pushbutton28, 'Enable', 'off') 
set(handles.pushbutton3, 'Enable', 'off')
set(handles.popupmenu10, 'Enable', 'off')
set(handles.pushbutton29, 'Enable', 'off')
set(handles.pushbutton4, 'Enable', 'off')
set(handles.checkbox4, 'Enable', 'off')
set(handles.pushbutton9, 'Enable', 'off')
set(handles.pushbutton10, 'Enable', 'off')
set(handles.popupmenu12, 'Enable', 'off')
set(handles.popupmenu11, 'Enable', 'off')
set(handles.pushbutton12, 'Enable', 'off')
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox4
hold(handles.ax,'on')
handles.DragInitGuess   = get(handles.checkbox4,'Value');
% make sure no toolbar buttons are enabled
children = get( get(handles.uitoggletool1,'Parent'),'Children');
set(children,'State','off');
zoom off
pan off
rotate3d off
datacursormode off
% Enable or disable click button callbacks
if handles.DragInitGuess
    % Change color of the plotted lines
    set(handles.lH,'Color','r');
    % Enable drag callbacks
    set(handles.lH,'hittest','off'); % so you can click on the Markers
    set(handles.f,'WindowButtonDownFcn',{@getCoord,handles}); % Defining what happens when clicking
    set(handles.f,'WindowButtonUpFcn',{@dragCoord,handles}); % Defining what happens when releasing the button
else
	% Change color of the plotted lines
    set(handles.lH,'Color','k');
    % Disable drag callbacks and update initial guess
    if isempty(getappdata(handles.ax,'InGuW_chk'))
    else   
        handles.InitGuess.W_chk    = getappdata(handles.ax,'InGuW_chk');
        handles.InitGuess.W_ZSL    = getappdata(handles.ax,'InGuW_ZSL');
        handles.InitGuess.PI_chk   = getappdata(handles.ax,'InGuPI_chk');
        handles.InitGuess.PI_ZSL   = getappdata(handles.ax,'InGuPI_ZSL');
    end
    set(handles.f,'WindowButtonDownFcn',''); % Defining what happens when clicking
    set(handles.f,'WindowButtonUpFcn','');
    if handles.holdPlot
    hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SpL_grid            = getEllipses(handles.InitGuess);
handles.Ell         = plot(handles.ax,SpL_grid.Wc,SpL_grid.PiC,'-b');hold(handles.ax,'on');
handles.lH(1)       = plot(handles.ax,handles.InitGuess.W_chk,handles.InitGuess.PI_chk,'-ok','LineWidth',2);
handles.lH(2)       = plot(handles.ax,handles.InitGuess.W_ZSL,handles.InitGuess.PI_ZSL,'-ok','LineWidth',2);
grid(handles.ax,'on');
axis(handles.ax,[0,max(handles.map.WcCorr_M(:))*1.05,0,max(handles.map.PiC_M(:))*1.05])
xlabel(handles.ax,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
ylabel(handles.ax,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);
if handles.holdPlot
    hold(handles.ax,'on')
else
    hold(handles.ax,'off')
end
% Enable GUI commands
set(handles.checkbox4, 'Enable', 'on') 
set(handles.pushbutton4, 'Enable', 'on') 
% Update handles structure
guidata(hObject, handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.k_q = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit2 = hObject;
set(handles.edit2,'String','0')
handles.k_q   	= str2double(get(handles.edit2,'String'));
% Update handles structure
guidata(hObject, handles)

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'Atan')
    handles.BaseFunc.Wch            = 'F_WChL_atan'; 
else
    handles.BaseFunc.Wch            = 'F_WChL_switch'; 
end
if handles.Initialized 
    % Disable buttons because the algorithm needs to be initialized again
    set(handles.popupmenu7, 'Enable', 'off')
    set(handles.popupmenu8, 'Enable', 'off')
    set(handles.pushbutton28, 'Enable', 'off') 
    set(handles.pushbutton3, 'Enable', 'off')
    set(handles.popupmenu10, 'Enable', 'off')
    set(handles.pushbutton29, 'Enable', 'off')
    set(handles.pushbutton4, 'Enable', 'off')
    set(handles.checkbox4, 'Enable', 'off')
    set(handles.pushbutton9, 'Enable', 'off')
    set(handles.pushbutton10, 'Enable', 'off')
    set(handles.popupmenu12, 'Enable', 'off')
    set(handles.popupmenu11, 'Enable', 'off')
    set(handles.pushbutton12, 'Enable', 'off')
    set(handles.pushbutton20, 'Enable', 'off') 
    set(handles.popupmenu13, 'Enable', 'off') 
    set(handles.pushbutton30, 'Enable', 'off') 
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.InitPanel = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uibuttongroup4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.InitBaseFunctions = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uibuttongroup3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.InitHeatParameters = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'SAE')
    handles.BaseFunc.A            = 'F_a_SAE'; 
else
    handles.BaseFunc.A            = 'F_a_ASME'; 
end
if handles.Initialized 
    % Disable buttons because the algorithm needs to be initialized again
    set(handles.popupmenu7, 'Enable', 'off')
    set(handles.popupmenu8, 'Enable', 'off')
    set(handles.pushbutton28, 'Enable', 'off') 
    set(handles.pushbutton3, 'Enable', 'off')
    set(handles.popupmenu10, 'Enable', 'off')
    set(handles.pushbutton29, 'Enable', 'off')
    set(handles.pushbutton4, 'Enable', 'off')
    set(handles.checkbox4, 'Enable', 'off')
    set(handles.pushbutton9, 'Enable', 'off')
    set(handles.pushbutton10, 'Enable', 'off')
    set(handles.popupmenu12, 'Enable', 'off')
    set(handles.popupmenu11, 'Enable', 'off')
    set(handles.pushbutton12, 'Enable', 'off')
    set(handles.pushbutton20, 'Enable', 'off') 
    set(handles.popupmenu13, 'Enable', 'off') 
    set(handles.pushbutton30, 'Enable', 'off') 
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'SAE')
    handles.BaseFunc.B            = 'F_b_SAE'; 
else
    handles.BaseFunc.B            = 'F_b_ASME'; 
end
if handles.Initialized 
    % Disable buttons because the algorithm needs to be initialized again
    set(handles.popupmenu7, 'Enable', 'off')
    set(handles.popupmenu8, 'Enable', 'off')
    set(handles.pushbutton28, 'Enable', 'off') 
    set(handles.pushbutton3, 'Enable', 'off')
    set(handles.popupmenu10, 'Enable', 'off')
    set(handles.pushbutton29, 'Enable', 'off')
    set(handles.pushbutton4, 'Enable', 'off')
    set(handles.checkbox4, 'Enable', 'off')
    set(handles.pushbutton9, 'Enable', 'off')
    set(handles.pushbutton10, 'Enable', 'off')
    set(handles.popupmenu12, 'Enable', 'off')
    set(handles.popupmenu11, 'Enable', 'off')
    set(handles.pushbutton12, 'Enable', 'off')
    set(handles.pushbutton20, 'Enable', 'off') 
    set(handles.popupmenu13, 'Enable', 'off') 
    set(handles.pushbutton30, 'Enable', 'off') 
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'SAE')
    handles.BaseFunc.K_loss          = 'F_kloss'; 
else
    handles.BaseFunc.K_loss          = 'F_kloss_none'; 
end
if handles.Initialized 
    % Disable buttons because the algorithm needs to be initialized again
    set(handles.popupmenu7, 'Enable', 'off')
    set(handles.popupmenu8, 'Enable', 'off')
    set(handles.pushbutton28, 'Enable', 'off') 
    set(handles.pushbutton3, 'Enable', 'off')
    set(handles.popupmenu10, 'Enable', 'off')
    set(handles.pushbutton29, 'Enable', 'off')
    set(handles.pushbutton4, 'Enable', 'off')
    set(handles.checkbox4, 'Enable', 'off')
    set(handles.pushbutton9, 'Enable', 'off')
    set(handles.pushbutton10, 'Enable', 'off')
    set(handles.popupmenu12, 'Enable', 'off')
    set(handles.popupmenu11, 'Enable', 'off')
    set(handles.pushbutton12, 'Enable', 'off')
    set(handles.pushbutton20, 'Enable', 'off') 
    set(handles.popupmenu13, 'Enable', 'off') 
    set(handles.pushbutton30, 'Enable', 'off') 
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uibuttongroup5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.ConstraintOptions = hObject;
% Update handles structure
guidata(hObject, handles);

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.Constraints.C_Wch.Th_lb   	= str2double(get(handles.edit8,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit8 = hObject;
set(handles.edit8,'String','0.4')
handles.Constraints.C_Wch.Th_lb   	= str2double(get(handles.edit8,'String'));
% Update handles structure
guidata(hObject, handles)

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
handles.Constraints.C_PIch.Th_lb   	= str2double(get(handles.edit9,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit9 = hObject;
set(handles.edit9,'String','0.4')
handles.Constraints.C_PIch.Th_lb   	= str2double(get(handles.edit9,'String'));
% Update handles structure
guidata(hObject, handles)

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
handles.Constraints.C_Wzs.Th_lb   	= str2double(get(handles.edit10,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit10 = hObject;
set(handles.edit10,'String','0.4')
handles.Constraints.C_Wzs.Th_lb   	= str2double(get(handles.edit10,'String'));
% Update handles structure
guidata(hObject, handles)

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
handles.Constraints.C_PIzs.Th_lb   	= str2double(get(handles.edit11,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit11 = hObject;
set(handles.edit11,'String','0.4')
handles.Constraints.C_PIzs.Th_lb   	= str2double(get(handles.edit11,'String'));
% Update handles structure
guidata(hObject, handles)

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
handles.Constraints.C_Wch.Th_ub   	= str2double(get(handles.edit12,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit12 = hObject;
set(handles.edit12,'String','0.4')
handles.Constraints.C_Wch.Th_ub   	= str2double(get(handles.edit12,'String'));
% Update handles structure
guidata(hObject, handles)

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
handles.Constraints.C_PIch.Th_ub   	= str2double(get(handles.edit13,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit13 = hObject;
set(handles.edit13,'String','0.4')
handles.Constraints.C_PIch.Th_ub   	= str2double(get(handles.edit13,'String'));
% Update handles structure
guidata(hObject, handles);

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
handles.Constraints.C_Wzs.Th_ub   	= str2double(get(handles.edit14,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit14 = hObject;
set(handles.edit14,'String','0.4')
handles.Constraints.C_Wzs.Th_ub   	= str2double(get(handles.edit14,'String'));
% Update handles structure
guidata(hObject, handles);

function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
handles.Constraints.C_PIzs.Th_ub   	= str2double(get(handles.edit15,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.edit15 = hObject;
set(handles.edit15,'String','0.4')
handles.Constraints.C_PIzs.Th_ub   	= str2double(get(handles.edit15,'String'));
% Update handles structure
guidata(hObject, handles)

% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'LSOptim')
    handles.LSOptim          = 1; 
    % Theshold Constraints tighter if LSOptim is used)
    set(handles.edit8,'String','0.3')
    set(handles.edit9,'String','0.3')
    set(handles.edit10,'String','0.3')
    set(handles.edit11,'String','0.3')
    set(handles.edit12,'String','0.3')
    set(handles.edit13,'String','0.3')
    set(handles.edit14,'String','0.3')
    set(handles.edit15,'String','0.3')
    handles.Constraints.C_Wch.Th_lb    	= str2double(get(handles.edit8,'String'));
    handles.Constraints.C_PIch.Th_lb   	= str2double(get(handles.edit9,'String'));
    handles.Constraints.C_Wzs.Th_lb   	= str2double(get(handles.edit10,'String'));
    handles.Constraints.C_PIzs.Th_lb   	= str2double(get(handles.edit11,'String'));
    handles.Constraints.C_Wch.Th_ub   	= str2double(get(handles.edit12,'String'));
    handles.Constraints.C_PIch.Th_ub   	= str2double(get(handles.edit13,'String'));
    handles.Constraints.C_Wzs.Th_ub   	= str2double(get(handles.edit14,'String'));
    handles.Constraints.C_PIzs.Th_ub   	= str2double(get(handles.edit15,'String'));
else
    handles.LSOptim          = 0; 
  	set(handles.edit8,'String','0.4')
    set(handles.edit9,'String','0.4')
    set(handles.edit10,'String','0.4')
    set(handles.edit11,'String','0.4')
    set(handles.edit12,'String','0.4')
    set(handles.edit13,'String','0.4')
    set(handles.edit14,'String','0.4')
    set(handles.edit15,'String','0.4')
    handles.Constraints.C_Wch.Th_lb    	= str2double(get(handles.edit8,'String'));
    handles.Constraints.C_PIch.Th_lb   	= str2double(get(handles.edit9,'String'));
    handles.Constraints.C_Wzs.Th_lb   	= str2double(get(handles.edit10,'String'));
    handles.Constraints.C_PIzs.Th_lb   	= str2double(get(handles.edit11,'String'));
    handles.Constraints.C_Wch.Th_ub   	= str2double(get(handles.edit12,'String'));
    handles.Constraints.C_PIch.Th_ub   	= str2double(get(handles.edit13,'String'));
    handles.Constraints.C_Wzs.Th_ub   	= str2double(get(handles.edit14,'String'));
    handles.Constraints.C_PIzs.Th_ub   	= str2double(get(handles.edit15,'String'));
end


% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu6 = hObject;
Solvers  = {'lsqnonlin','LSOptim'};
set(handles.popupmenu6,'String',Solvers);
guidata(hObject, handles);


% % --- Executes on key press with focus on pushbutton12 and none of its controls.
% function pushbutton12_KeyPressFcn(hObject, eventdata, handles)
% % hObject    handle to pushbutton12 (see GCBO)
% % eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
% %	Key: name of the key that was pressed, in lower case
% %	Character: character interpretation of the key(s) that was pressed
% %	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% % --- Executes on key press with focus on popupmenu7 and none of its controls.
% function popupmenu7_KeyPressFcn(hObject, eventdata, handles)
% % hObject    handle to popupmenu7 (see GCBO)
% % eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
% %	Key: name of the key that was pressed, in lower case
% %	Character: character interpretation of the key(s) that was pressed
% %	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% % handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
contents = cellstr(get(handles.popupmenu8,'String'));
handles.curr_PlotType = contents{get(handles.popupmenu8,'Value')};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.popupmenu8 = hObject;
handles.curr_PlotType = 'Ellipse Model';
handles.PlotTypes_sh    = {'Ellipse Model','Base Functions'};
handles.PlotTypes       = {'Ellipse Model','Efficiency Model','Base Functions','Eff. Iso Lines','3D plot','Enthalpy gain','Work Coefficient'};
set(handles.popupmenu8,'String',handles.PlotTypes_sh);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu8 = hObject;
set(handles.popupmenu8, 'Enable', 'off') 
guidata(gcbo, handles);
guidata(hObject, handles);

% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select parameter struct to use
% Update current selected PlotType
contents = cellstr(get(handles.popupmenu8,'String'));
handles.curr_PlotType = contents{get(handles.popupmenu8,'Value')};
% Select type of plot
if strcmp(handles.curr_PlotType,'Ellipse Model')
    plot_modelMassFlow(handles.map,handles.Options,handles.Curr_param,handles.ax,handles.Plot_color); 
    if handles.holdPlot
        hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end 
elseif strcmp(handles.curr_PlotType,'Efficiency Model')
    plot_modelEfficiency(handles.map,handles.Options,handles.Curr_param,handles.ax2,handles.Plot_color)
    if handles.holdPlot2
        hold(handles.ax2,'on')
    else
        hold(handles.ax2,'off')
    end
elseif strcmp(handles.curr_PlotType,'Eff. Iso Lines')
    plot_modelContour(handles.map,handles.Options,handles.Curr_param,handles.ax)
    if handles.holdPlot
        hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end
elseif strcmp(handles.curr_PlotType,'3D plot')
    plot_model3D(handles.map,handles.Options,handles.Curr_param,handles.ax)
    if handles.holdPlot
        hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end
elseif strcmp(handles.curr_PlotType,'Enthalpy gain')
    plot_lambda(handles.map,handles.Options,handles.Curr_param,handles.ax2,2)
    if handles.holdPlot2
        hold(handles.ax2,'on')
    else
        hold(handles.ax2,'off')
    end
elseif strcmp(handles.curr_PlotType,'Work Coefficient')
    plot_lambda(handles.map,handles.Options,handles.Curr_param,handles.ax,1)
    if handles.holdPlot
        hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end
elseif strcmp(handles.curr_PlotType,'Base Functions')
    plot_baseFunctions(handles.map,handles.Options,handles.Curr_param)
else
    %Error
end
guidata(hObject, handles);

% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7
handles.popupmenu7 = hObject;
set(handles.popupmenu7,'String',handles.popupmenu_Plotparameters);
contents = cellstr(get(handles.popupmenu7,'String'));
if strcmp(contents{get(handles.popupmenu7,'Value')},'paramEll_i')
    handles.Curr_param = handles.paramEll_i;
    handles.Plot_color = '-g';
    set(handles.popupmenu8,'String',handles.PlotTypes_sh);
    set(handles.popupmenu8,'Value',1);
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramEll')
    handles.Curr_param = handles.paramEll;
    handles.Plot_color = '-r';
    set(handles.popupmenu8,'String',handles.PlotTypes_sh);
    set(handles.popupmenu8,'Value',1);
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramT_i')
    handles.Curr_param = handles.paramT_i;
    handles.Plot_color = '-k';
    set(handles.popupmenu8,'String',handles.PlotTypes);
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramT')
    handles.Plot_color = '-b';
    handles.Curr_param = handles.paramT;
    set(handles.popupmenu8,'String',handles.PlotTypes);
else
    %error
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu7 = hObject;
set(handles.popupmenu7, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.popupmenu7,'String'));
if strcmp(contents{get(handles.popupmenu7,'Value')},'paramEll_i')
    assignin('base','paramEll_i',handles.paramEll_i)
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramEll')
    assignin('base','paramEll',handles.paramEll)
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramT_i')
    assignin('base','paramT_i',handles.paramT_i)
elseif strcmp(contents{get(handles.popupmenu7,'Value')},'paramT')
    assignin('base','paramT',handles.paramT)
else
    %error
end
assignin('base','Options',handles.Options)
assignin('base','map',handles.map)
disp ('Parameters correctly extracted')


% --- Executes on selection change in popupmenu10.
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10
handles.popupmenu10 = hObject;
set(handles.popupmenu10,'String',handles.popupmenu_EllipseInitparameters);
contents = cellstr(get(handles.popupmenu10,'String'));
if strcmp(contents{get(handles.popupmenu10,'Value')},'paramEll_i')
    handles.InitEllipseParam = handles.paramEll_i;
elseif strcmp(contents{get(handles.popupmenu10,'Value')},'paramEll')
    handles.InitEllipseParam = handles.paramEll;
else
    %error
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu10 = hObject;
set(handles.popupmenu10, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11
handles.popupmenu11 = hObject;
set(handles.popupmenu11,'String',handles.popupmenu_TotalInitparameters);
contents = cellstr(get(handles.popupmenu11,'String'));
if strcmp(contents{get(handles.popupmenu11,'Value')},'paramT_i')
    handles.InitTotalParam = handles.paramT_i;
elseif strcmp(contents{get(handles.popupmenu11,'Value')},'paramT')
    handles.InitTotalParam = handles.paramT;
else
    %error
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu11 = hObject;
set(handles.popupmenu11, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12
handles.popupmenu12 = hObject;
set(handles.popupmenu12,'String',handles.popupmenu_EffInitparameters);
contents = cellstr(get(handles.popupmenu12,'String'));
if strcmp(contents{get(handles.popupmenu12,'Value')},'paramEll_i')
    handles.InitEffParam = handles.paramEll_i;
elseif strcmp(contents{get(handles.popupmenu12,'Value')},'paramEll')
    handles.InitEffParam = handles.paramEll;
else
    %error
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.popupmenu12 = hObject;
set(handles.popupmenu12, 'Enable', 'off') 
guidata(gcbo, handles)

% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.popupmenu13,'String'));
handles.curr_PlotTypeRaw = contents{get(handles.popupmenu13,'Value')};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
 
handles.popupmenu13 = hObject;
set(handles.popupmenu13, 'Enable', 'off')
handles.curr_PlotTypeRaw = 'Pressure_Ratio';
PlotTypes  = {'Pressure_Ratio','Efficiency','Enthalphy','Work Coefficient'};
set(handles.popupmenu13,'String',PlotTypes);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select type of plot
if strcmp(handles.curr_PlotTypeRaw,'Pressure_Ratio')
    %set(handles.ax,'ColorOrder',parula(length(handles.map.Nc_vec)));
    plot(handles.ax,handles.map.WcCorr_M,handles.map.PiC_M,'-x','LineWidth',2);
    axis(handles.ax,[0,max(handles.map.WcCorr_M(:))*1.05,0,max(handles.map.PiC_M(:))*1.05])
    grid(handles.ax,'on');
    xlabel(handles.ax,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
    ylabel(handles.ax,'$\Pi_c \quad [-]$','interpreter', 'latex','FontSize',13);
elseif strcmp(handles.curr_PlotTypeRaw,'Efficiency')
    plot(handles.ax2,handles.map.WcCorr_M,handles.map.etaC_M,'-x','LineWidth',2);
    %set(handles.ax2,'ColorOrder',parula(length(handles.map.Nc_vec)));
    axis(handles.ax2,[0,max(handles.map.WcCorr_M(:))*1.05,min(handles.map.etaC_M(:))*0.7,max(handles.map.etaC_M(:))*1.05])
    grid(handles.ax2,'on');
    xlabel(handles.ax2,'$\bar{W}_{c} \quad [kg/s]$','interpreter', 'latex','FontSize',13);
    ylabel(handles.ax2,'$\eta_c \quad [-]$','interpreter', 'latex','FontSize',13);
elseif strcmp(handles.curr_PlotTypeRaw,'Enthalphy')
    plot_lambda(handles.map,handles.Options,[],handles.ax2,2)
    %set(handles.ax2,'ColorOrder',parula(length(handles.map.Nc_vec)));
    if handles.holdPlot2
        hold(handles.ax2,'on')
    else
        hold(handles.ax2,'off')
    end
elseif strcmp(handles.curr_PlotTypeRaw,'Work Coefficient')
    plot_lambda(handles.map,handles.Options,[],handles.ax,1)
    %set(handles.ax,'ColorOrder',parula(length(handles.map.Nc_vec)));
    if handles.holdPlot
        hold(handles.ax,'on')
    else
        hold(handles.ax,'off')
    end
else
    %Error
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu14.
function popupmenu14_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu14


% --- Executes during object creation, after setting all properties.
function popupmenu14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton30 = hObject;
set(handles.pushbutton30, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton10 = hObject;
set(handles.pushbutton10, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton9 = hObject;
set(handles.pushbutton9, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton20 = hObject;
set(handles.pushbutton20, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function checkbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.checkbox4       = hObject;
handles.DragInitGuess   = 0;
set(handles.checkbox4, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton4 = hObject;
set(handles.pushbutton4, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton3 = hObject;
set(handles.pushbutton3, 'Enable', 'off') 
guidata(gcbo, handles)


% --- Executes during object creation, after setting all properties.
function pushbutton12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton12 = hObject;
set(handles.pushbutton12, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton29 = hObject;
set(handles.pushbutton29, 'Enable', 'off') 
guidata(gcbo, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton28 = hObject;
set(handles.pushbutton28, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.pushbutton16 = hObject;
set(handles.pushbutton16, 'Enable', 'off') 
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function checkbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.holdPlot = get(hObject,'Value');
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function checkbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.holdPlot2 = get(hObject,'Value');
guidata(gcbo, handles);

% Subfunctions
function [map,requiredMapData,MapFormat] = CheckMapData(map)
FoundNcCorr     = isfield(map,'NcCorr');
FoundWcCorr     = isfield(map,'WcCorr');
FoundPiC        = isfield(map,'PiC');
FoundetaC       = isfield(map,'etaC');

FoundTCref    	= isfield(map,'TCref');
if ~FoundTCref
    map.TCref   = 298;
    FoundTCref  = 1;
end
FoundpCref      = isfield(map,'pCref');
if ~FoundpCref
    map.pCref   = 1e5;
    FoundpCref  = 1;
end
% FoundMODEL      = isfield(map,'MODEL');
% if ~FoundMODEL
%     map.MODEL = 'Unknown_Compressor';
%     FoundMODEL = 1;
% end
requiredMapData = logical(FoundNcCorr.*FoundWcCorr.*FoundPiC.*FoundetaC.*FoundTCref.*FoundpCref);
if requiredMapData
    % Also Map signals have to be column vector with same size the four of them
    NcCorrSize      = size(map.NcCorr);
    NcCorrOk        = (NcCorrSize(1)>1)&&(NcCorrSize(2)==1);
    WcCorrSize      = size(map.WcCorr);
    WcCorrOk        = (WcCorrSize(1)>1)&&(WcCorrSize(2)==1);
    PiCSize         = size(map.PiC);
    PiCOk           = (PiCSize(1)>1)&&(PiCSize(2)==1);
    etaCSize        = size(map.etaC);
    etaCOk          = (etaCSize(1)>1)&&(etaCSize(2)==1);
    % first value has to be equal, also the second value, to be column vectors.
    lengthEqual     = prod([NcCorrSize(1),NcCorrSize(1),NcCorrSize(1),NcCorrSize(1)] == [NcCorrSize(1),WcCorrSize(1),PiCSize(1),etaCSize(1)]);
    ColumnVectors 	= prod([NcCorrSize(2),NcCorrSize(2),NcCorrSize(2),NcCorrSize(2)] == [NcCorrSize(2),WcCorrSize(2),PiCSize(2),etaCSize(2)]);
    MapFormat       = logical(lengthEqual.*ColumnVectors.*NcCorrOk.*WcCorrOk.*PiCOk.*etaCOk);
else
    MapFormat       = 0;
end

function [grid] = getEllipses(InitGuess)
Nc_vec         	= InitGuess.N_Ch;
points_in_SpL  	= 80;
Nc_vec_n        = InitGuess.N_Ch_n;
Wch             = InitGuess.W_chk;
WZSL            = InitGuess.W_ZSL;
PIch            = InitGuess.PI_chk;
PIZSL           = InitGuess.PI_ZSL;
CUR_vec         = InitGuess.CUR_i;
grid.Wc         = zeros(sum(points_in_SpL),length(Nc_vec));
grid.Nc         = repmat(Nc_vec',points_in_SpL,1);
for i = 1:length(Nc_vec_n)
    Wc_ellipse_1    = linspace(WZSL(i),0.98.*Wch(i),points_in_SpL-ceil(0.2*points_in_SpL))';
    Wc_ellipse_2    = linspace(0.98.*Wch(i),Wch(i),ceil(0.2*points_in_SpL))';
    Wc_ellipse      = [Wc_ellipse_1;Wc_ellipse_2];
    grid.Wc(:,i)    = Wc_ellipse;
end
% Compute Pressure Ratio
Wch         = repmat(Wch',points_in_SpL,1);
PIch        = repmat(PIch',points_in_SpL,1);
WZSL        = repmat(WZSL',points_in_SpL,1);
PIZSL       = repmat(PIZSL',points_in_SpL,1);
CUR         = repmat(CUR_vec',points_in_SpL,1);
grid.PiC    = zeros(size(grid.Wc));
grid.PiC = (((1-(((grid.Wc-WZSL)./(Wch-WZSL)).^CUR))).^(1./CUR)).*(PIZSL-PIch)+PIch;

function getCoord(f,evnt,handles)
drawnow
lH = handles.lH;
click_type = get(f,'SelectionType');
aH  = get(f,'CurrentAxes');
ptH = getappdata(aH,'CurrentPoint');
delete(ptH)
if strcmp(click_type,'normal')
    %Finding the closest point and highlighting it  
    minDist = realmax;
    finalIdx = NaN;
    finalH = NaN;
    pt = get(aH,'CurrentPoint'); %Getting click position
    for ii = 1:length(lH)
        lHcurr = lH(ii);
        xp=get(lHcurr,'Xdata'); %Getting coordinates of line object
        yp=get(lHcurr,'Ydata');
        dx=daspect(aH);      %Aspect ratio is needed to compensate for uneven axis when calculating the distance
        [newDist idx] = min( ((pt(1,1)-xp).*dx(2)).^2 + ((pt(1,2)-yp).*dx(1)).^2 );
        if (newDist < minDist)
            finalH = lHcurr;
            finalIdx = idx;
            minDist = newDist;
        end
    end
    xp=get(finalH,'Xdata'); %Getting coordinates of line object
    yp=get(finalH,'Ydata');
    ptH = plot(aH,xp(finalIdx),yp(finalIdx),'r*','MarkerSize',20);
    setappdata(aH,'SelPoint',ptH);
    setappdata(aH,'ModBase',finalH);
    setappdata(aH,'ModPoint',finalIdx);
elseif strcmp(click_type,'alt')
    disp('Done clicking!');
    uiresume(f);
end

function dragCoord(f,evnt,handles)
drawnow
click_type = get(f,'SelectionType');
aH  = get(f,'CurrentAxes');
ptH = getappdata(aH,'SelPoint'); % Remove selected point
delete(ptH)
if strcmp(click_type,'normal')
	ModBase = getappdata(aH,'ModBase');
    xp = get(ModBase,'Xdata'); %Getting coordinates of line object
    yp = get(ModBase,'Ydata');
	ModPoint = getappdata(aH,'ModPoint');
 	pt  = get(aH,'CurrentPoint'); %Getting click position   
    xp(ModPoint) = pt(1,1);
    yp(ModPoint) = pt(1,2);
    set(ModBase,'Xdata',xp); %Getting coordinates of line object
    set(ModBase,'Ydata',yp); 
    W_chk    = get(handles.lH(1),'Xdata')';
    W_ZSL    = get(handles.lH(2),'Xdata')';
    PI_chk   = get(handles.lH(1),'Ydata')';
    PI_ZSL   = get(handles.lH(2),'Ydata')';
    Curr_InitGuess.W_chk    = W_chk(ModPoint); % Get current ellipse, recompute and replot
    Curr_InitGuess.W_ZSL    = W_ZSL(ModPoint);
    Curr_InitGuess.PI_chk   = PI_chk(ModPoint);
    Curr_InitGuess.PI_ZSL   = PI_ZSL(ModPoint);
	Curr_InitGuess.N_Ch     = handles.InitGuess.N_Ch(ModPoint);
    Curr_InitGuess.N_Ch_n   = handles.InitGuess.N_Ch_n(ModPoint);
    Curr_InitGuess.CUR_i    = handles.InitGuess.CUR_i(ModPoint);
    SpL_grid                = getEllipses(Curr_InitGuess);
	set(handles.Ell(ModPoint),'Xdata',SpL_grid.Wc); %Getting coordinates of line object
    set(handles.Ell(ModPoint),'Ydata',SpL_grid.PiC);
    drawnow
	setappdata(aH,'InGuW_chk',W_chk);
	setappdata(aH,'InGuW_ZSL',W_ZSL);
	setappdata(aH,'InGuPI_chk',PI_chk);
	setappdata(aH,'InGuPI_ZSL',PI_ZSL);
elseif strcmp(click_type,'alt')
    disp('Done clicking!');
    uiresume(f);
end


% --- Executes during object creation, after setting all properties.
function uitoggletool5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uitoggletool5 = hObject;
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function uitoggletool1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uitoggletool1 = hObject;
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function uitoggletool2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uitoggletool2 = hObject;
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function uitoggletool3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uitoggletool3 = hObject;
guidata(gcbo, handles);


% --- Executes during object creation, after setting all properties.
function uitoggletool4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.uitoggletool4 = hObject;
guidata(gcbo, handles);


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.ax)
legend(handles.ax,'hide')
guidata(gcbo, handles);

% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.ax2)
legend(handles.ax2,'hide')
guidata(gcbo, handles);


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Create figure
% scrsz = get(0,'ScreenSize');
%     figure('Position',[1 0 scrsz(3) scrsz(4)],'numbertitle','off','name', 'Base Functions') 
%     set(gcf, 'PaperType', 'A3','PaperOrientation','landscape','PaperPositionMode','auto')
%     axesHandle = gca;
%figure(1)
figure('name', 'Axes 1')
NewAx1 = gca;

copyobj(get(handles.ax,'children'),NewAx1); 
colormap(NewAx1,get(handles.ax,'DefaultAxesColorOrder'))
set(NewAx1,'Xlim',get(handles.ax,'XLim')) 
set(NewAx1,'Ylim',get(handles.ax,'YLim')) 
set(NewAx1,'Zlim',get(handles.ax,'ZLim')) 

set(NewAx1,'XGrid',get(handles.ax,'XGrid')) 
set(NewAx1,'YGrid',get(handles.ax,'YGrid')) 
set(NewAx1,'ZGrid',get(handles.ax,'ZGrid')) 

[az,el] = view(handles.ax); 
view(NewAx1,[az,el]) 

set(NewAx1,'XLabel',copyobj(get(handles.ax,'XLabel'),NewAx1));
set(NewAx1,'YLabel',copyobj(get(handles.ax,'YLabel'),NewAx1));
set(NewAx1,'ZLabel',copyobj(get(handles.ax,'ZLabel'),NewAx1));

% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure('name', 'Axes 2')
NewAx2 = gca;
copyobj(get(handles.ax2,'children'),NewAx2); 
colormap(NewAx2,get(handles.ax2,'DefaultAxesColorOrder'))
set(NewAx2,'Xlim',get(handles.ax2,'XLim')) 
set(NewAx2,'Ylim',get(handles.ax2,'YLim')) 
set(NewAx2,'Zlim',get(handles.ax2,'ZLim')) 

set(NewAx2,'XGrid',get(handles.ax2,'XGrid')) 
set(NewAx2,'YGrid',get(handles.ax2,'YGrid')) 
set(NewAx2,'ZGrid',get(handles.ax2,'ZGrid')) 

[az,el] = view(handles.ax2); 
view(NewAx2,[az,el]) 

set(NewAx2,'XLabel',copyobj(get(handles.ax2,'XLabel'),NewAx2));
set(NewAx2,'YLabel',copyobj(get(handles.ax2,'YLabel'),NewAx2));
set(NewAx2,'ZLabel',copyobj(get(handles.ax2,'ZLabel'),NewAx2));
