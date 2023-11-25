%% exampleLiU_CPgui - Template to initialize, load and run the LiU CPgui with the selected Compressor Map
%       Version 1.2: 2017-12-20
%       Copyright (C) 2017, Xavier Llamas
%
%       This file is part of LiU CPgui.
%
%       LiU CPgui is free software: you can redistribute it and/or modify
%       it under the terms of the GNU Lesser General Public License as
%       published by the Free Software Foundation, version 3 of the License.
%
%       This package is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%       GNU Lesser General Public License for more details.
%
%       You should have received a copy of the GNU Lesser General Public License
%       along with LiU CPgui.  If not, see <http://www.gnu.org/licenses/>.

%% Add required folders to path
init_LiU_CPgui
if ~found_InitPath
    return;
end
%% Load here your compressor map signals
load ExampleMap  % (ExampleMap is not a real measured compressor map)
% Required signals:
NcCorr      	= Nc;           % [rpm]     Column vector with the corrected compressor speed
WcCorr       	= Wc;           % [kg/s]    Column vector with the corrected compressor mass flow
PiC            	= PI;           % [-]       Column vector with the compressor pressure ratio
etaC         	= eff;          % [-]       Column vector with the compressor efficiency
TCref        	= Treference;   % [K]       Compressor map reference temperature, scalar
pCref       	= preference;   % [Pa]      Compressor map reference pressure, scalar
% Optional signals:
D_2           	= diameter;    	% [m]       Compressor outer impeller diameter, scalar
% T01        	= [];         	% [K]       Column vector with the compressor inlet temperature
% p01         	= [];           % [Pa]   	Column vector with the compressor inlet pressure
%% The provided get_mapStruct, collects the compressor map signals in the map struct that the GUI requires as input.
% Example call:
% [map]  = get_mapStruct(NcCorr,WcCorr,PiC,etaC,TCref,pCref,dC2,T01,p01)
map          	= get_mapStruct(NcCorr,WcCorr,PiC,etaC,TCref,pCref,D_2);
%% Call the Compressor Parameterization Package GUI with the created compressor map struct
LiU_CPgui(map)