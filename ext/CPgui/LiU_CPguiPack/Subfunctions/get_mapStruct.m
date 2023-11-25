function [map]  = get_mapStruct(NcCorr,WcCorr,PiC,etaC,TCref,pCref,D_2,T01,p01)
% [map]  = get_mapStruct(NcCorr,WcCorr,PiC,etaC,TCref,pCref,D_2,T01,p01)
%
% Creates the map structure required by the parameterization GUI as input,
% first six inputs are mandatory
%
%   NcCorr        	- Column vector with the corrected compressor speed
%   WcCorr         	- Column vector with the corrected compressor mass flow
%   PiC         	- Column vector with the compressor pressure ratio
%   etaC         	- Column vector with the compressor efficiency
%   TCref         	- Compressor map reference temperature
%   pCref         	- Compressor map reference pressure
%   D_2         	- Compressor outer impeller diameter
%   T01         	- Column vector with the compressor inlet temperature
%   p01         	- Column vector with the compressor inlet pressure
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
switch nargin
  case 0
end
%% Create map strcuct with the minimum required data in column vectors
map.NcCorr      = NcCorr;
map.WcCorr      = WcCorr;
map.PiC         = PiC;
map.etaC        = etaC;

%% Include extra inputs in the map structure
if exist('TCref','var')
    if ~isempty(TCref)        
        map.TCref	= TCref;
    end
end

if exist('pCref','var')
	if ~isempty(pCref)
        map.pCref	= pCref;
    end       
end 

if exist('D_2','var')
	if ~isempty(D_2)
        map.D_2	= D_2;
    end 
end

if exist('T01','var')
	if ~isempty(T01)
        map.T01	= T01;
    end 
end

if exist('p01','var')
	if ~isempty(p01)
        map.p01	= p01;
    end 
end