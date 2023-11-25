function [k_loss] = F_kloss_none(param,Nc_corr,Wc,map)
% [k_loss] = F_kloss_none(param,Nc_corr,Wc,map)
%
% Computes k_loss values given the parameters and the corrected speed and mass flow values
%
%   param        	- Parameters of the subfunction c
%   Nc_corr         - Scalar or vector containing Speed values
%   Wc              - Scalar or vector containing Mass flow values
%   map             - Struct containing the map signals
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
%% Return ones of the same size of the input vectors
k_loss  = ones(size(Wc));
end