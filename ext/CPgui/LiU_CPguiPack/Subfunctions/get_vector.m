function vec = get_vector(M)
%vec = get_vector(M)
%
% Returns a vector of the mean of each column in the Matrix M
%
%   M               - Matrix with n columns
%   vec         	- Column vector of length n with the mean value of each M column
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
%% Compute map signals
[~,col] = size(M);
vec = zeros(col,1);
for jj = 1:col
 	extract = M(:,jj);
    vec(jj) = mean(extract(isnan(extract)<1));
end