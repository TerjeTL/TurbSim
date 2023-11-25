%% initLiU_CPgui - Initializes the path for the LiU CPgui
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

%% Find where LiU CPgui package is located
found_InitPath = ~isempty(which('init_LiU_CPgui'));
if found_InitPath
    TgFolder = pwd;
 	cd(TgFolder)
    %addpath(genpath(TgFolder))
    addpath(TgFolder)
    addpath(genpath([TgFolder,'/Subfunctions']))
else
    message = sprintf('The initialization script cannot locate and add the folder LiU_CPguiPack to the path. Try to run the script from the mentioned folder or add it manually to the path.');
	questdlg(message, 'Toolbox missing', 'OK', 'OK');
    return;
end
% Clear variables
clear TgFolder 