function [map]  = get_signals(map,k_q_in)
%[map]  = get_signals(map,k_q_in)
%
% Creates all initial signals required for the given compressor map, in the
% given map struct
%
%   map             - Struct containing the map signals
%   k_q_in         	- Heat transfer correction parameter
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
manual_Cp       = 0;
manual_Kfric  	= 0;
manual_b2       = 0;

switch nargin
  case 0
        % Error;
  case 1
    manual_Cp       = 0;
    manual_Kfric  	= 0;
    manual_b2       = 0;
  case 2
  otherwise
    manual_Cp       = 1;
    manual_Kfric  	= 1;
    manual_b2       = 1;
end
dC2Bool             = 1;
map                 = get_cGeometry(map);
map                 = get_matrices(map);
map                 = get_gammaC(map);
map                 = get_cpC(map);
map                 = get_T01(map);
map                 = get_p01(map);
map                 = get_rhoC1(map);
map.Nc_max_map      = max(map.NcCorr);
map.Wc_max_map      = max(map.WcCorr);
map.PI_max_map      = max(map.PiC);
map.Nc_vec          = get_vector(map.NcCorr_M);
map.Nc_vec_n        = map.Nc_vec./map.Nc_max_map;
map.n_SpL           = length(map.Nc_vec);
map.n_measured      = length(map.WcCorr);
map.meas_SpL        = zeros(map.n_SpL,1);
map.gamma_air       = mean(map.gammaC);
map.Cp_air          = mean(map.cpC);
map.R_air           = mean(map.Rc);
map.T_1             = mean(map.T01);
map.rho1            = mean(map.p01)./(map.R_air.*map.T_1);
map.T01_M           = write_matrix(map.T01,map);
map.p01_M           = write_matrix(map.p01,map);
map.T_ref           = map.TCref;
map.p_ref           = map.pCref;
for i = 1:map.n_SpL 
    map.meas_SpL(i) = sum(isnan(map.WcCorr_M(:,i))<1);
end
% Check if Diameter is provided with the map, otherwise provide a Guess
if ~isfield(map,'D_2')
    U_2_max = 516.1; %[m/s] rule of thumb from SAE Technical paper 2016-01-1037
    map.D_2 = U_2_max.*60./map.Nc_max_map./pi;
    dC2Bool = 1;
end
map                 = get_Ucm(map);
map                 = get_lambda(map);
map                 = get_PhiCm(map);
map                 = get_PsiC(map);
map                 = get_MachC(map);
map.PhiC_M          = write_matrix(map.PhiC,map);
map.PsiC_M          = write_matrix(map.PsiC,map);
map.hoC_M           = write_matrix(map.hoC,map);
map.hosC_M          = write_matrix(map.hosC,map);
map.lambda_M        = write_matrix(map.lambda,map);
map.MachC_M         = write_matrix(map.MachC,map);
map.MachC_vec       = get_vector(map.MachC_M);
map.Delta_h_max     = max(map.hoC_M(:));
%% Correction method with constant kp?
if isfield(map,'etaC_heat') % restart efficiency if it has been already corrected before applying the heat correction again.
    map.etaC	= map.etaC_heat;  
else
    map.etaC    = map.etaC;
end

if k_q_in == 0
else
    if dC2Bool>0
        %% Correct Efficiency
        % Parameters for correction
        if manual_Cp;       Cp_diff = Cp_diff_in; 	else   Cp_diff	= 0.55; % (From Casey and Schlegel 2009)     
        end;
        if manual_Kfric;    k_fric  = k_fric_in;  	else   k_fric 	= 0.004; % (From Casey and Schlegel 2009)   
        end;
        if manual_b2;       b2  = b2_in;  	elseif isfield(map,'b2'); b2 = map.b2;   else b2      = 0.1*map.D_2;    end;
        % Correct Efficiency
        % Parameters for correction
        % Get and reshape signals
        map                     = get_lambda_euler(map,k_fric);
        map.lambda_euler_M      = write_matrix(map.lambda_euler,map);
        map                     = iterate_rhoc2m(map,k_fric,Cp_diff,b2);
        map.PhiC2_M             = write_matrix(map.PhiC2,map);

        if length(k_q_in)>1 % If there are different k_q for each speed line
            map                     = get_efficiency_correction(map,k_fric,k_q_in);
        else
            k_q                     = max(0,k_q_in);
            map                     = get_adiabatic_efficiency(map,k_fric,k_q);
        end

        map.etaC_heat           = map.etaC;
        map.etaC                = map.etaC_ad;
        % Update Signals
        map             = get_matrices(map);
        map             = get_gammaC(map);
        map          	= get_cpC(map);
        map           	= get_T01(map);
        map          	= get_p01(map);
        map          	= get_rhoC1(map);
        map             = get_Ucm(map);
        map             = get_lambda(map);
        map             = get_PhiCm(map);
        map             = get_PsiC(map);
        map             = get_MachC(map);
        map.PhiC_M  	= write_matrix(map.PhiC,map);
        map.PsiC_M    	= write_matrix(map.PsiC,map);
        map.hoC_M      	= write_matrix(map.hoC,map);
        map.hosC_M    	= write_matrix(map.hosC,map);
        map.lambda_M   	= write_matrix(map.lambda,map);
        map.MachC_M    	= write_matrix(map.MachC,map);
        map.Delta_h_max	= max(map.hoC_M(:));
    end
end
end

%% Subfunctions
function [var_M] = write_matrix(var,map)
% Separate per SpL
diff = max(map.NcCorr).*0.01;
[allN, sortInd, NcCorrfound] = get_SpL(map,diff);
if ~NcCorrfound 
  disp('Not possible to write in matrix format. No NcCorr'); return; 
end
cols    = length(allN);
rows    = length(sortInd{1});
for ii=2:length(sortInd),
  rows = max(length(sortInd{ii}),rows);
end
var_M   = NaN*ones(rows,cols);
%% Rewrite the vector variable into matrix format
for ii=1:cols,
  var_M(1:length(sortInd{ii}),ii)   = var(sortInd{ii})';
end
end

function [map] = iterate_rhoc2m(map,k_fric,Cp_diff,b2)
% switch nargin
%   case 0,
%     clc; clear all;
%     load mapDB;
%     map = allMAPS{22};
%     doPlot = 1;
%   otherwise
%     doPlot = 0;
% end
% get lambda euler
[map, lambdaFound]  = get_lambda(map);
[map, PhiFound]     = get_PhiCm(map);
[map, rhoC1found]   = get_rhoC1(map);
[map, Ucfound]      = get_Ucm(map);
lambda_euler = map.lambda./(1+k_fric./map.PhiC);
% Initial guess for rho
rho2 = map.rhoC.*map.PiC;
iter_max = 200;
% compute cm2
D2          = map.D_2;
%b2          = 0.1*D2; % Assumption
%b2          = map.b2;
m_dot       = map.Wc;
threshold   = 1e-6;
delta_rho   = 10.*ones(length(map.PiC),1);
u2          = map.Uc;
pt3         = map.p01.*map.PiC; 
cp          = map.cpC;
Tt3         = map.T01.*(1+(map.PiC.^((map.gammaC-1)./map.gammaC)-1)./map.etaC);

cm2         = zeros(length(map.PiC),1);
cu2         = zeros(length(map.PiC),1);
c2          = zeros(length(map.PiC),1);
p2          = zeros(length(map.PiC),1);
T2          = zeros(length(map.PiC),1);
rho2_new    = zeros(length(map.PiC),1);
rho2_ok     = zeros(length(map.PiC),1);
Cp_diff_run = Cp_diff.*ones(length(map.PiC),1);
for i=1:length(map.PiC)
    iter2 = 0;
    while rho2_ok(i)<1
        iter = 0;
        while (delta_rho(i)>threshold)
            cm2(i) = m_dot(i)./(pi.*D2.*b2.*rho2(i));
            cu2(i) = lambda_euler(i).*u2(i);
            c2(i) = sqrt(cm2(i).^2+cu2(i).^2);
            p2(i) = pt3(i)-0.5.*Cp_diff_run(i).*rho2(i).*c2(i).^2;
            T2(i) = Tt3(i)-(c2(i).^2)./(2.*cp(i));
            rho2_new(i) = p2(i)./(map.Rc(i).*T2(i));
            delta_rho(i) = abs(rho2_new(i)-rho2(i));
            rho2(i) =rho2_new(i);
            iter = iter+1;
            if iter>iter_max
                disp('Max iteration reached in rho2')
                return
            end
        end
        if rho2(i)<map.rhoC(i)
            % Might happen at high flows that the given Cp_diff gives a too low outlet density, in that case, Cp_diff is lowered and the iteration starts again
            Cp_diff_run(i) = Cp_diff_run(i)-0.02; 
            delta_rho(i) = 10;
        else
            rho2_ok(i) = 1;
        end
        iter2 = iter2+1;
        if iter2>iter_max
            disp('Max iteration reached in rho2')
            return
        end
    end
end
map.rhoC2   = rho2;
map.T2      = T2;
map.p2      = p2;
map.cu2   	= cu2;
map.cm2   	= cm2;
map.c2      = c2;
map.PhiC2	= map.PhiC.*(map.rhoC.*map.D_2)./(pi.*b2.*map.rhoC2);
% switch nargout
%   case 0,
%     VcCorrFound
%     etaCfound
%     PiCfound
%     T01found
%     gammaCfound
%     VcCorr2Found
%   otherwise
end

function [map, Ucfound] = get_Ucm(map, dC)
switch nargin
  case 0
  case 2 % Assume get_Uc(N_tc,dC)
    clear tempMap;
    tempMap.Ntc = map;
    tempMap.dC = dC;
    map = tempMap;
  otherwise
end
f_Uc = @(N, d)(N/60 * pi*d  );
Ucfound = 1;

if isfield(map,'Uc')
elseif  isfield(map,'Ntc') && isfield(map,'dC')
  map.Uc = f_Uc(map.Ntc, map.D_2);
else
  [map, Ntcfound] = get_Ntc(map);
  [map, cGeomFound] = get_cGeometry(map);
  if Ntcfound && cGeomFound
    map.Uc = f_Uc(map.Ntc, map.D_2);    
  else
    Ucfound = 0;
  end
end

switch nargout
  case 0    
    map.Uc
  otherwise
end
end

function [map, T01found] = get_T01(map)

switch nargin
  case 0
    clc; clear all;
    load mapDB;
    map = allMAPS{30};
  otherwise
end

T01Ass = 298;
T01found = 1;
if isfield(map,'T01')
  
elseif isfield(map,'TCref') && isfield(map,'PiC')
  map.T01 = map.TCref*ones(size(map.PiC));   
elseif isfield(map,'PiC')
  map.T01 = T01Ass*ones(size(map.PiC));   
else
  T01found = 0;  
end
end

function [allNcCorr, sortInd, NcCorrFound] = get_SpL(map, dif)
% Modified function from get_allNcCorr (Oskar Leufvén)
switch nargin
  case 0
    dif = 1e3;
  case 1
    dif = .99e3;
  otherwise
end

diffNabs = dif;
  
allNcCorr = [];
sortInd = [];
NcCorrFound = 1;

%% Check if speed and mass flow are available
[map, NcCorrFound] = get_NcCorr(map);
if ~NcCorrFound, NcCorrFound = 0; return, end;
[map, WcCorrfound] = get_WcCorr(map);
if ~WcCorrfound, NcCorrFound = 0; return, end;


if isfield(map,'allNtc') % Even more special field! =) Since it is commonly Ntc that is measured.
  % Sort in Ntc, and store away allNcCorr later
  dif = 2e3;
  allNtc = map.allNtc;
  sortInd = cell(length(allNtc),1);
  TO1mean = NaN*ones(length(allNtc),1);
  for kkk = 1:length(allNtc)
    iii = find(abs(map.Ntc - allNtc(kkk)) < dif);
    if ~isempty(iii)
      sortInd(kkk,:) = {iii'};
      T01mean(kkk,:) = mean(map.T01(iii));
    end
  end
  for kkk = 1:length(allNtc)
    sortedIndices = sortInd{kkk};  
    PAct = map.PiC(sortedIndices);
    WAct = map.WcCorr(sortedIndices);
    theta = atan(PAct./WAct);
    [foo, iii] = sort(theta,'descend');  
    sortInd{kkk} = sortedIndices(iii);
  end
  allNcCorr = allNtc./(sqrt(T01mean/298));

elseif isfield(map,'allNcCorr') % Special field added to struct for especially difficult maps, where the different desired speeds are extracted separetly  
  dif = min(diff(map.allNcCorr)).*0.5; % Sharp limit
  allNcCorr = map.allNcCorr;  
  allNcCorrTemp = allNcCorr;
  sortIndTemp  = cell(length(allNcCorr),1);
  
  for kkk = 1:length(allNcCorr)
    iii = find(abs(map.NcCorr - allNcCorr(kkk)) < dif); 
    if ~isempty(iii)
      sortIndTemp(kkk,:) = {iii'};      
    else
      sortIndTemp(kkk,:) = {NaN};
      allNcCorrTemp(kkk) = NaN;
    end
  end
  allNcCorrTemp = allNcCorrTemp(~isnan(allNcCorrTemp));
  
  sortInd = cell(length(allNcCorrTemp),1);
  counter = 1;
  for kkk = 1:length(allNcCorr)
    iii = find(abs(map.NcCorr - allNcCorr(kkk)) < dif);     
    if ~isempty(iii)
      sortInd(counter,:) = sortIndTemp(kkk,:);      
      counter = counter+1;
    end    
  end
  
  allNcCorr = allNcCorrTemp;
  
  for kkk = 1:length(allNcCorr)
    sortedIndices = sortInd{kkk};  
    PAct = map.PiC(sortedIndices);
    WAct = map.WcCorr(sortedIndices);
    theta = atan(PAct./WAct);
    [foo, iii] = sort(theta,'descend');  
    sortInd{kkk} = sortedIndices(iii);
    
  end
  
  
  
else
  allNcCorr = unique(floor(map.NcCorr/diffNabs))*diffNabs + diffNabs/4;
  
  for kkk = 1:length(allNcCorr)
    iii = find(abs(allNcCorr(kkk) - map.NcCorr) < diffNabs);
    allNcCorr(kkk) = mean(map.NcCorr(iii));
  end
  
  allNcCorr = unique(floor(allNcCorr/diffNabs))*diffNabs + diffNabs/2;
  
  for kkk = 1:length(allNcCorr)
    iii = find(abs(allNcCorr(kkk) - map.NcCorr) < diffNabs);
    allNcCorr(kkk) = mean(map.NcCorr(iii));
  end
  
  allNcCorr = round(allNcCorr/diffNabs)*diffNabs;
  allNcCorr = unique(allNcCorr);
  allNcCorr = unique(floor(allNcCorr/(diffNabs/2)))*diffNabs/2;
  
  for kkk = 1:length(allNcCorr)
    iii = find(abs(allNcCorr(kkk) - map.NcCorr) < diffNabs);
    allNcCorr(kkk) = mean(map.NcCorr(iii));
    PAct = map.PiC(iii);
    WAct = map.WcCorr(iii);
    theta = atan(PAct./WAct);
    [foo, jjj] = sort(theta,'descend'); 
    %[foo, jjj] = sort(map.WcCorr(iii),'ascend');
    sortInd{kkk} = iii(jjj);
  end
  
  sortInd = sortInd';
  allNcCorr = round(allNcCorr/(diffNabs/2))*(diffNabs/2);
  
  maxIteration = 10;
  actIteration = 0;
  absDiff = abs(allNcCorr(1:end-1) - allNcCorr(2:end));
  [minAbsDiff, iii] = min(absDiff);
  while minAbsDiff < diffNabs & actIteration < maxIteration
    actIteration = actIteration + 1;
    
    lowerIndicies  = find(abs(allNcCorr(iii)- map.NcCorr) < diffNabs);
    higherIndicies = find(abs(allNcCorr(iii+1)- map.NcCorr) < diffNabs);
    allNcCorr(iii) = mean([map.NcCorr(lowerIndicies); map.NcCorr(higherIndicies)]);
    allNcCorr(iii+1) = mean([map.NcCorr(lowerIndicies); map.NcCorr(higherIndicies)]);
    
    lowerIndicies  = find(abs(allNcCorr(iii)- map.NcCorr) < diffNabs);
    higherIndicies = find(abs(allNcCorr(iii+1)- map.NcCorr) < diffNabs);
    allNcCorr(iii) = mean([map.NcCorr(lowerIndicies); map.NcCorr(higherIndicies)]);
    allNcCorr(iii+1) = mean([map.NcCorr(lowerIndicies); map.NcCorr(higherIndicies)]);
    
    allNcCorr = unique(allNcCorr);
    absDiff = abs(allNcCorr(1:end-1) - allNcCorr(2:end));
    [minAbsDiff, iii] = min(absDiff);
  end

end


%allNcCorr = round(allNcCorr/1e2)*1e2;


switch nargout
    case 0
    otherwise
end
end

function [map, rhoC1found] = get_rhoC1(map)
switch nargin
  case 0
  otherwise
end


f_rho = @(p01, T01, Rair)(p01./Rair./T01);

rhoC1found = 1;

[map, p01found]  = get_p01(map);
[map, T01found]  = get_T01(map);
[map, Rcfound] = get_Rc(map);


if p01found && T01found && Rcfound
  map.rhoC = f_rho(map.p01, map.T01, map.Rc);
else
  rhoC1found = 0;
end

switch nargout
  case 0
  otherwise
end
end

function [map, RcFound] = get_Rc(map)
% http://en.wikipedia.org/wiki/Gas_constant @121008
%
%          R = 8.314462175 J / K / mol
% M_dry_air  = 28.97;
%
% R_specific = R / M  for dry air: 287.04 J/kg/K
%
% http://www.tis-gdv.de/tis_e/misc/klima.htm
%
% R = cp - cv, gamma = cp/cv
% cp = R + cv = R + cp/gamma => cp (1 - 1/gamma) = R
%
%%
switch nargin
  case 0
  otherwise
end
RairAssumed = 287.04; % Dry gas!
RcFound = 1;
%%

if isfield(map,'Rc')
  return;

elseif isfield(map,'gammaC') && isfield(map,'cpC')
  map.Rc = map.cpC .* (1 - 1./map.gammaC);  
  
% elseif isfield(map,'T01')
%   map = get_gammaC(map);
%   map = get_cpC(map);
%   map = get_Rc(map);
  
elseif isfield(map,'Ntc')
  map.Rc = RairAssumed*ones(size(map.Ntc));
elseif isfield(map,'NcCorr')
  map.Rc = RairAssumed*ones(size(map.NcCorr));  
elseif isfield(map,'VcCorr')
  map.Rc = RairAssumed*ones(size(map.VcCorr));  
else
  RcFound = 0;
end
end

function [map, PsiFound] = get_PsiC(map)
switch nargin
  case 0
  otherwise
end

% 1991:Jensen
f_Psi = @(cpC, T01, PiC, gammaC, Uc)( ( cpC .* T01 .* ( PiC .^ ((gammaC-1)./gammaC) -1 ) ) ./ (0.5.*Uc.^2)  );

PsiFound = 1;

[map, PiCfound]    = get_PiC(map);
[map, T01found]    = get_T01(map);
[map, cpCfound]    = get_cpC(map);
[map, gammaCfound] = get_gammaC(map);
[map, Ucfound]     = get_Ucm(map);


if PiCfound && T01found && cpCfound && gammaCfound && Ucfound
  map.PsiC = f_Psi(map.cpC, map.T01, map.PiC, map.gammaC, map.Uc);
else
  PsiFound = 0;
end



switch nargout
  case 0
  otherwise
end
end

function [map, PiCfound] = get_PiC(map)

switch nargin
  case 0
  otherwise
end

PiCfound = 1;

p01Ass = 100e3; % Default ASSumed compressor inlet pressure

if isfield(map,'PiC')
  
elseif isfield(map,'p02')
  if isfield(map,'p01')
    map.PiC = map.p02./map.p01;
  else
    map.Pic = map.p02./p01Ass;
  end
else
  PiCfound = 0;
  %disp('Not possible to calculate PiC');    
end
end

function [map, PhiFound] = get_PhiCm(map)
switch nargin
  case 0
    doPlot = 1;
  otherwise
    doPlot = 0;
end


f_Phi = @(Wc, rhoC, dC, Uc)(Wc./rhoC./dC^2./Uc);

PhiFound = 1;
T01Ass = 298;
p01Ass = 100e3;
RairAss = 287;

[map, Wcfound]    = get_Wc(map);
[map, rhoC1found] = get_rhoC1(map);
[map, cGeomFound] = get_cGeometry(map);
[map, Ucfound]    = get_Ucm(map);



if Wcfound && Ucfound && cGeomFound && rhoC1found
  map.PhiC = f_Phi(map.Wc, map.rhoC, map.D_2, map.Uc);
elseif doPlot  
  Wcfound 
  Ucfound 
  cGeomFound 
  rhoC1found
  PhiFound = 0;
end

switch nargout
  case 0
  otherwise
end
end

function [map, p01found] = get_p01(map)

switch nargin
  case 0
  otherwise
end

p01Ass = 100e3;

p01found = 1;
if isfield(map,'p01')
  
elseif isfield(map,{'PiC','p02'})
  map.p01 = map.p02./map.PiC;  
elseif isfield(map,'PiC')
  map.p01 = p01Ass*ones(size(map.PiC));  
  p01found = 1;
end
end

function [map, Ntcfound] = get_Ntc(map)

switch nargin
  case 0
  otherwise
end

Ntcfound = 1;
T01Ass = 298;
T03Ass = 873;
TCrefAss = 298;

if isfield(map,'Ntc')
  
elseif isfield(map,{'NcCorr', 'T01', 'TCref'})
  map.Ntc = map.NcCorr.*sqrt(map.T01./map.TCref);
  
elseif isfield(map,{'NcCorr', 'TCref'})
  map.Ntc = map.NcCorr.*sqrt(T01Ass./map.TCref);

elseif isfield(map,{'NcCorr'})
  map.Ntc = map.NcCorr*sqrt(T01Ass/TCrefAss);

elseif isfield(map,{'Uc', 'D_2'})
  map.Ntc = map.Uc/pi/map.D_2*60;% Uc = Ntc/60*pi*D_2
  
elseif isfield(map,{'UcCorr', 'D_2', 'TCref', 'T01'})
  map.Ntc = map.UcCorr ./ sqrt(map.T01 ./ map.TCref) / pi / map.D_2 * 60;% UcCorr = Ntc / sqrt(T01/TCref) / 60 * pi * D_2
  
elseif isfield(map,{'UcCorr', 'D_2', 'TCref'})
  map.Ntc = map.UcCorr ./ sqrt(T01Ass ./ map.TCref) / pi / map.D_2 * 60;% UcCorr = Ntc / sqrt(T01/TCref) / 60 * pi * D_2
  
elseif isfield(map,{'UcCorr', 'D_2', 'T01'})
  map.Ntc = map.UcCorr ./ sqrt(map.T01 ./ TCrefAss) / pi / map.D_2 * 60;% UcCorr = Ntc / sqrt(T01/TCref) / 60 * pi * D_2  

elseif isfield(map,{'UcCorr', 'D_2'})
  map.Ntc = map.UcCorr ./ sqrt(T01Ass ./ TCrefAss) / pi / map.D_2 * 60;% UcCorr = Ntc / sqrt(T01/TCref) / 60 * pi * D_2    
  

elseif isfield(map,'TSP') 
  if isfield(map,'T03')
    map.Ntc = map.TSP*sqrt(map.T03);
  else
    map.Ntc = map.TSP.*sqrt(T03Ass);
  end

elseif isfield(map,{'NtCorr', 'T03', 'TTref'})
  map.Ntc = map.NtCorr.*sqrt(map.T03./map.TTref);
  
elseif isfield(map,{'Ut', 'dT1'}) 
  map.Ntc = map.Ut / pi / map.dT1 * 60;
else
  Ntcfound = 0;  
end

switch nargout
  case 0
  otherwise
end
end

function [map, NcCorrFound] = get_NcCorr(map)

switch nargin
  case 0
  otherwise    
end

NcCorrFound = 1;
TCrefAss = 298;   % Default ASSumed compressor map reference temperature
T01Ass = 298;

f_NcCorr = @(Ntc, T01, TCref)(Ntc./sqrt(T01./TCref));

[map, Ntcfound] = get_Ntc(map);

if isfield(map, 'NcCorr') 
  
elseif Ntcfound
  if isfield(map,'T01')
    if isfield(map,'TCref');
      map.NcCorr = f_NcCorr(map.Ntc, map.T01, map.TCref);
    else      
      map.NcCorr = f_NcCorr(map.Ntc, map.T01, TCrefAss);
    end
  else
    if isfield(map,'TCref');
      map.NcCorr = f_NcCorr(map.Ntc, T01Ass, map.TCref);
    else      
      map.NcCorr = f_NcCorr(map.Ntc, T01Ass, TCrefAss);
    end    
  end
else
  map;
  NcCorrFound = 0;
  %disp('No NcCorr could be found');
end
end

function [map] = get_matrices(map,diff_input)

switch nargin
  case 0
  case 1    
    manual_diff = 0;
  case 2    
    manual_diff = 1;
  otherwise    
end
%% Handle Variable geometry compressors
if isfield(map, 'posVGC')
  title('VGC map');
  uniqueVGCpos = unique(map.posVGC);
  for kkk = 1:length(uniqueVGCpos)
    iii = find(map.posVGC == uniqueVGCpos(kkk));
    tempMap = extractSubMap(map, iii);    
    tempMap = rmfield(tempMap, 'posVGC');
    plot_compMapFlow(tempMap, FlowLineStyle, plotHandle,infoPlot);    
  end  
  return
end
%% Check if the variables are available
[map, PiCfound] = get_PiC(map); 
if ~PiCfound
  disp('Not possible to plot compressor flow map. No PiC'); return; 
end
[map, NcCorrfound] = get_NcCorr(map);
if ~NcCorrfound 
  disp('Not possible to plot compressor flow map. No NcCorr'); return; 
end
if manual_diff
    diff = max(map.NcCorr).*diff_input;
else
    diff = max(map.NcCorr).*0.01;
end
[allN, sortInd, NcCorrfound] = get_SpL(map,diff);
%[allN, sortInd, map, NcCorrfound] = get_allNcCorr(map);
if ~NcCorrfound 
  disp('Not possible to plot compressor flow map. No NcCorr'); return; 
end
[map, WcCorrfound] = get_WcCorr(map);
if ~WcCorrfound
  disp('Not possible to plot compressor flow map. No WcCorr'); return; 
end
cols    = length(allN);
rows    = length(sortInd{1});
for ii=2:length(sortInd)
  rows = max(length(sortInd{ii}),rows);
end
NaNMap  = NaN*ones(rows,cols);
PiC    	= NaNMap;
WcCorr  = NaNMap;
etaC    = NaNMap;
NcCorr  = NaNMap;
%% Fill in the matrices
for ii=1:cols
  PiC(1:length(sortInd{ii}),ii)     = map.PiC(sortInd{ii})';
  etaC(1:length(sortInd{ii}),ii)    = map.etaC(sortInd{ii})';
  WcCorr(1:length(sortInd{ii}),ii)  = map.WcCorr(sortInd{ii})';
  NcCorr(1:length(sortInd{ii}),ii)  = map.NcCorr(sortInd{ii})';
end
% Output the matrices
map.PiC_M       = PiC;
etaC(etaC==0)   = NaN;
map.etaC_M      = etaC;
map.WcCorr_M    = WcCorr;
map.NcCorr_M    = NcCorr;
end

function [map, Machfound] = get_MachC(map)

switch nargin
  case 0
  case 2 % Assume get_Uc(N_tc,dC)
    clear tempMap;
    tempMap.Ntc = map;
    tempMap.dC = dC;
    map = tempMap;
  otherwise
end

f_Mach = @(Uc,gammaC,Rair,T01)(Uc./sqrt(gammaC.*Rair.*T01));
Machfound = 1;

if isfield(map,'MachC')
  
elseif  isfield(map,'Uc') && isfield(map,'gammaC') && isfield(map,'Rair') && isfield(map,'T01')
  map.MachC = f_Mach(map.Uc, map.gammaC, map.Rair, map.T01);
  
else
  [map, RcFound] = get_Rc(map);
  [map, Ucfound] = get_Ucm(map);
  [map, gammaCfound]  = get_gammaC(map);
  [map, T01found]    = get_T01(map);
  if Ucfound && gammaCfound && RcFound && T01found
    map.MachC = f_Mach(map.Uc, map.gammaC, map.Rc, map.T01);    
  else
    Machfound = 0;
  end
end

switch nargout
  case 0
  otherwise
end
end

function [map,lambdaEulerFound] = get_lambda_euler(map,k_fric)
switch nargin
  case 0
    k_fric = 0.004;
  case 1
    k_fric = 0.004;
  otherwise
end
%  Casey
f_lambda_euler      = @(lambda,K,PhiC)(lambda./(1+K./PhiC));
lambdaEulerFound    = 1;
[map, lambdaFound]  = get_lambda(map);
[map, PhiFound]     = get_PhiCm(map);

if PhiFound && lambdaFound
  map.lambda_euler    = f_lambda_euler(map.lambda, k_fric, map.PhiC);
else
  lambdaEulerFound = 0;
end


switch nargout
  case 0
    PiCfound
    T01found
    cpCfound    
    gammaCfound
    Ucfound
  otherwise
end
end

function [map, lambdaFound] = get_lambda(map)
switch nargin
  case 0
  otherwise
end

% Dixon // Casey
f_hosC      = @(cpC, T01, PiC, gammaC)( ( cpC .* T01 .* ( PiC .^ ((gammaC-1)./gammaC) -1 ) ) );
f_hoC       = @(hosC, etaC)(hosC./etaC);
f_lambda    = @(hoC,Uc)(hoC./Uc.^2);

lambdaFound = 1;

[map, PiCfound]     = get_PiC(map);
[map, T01found]     = get_T01(map);
[map, cpCfound]     = get_cpC(map);
[map, gammaCfound]  = get_gammaC(map);
[map, Ucfound]      = get_Ucm(map);
[map, etaCfound]    = get_etaC(map);

if PiCfound && T01found && cpCfound && gammaCfound && Ucfound && etaCfound
  map.hosC      = f_hosC(map.cpC, map.T01, map.PiC, map.gammaC);
  map.hoC       = f_hoC(map.hosC, map.etaC);
  map.lambda    = f_lambda(map.hoC, map.Uc);
else
  lambdaFound = 0;
end



switch nargout
  case 0
  otherwise
end
end

function [map, gammaCfound] = get_gammaC(map)
% Interpolation to get gammaC from T01

% Data from Swiss Auto Wenko-maps
Tint     =  [296.340	297.906	299.472	301.038	302.604	304.170	305.456	306.742	308.028	309.314	310.600	312.826	315.052	317.278	319.504	321.730	324.020	326.310	328.600	330.890	333.180	334.586	335.992	337.398	338.804	340.210	342.032	343.854	345.676	347.498	349.320	350.514	351.708	352.902	354.096	355.290	357.414	359.538	361.662	363.786	365.910	367.972	370.034	372.096	374.158	376.220	378.164	380.108	382.052	383.996	385.940	387.216	388.492	389.768	391.044	392.320	393.688	395.056	396.424	397.792	399.160	400.620	402.080	403.540	405.000	406.460	408.188	409.916	411.644	413.372	415.100	416.774	418.448	420.122	421.796	423.470	425.146	426.822	428.498	430.174	431.850	434.164	436.478	438.792	441.106	443.420	444.790	446.160	447.530	448.900	450.270	452.478	454.686	456.894	459.102	461.310	462.916	464.522	466.128	467.734	469.340	471.568	473.796	476.024	478.252	480.480]';
gammaInt =  [1.3977	1.3977	1.3976	1.3976	1.3975	1.3975	1.3975	1.3975	1.3974	1.3974	1.3974	1.3973	1.3972	1.3972	1.3971	1.3970	1.3969	1.3968	1.3967	1.3966	1.3965	1.3964	1.3964	1.3963	1.3963	1.3962	1.3961	1.3960	1.3960	1.3959	1.3958	1.3957	1.3957	1.3956	1.3956	1.3955	1.3954	1.3953	1.3951	1.3950	1.3949	1.3948	1.3947	1.3945	1.3944	1.3943	1.3941	1.3940	1.3938	1.3937	1.3935	1.3935	1.3934	1.3934	1.3933	1.3933	1.3932	1.3931	1.3929	1.3928	1.3927	1.3926	1.3925	1.3924	1.3923	1.3922	1.3921	1.3919	1.3918	1.3916	1.3915	1.3914	1.3913	1.3911	1.3910	1.3909	1.3908	1.3907	1.3905	1.3904	1.3903	1.3901	1.3899	1.3896	1.3894	1.3892	1.3891	1.3890	1.3889	1.3888	1.3887	1.3885	1.3883	1.3881	1.3879	1.3877	1.3875	1.3874	1.3872	1.3871	1.3869	1.3867	1.3865	1.3863	1.3861	1.3859]';


gammaCfound = 1;

if ~isstruct(map) % In data is temperature, return gamma as first argument!
  map = interp1(Tint, gammaInt, map,'linear','extrap');    

elseif isfield(map,'gammaC')
  if length(map.gammaC) == 1
    map.gammC = map.gammaC*ones(length(map.NcCorr),1);
  end
  
elseif isfield(map,'T01')
  map.gammaC = interp1(Tint, gammaInt, map.T01,'linear','extrap');
  
elseif isstruct(map)
  T01assumed     = 298;
  gammaCassumed  = interp1(Tint, gammaInt, T01assumed,'linear','extrap');
  map.gammaC = gammaCassumed*ones(size(map.NcCorr));   
else
  gammaCfound = 0;
end
end

function [map, cpCfound] = get_cpC(map)
switch nargin
  case 0
  otherwise
end

cpCfound = 1;


[map, T01found] = get_T01(map);

if isfield(map,'cpC')
  
elseif isfield(map,'gammaC') && isfield(map,'Rc')
  % R = cp - cv, gamma = cp/cv
  % cp = R + cv = R + cp/gamma => cp (1 - 1/gamma) = R
  % cp = R*gamma/(gamma-1)
  map.cpC = map.Rc.*map.gammaC./(map.gammaC-1);  
  
elseif isfield(map,'T01')  
  [map, gammaCfound] = get_gammaC(map);
  
  [map, RcFound]     = get_Rc(map);
  if gammaCfound && RcFound
    map.cpC = map.Rc./(1-1./map.gammaC);
  end
else
  cpCfound = 0;
end


switch nargout
  case 0
    otherwise      
end
end

function [map, cGeomFound] = get_cGeometry(map)
switch nargin
  case 0
  otherwise
    doPlot = 0;
end
cGeomFound = 1;
%% Function declarations
% dC-declaration
f_dC_of_dC1dCh = @(dC1,dCh)( sqrt(1/2 * (dC1^2 + dCh^2)) );
% Function of dC1
f_dCs_of_dC1 = @(dC1)( ( 1.94 +  90.40*dC1 ) * 1e-3 );
% Functions of D_2
f_dC1_of_dC2 = @(D_2)( ( 0.51 + 709.68*D_2 ) * 1e-3 );
f_dCs_of_dC2 = @(D_2)( ( 1.86 +  64.75*D_2 ) * 1e-3 );
% Function of dCs
f_dCh_of_dCs = @(dCs)( 3*dCs );
% Functions of maxN
f_dC1_of_maxN = @(maxN)(  4.16e3 / maxN^0.95 );
f_dC2_of_maxN = @(maxN)( 26.29e3 / maxN^1.08 );
% Function of maxW
f_dC1_of_maxW = @(maxW)( ( 23.22 +  73.99*maxW) * 1e-3 );
f_dC2_of_maxW = @(maxW)( ( 30.65 + 109.11*maxW) * 1e-3 );


[map, NcCorrFound] = get_NcCorr(map);
if ~NcCorrFound, cGeomFound = 0; return, end;
[map, WcCorrfound] = get_WcCorr(map);
if ~WcCorrfound, cGeomFound = 0; return, end;
[map, PiCfound] = get_PiC(map);
if ~PiCfound, cGeomFound = 0; return, end;


%% Modeling efforts...
if isfield(map,'dC1') && isfield(map,'dCh') % Try to estimate dC from available geometry data of map
  map.dC = f_dC_of_dC1dCh(  map.dC1, map.dCh );
  
elseif isfield(map,'dC1')
  if doPlot
    warning('No dCh, using model for dCs(dC1), and dCh=3*dCs');
  end
  map.dCs = f_dCs_of_dC1(   map.dC1 );
  map.dCh = f_dCh_of_dCs(   map.dCs );
  map.dC  = f_dC_of_dC1dCh( map.dC1, map.dCh );
  
elseif isfield(map,'D_2')
  if doPlot    
    warning('No dCh and dC1, using two models: dCs(D_2), dC1(D_2); and finally dCh=3*dCs');
  end
  map.dC1 = f_dC1_of_dC2(   map.D_2 );
  map.dCs = f_dCs_of_dC2(   map.D_2 );
  map.dCh = f_dCh_of_dCs(   map.dCs );
  map.dC  = f_dC_of_dC1dCh( map.dC1, map.dCh );  

else % If no geometry data is available; use estimates of maximum N,W,P
  %[maxN, maxW, maxP] = estimate_maxCompNWP(map);
  
  maxN = max(map.NcCorr);
  maxW = max(map.WcCorr);
  
  
  dC1_N = f_dC1_of_maxN( maxN );
  %dC2_N = f_dC2_of_maxN( maxN );
  dCs_N = f_dCs_of_dC1( dC1_N );
  dCh_N = f_dCh_of_dCs( dCs_N );
  dC_N  = f_dC_of_dC1dCh( dC1_N, dCh_N );
  
  dC1_W = f_dC1_of_maxW( maxW );
  %dC2_W = f_dC2_of_maxW( maxW );
  dCs_W = f_dCs_of_dC1( dC1_W );
  dCh_W = f_dCh_of_dCs( dCs_W ); 
  dC_W  = f_dC_of_dC1dCh( dC1_W, dCh_W );
  
  if doPlot
    [map, Ntcfound] = get_Ntc(map);  
    estimateFromN = [[dC1_N dC2_N dCs_N dCh_N dC_N]*1e3 max(map.Ntc)/1e3 max(map.Ntc)*pi/30*dC2_N/2];
    estimateFromW = [[dC1_W dC2_W dCs_W dCh_W dC_W]*1e3 max(map.Ntc)/1e3 max(map.Ntc)*pi/30*dC2_W/2];
  end
  
  
  [foo, whichMax] = max( [dC_N; dC_W]);
  
  if whichMax == 1
    map.dC1 = dC1_N;    %map.D_2 = dC2_N;    
    map.dCs = dCs_N;    map.dCh = dCh_N;    map.dC  = dC_N;
  else
    map.dC1 = dC1_W;    %map.D_2 = dC2_W;    
    map.dCs = dCs_W;    map.dCh = dCh_W;    map.dC  = dC_W;
  end
  
%  dCfound = 0;
end


%% Plot things?
if doPlot
  subplot(2,1,1);
  plot_compMapFlow(map);  
end


%% Output handling
switch nargout
  case 0
  otherwise
end
end

function [map] = get_adiabatic_efficiency(map,k_fric,k_q)
switch nargin
  case 0
  otherwise
    doPlot = 0;
end

[map, PhiFound]  	= get_PhiCm(map);
[map, etaCfound]  	= get_etaC(map);
[map, PiCfound]     = get_PiC(map);
[map, T01found]     = get_T01(map);
[map, gammaCfound]  = get_gammaC(map);
if PhiFound && PiCfound && T01found && gammaCfound
    map.lambda_heat     = k_q./(map.PhiC2.*map.MachC.^3);
    map.lambda_ad       = map.lambda-map.lambda_heat;
    map.lambda_euler_ad = map.lambda_ad./(1+k_fric./map.PhiC);
    map.hoC_ad          = map.lambda_ad.*map.Uc.^2;
    map.etaC_ad         = map.hosC./map.hoC_ad;
else  

end

switch nargout
  case 0
  otherwise
end
end

function [map, WcCorrfound] = get_WcCorr(map)
switch nargin
  case 0   
    otherwise
end
f_WcCorr = @(Wc, p01, T01, pCref, TCref)( Wc .* sqrt( T01 ./ TCref ) ./ ( p01 ./ pCref ) );

WcCorrfound = 1;
pCrefAss = 100e3; % Default ASSumed compressor map reference pressure
TCrefAss = 298;   % Default ASSumed compressor map reference temperature
p01Ass   = 100e3; % Default ASSumed compressor inlet pressure
T01Ass   = 298;   % Default ASSumed compressor inlet temperature
RairAss  = 287;   %
%% <--------------- Korrigerat massflöde
if isfield(map, 'WcCorr') 
  
%% <--------------- Massflöde  
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'T01') && isfield(map, 'pCref') && isfield(map, 'TCref') 
  map.WcCorr = f_WcCorr(map.Wc, map.p01, map.T01, map.pCref, map.TCref);
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'T01') && isfield(map, 'pCref')
  map.WcCorr = f_WcCorr(map.Wc, map.p01, map.T01, map.pCref, TCrefAss);
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'T01') && isfield(map, 'TCref')  
  map.WcCorr = f_WcCorr(map.Wc, map.p01, map.T01, pCrefAss, map.TCref);
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'T01') 
  map.WcCorr = f_WcCorr(map.Wc, map.p01, map.T01, pCrefAss, TCrefAss);    
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'pCref') && isfield(map, 'TCref') 
  map.WcCorr = f_WcCorr(map.Wc, map.p01, T01Ass, map.pCref, map.TCref);    
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'pCref') 
  map.WcCorr = f_WcCorr(map.Wc, map.p01, T01Ass, map.pCref, TCrefAss);
elseif isfield(map, 'Wc') && isfield(map, 'p01') && isfield(map, 'TCref') 
  map.WcCorr = f_WcCorr(map.Wc, map.p01, T01Ass, pCrefAss, map.TCref);  
elseif isfield(map, 'Wc') && isfield(map, 'T01') && isfield(map, 'pCref') 
  map.WcCorr = f_WcCorr(map.Wc, p01Ass, map.T01, map.pCref, TCrefAss);
elseif isfield(map, 'Wc') && isfield(map, 'T01') && isfield(map, 'TCref') 
  map.WcCorr = f_WcCorr(map.Wc, p01Ass, map.T01, pCrefAss, map.TCref);    
elseif isfield(map, 'Wc') && isfield(map, 'pCref') && isfield(map, 'TCref') 
  map.WcCorr = f_WcCorr(map.Wc, p01Ass, T01Ass, map.pCref, map.TCref);  

%% <--------------- Korrigerat volymsflöde  
elseif isfield(map,'VcCorr') 
  if isfield(map, 'p01')     
    if isfield(map,'T01')
      if isfield(map,'pCref') 
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./map.T01).*map.p01/RairAss./map.T01, map.p01, map.T01, map.pCref, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./map.T01).*map.p01/RairAss./map.T01, map.p01, map.T01, map.pCref, TCrefAss);   end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./map.T01).*map.p01/RairAss./map.T01, map.p01, map.T01, pCrefAss, map.TCref);  
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./map.T01).*map.p01/RairAss./map.T01, map.p01, map.T01, pCrefAss, TCrefAss); end
      end
    else
      if isfield(map,'pCref')
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr*sqrt(map.TCref./T01Ass).*map.p01/RairAss./T01Ass, map.p01, T01Ass, map.pCref, map.TCref);  
        else                     map.WcCorr = f_WcCorr(map.VcCorr*sqrt(TCrefAss./T01Ass).*map.p01/RairAss./T01Ass, map.p01, T01Ass, map.pCref, TCrefAss);     end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr*sqrt(map.TCref./T01Ass).*map.p01/RairAss./T01Ass, map.p01, T01Ass, pCrefAss, map.TCref);   
        else                     map.WcCorr = f_WcCorr(map.VcCorr*sqrt(TCrefAss./T01Ass).*map.p01/RairAss./T01Ass, map.p01, T01Ass, pCrefAss, TCrefAss);     end
      end
    end
  else
    if isfield(map,'T01')
      if isfield(map,'pCref') 
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./map.T01).*p01Ass/RairAss./map.T01, p01Ass, map.T01, map.pCref, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./map.T01).*p01Ass/RairAss./map.T01, p01Ass, map.T01, map.pCref, TCrefAss);  end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./map.T01).*p01Ass/RairAss./map.T01, p01Ass, map.T01, pCrefAss, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./map.T01).*p01Ass/RairAss./map.T01, p01Ass, map.T01, pCrefAss, TCrefAss);  end
      end
    else
      if isfield(map,'pCref')
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./T01Ass).*p01Ass/RairAss./T01Ass, p01Ass, T01Ass, map.pCref, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./T01Ass).*p01Ass/RairAss./T01Ass, p01Ass, T01Ass, map.pCref, TCrefAss);  end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(map.TCref./T01Ass).*p01Ass/RairAss./T01Ass, p01Ass, T01Ass, pCrefAss, map.TCref);
        else                     map.WcCorr = f_WcCorr(map.VcCorr.*sqrt(TCrefAss./T01Ass).*p01Ass/RairAss./T01Ass, p01Ass, T01Ass, pCrefAss, TCrefAss); end
      end
    end    
  end  
  
%<--------------- Volymsflöde  
elseif isfield(map,'Vc') 
  if isfield(map, 'p01')     
    if isfield(map,'T01')
      if isfield(map,'pCref') 
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./map.T01, map.p01, map.T01, map.pCref, map.TCref);
        else                     map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./map.T01, map.p01, map.T01, map.pCref, TCrefAss);  end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./map.T01, map.p01, map.T01, pCrefAss, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./map.T01, map.p01, map.T01, pCrefAss, TCrefAss); end
      end
    else
      if isfield(map,'pCref')
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./T01Ass, map.p01, T01Ass, map.pCref, map.TCref);  
        else                     map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./T01Ass, map.p01, T01Ass, map.pCref, TCrefAss);   end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./T01Ass, map.p01, T01Ass, pCrefAss, map.TCref);   
        else                     map.WcCorr = f_WcCorr(map.Vc.*map.p01./RairAss./T01Ass, map.p01, T01Ass, pCrefAss, TCrefAss);    end
      end
    end
  else
    if isfield(map,'T01')
      if isfield(map,'pCref') 
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./map.T01, p01Ass, map.T01, map.pCref, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./map.T01, p01Ass, map.T01, map.pCref, TCrefAss); end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./map.T01, p01Ass, map.T01, pCrefAss, map.TCref);
        else                     map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./map.T01, p01Ass, map.T01, pCrefAss, TCrefAss);  end
      end
    else
      if isfield(map,'pCref')
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./T01Ass, p01Ass, T01Ass, map.pCref, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./T01Ass, p01Ass, T01Ass, map.pCref, TCrefAss);  end
      else
        if isfield(map,'TCref'), map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./T01Ass, p01Ass, T01Ass, pCrefAss, map.TCref); 
        else                     map.WcCorr = f_WcCorr(map.Vc*p01Ass./RairAss./T01Ass, p01Ass, T01Ass, pCrefAss, TCrefAss); end
      end
    end    
  end   
  
else  
  WcCorrfound = 0;    
end


 switch nargout
   case 0     
   otherwise     
 end
end

function [map, etaCfound] = get_etaC(map)

switch nargin
  case 0
  otherwise
end

etaCfound = 1;

if isfield(map,'etaC')
  map.etaC = min(1, max(0, map.etaC));
else
  etaCfound = 0;  
end
end

function [map, Wcfound] = get_Wc(map)

switch nargin
  case 0
  otherwise
end

f_Wc = @(WcCorr, p01, T01, pCref, TCref)(WcCorr.*(p01./pCref)./sqrt(T01./TCref) );

Wcfound = 1;
p01Ass = 100e3;
T01Ass = 298;
TCrefAss = 298;
pCrefAss = 100e3;


if isfield(map,'Wc')
  
else
  [map, WcCorrfound] = get_WcCorr(map);
  if WcCorrfound
    if isfield(map,'p01')
      if isfield(map,'T01')
        if isfield(map,'TCref')
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, map.p01, map.T01, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, map.p01, map.T01, pCrefAss, map.TCref);
          end
        else
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, map.p01, map.T01, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, map.p01, map.T01, pCrefAss, TCrefAss);
          end
        end
      else
        if isfield(map,'TCref')
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, map.p01, T01Ass, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, map.p01, T01Ass, pCrefAss, map.TCref);
          end
        else
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, map.p01, T01Ass, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, map.p01, T01Ass, pCrefAss, map.TCref);
          end
        end
      end
    else
      if isfield(map,'T01')
        if isfield(map,'TCref')
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, p01Ass, map.T01, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, p01Ass, map.T01, pCrefAss, map.TCref);
          end
        else
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, p01Ass, map.T01, map.pCref, TCrefAss);
          else
            map.Wc = f_Wc(map.WcCorr, p01Ass, map.T01, pCrefAss, TCrefAss);
          end
        end
      else
        if isfield(map,'TCref')
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, p01Ass, T01Ass, map.pCref, map.TCref);
          else
            map.Wc = f_Wc(map.WcCorr, p01Ass, T01Ass, pCrefAss, map.TCref);
          end
        else
          if isfield(map,'pCref')
            map.Wc = f_Wc(map.WcCorr, p01Ass, T01Ass, map.pCref, TCrefAss);
          else
            map.Wc = f_Wc(map.WcCorr, p01Ass, T01Ass, pCrefAss, TCrefAss);
          end
        end
      end
    end
  else
    Wcfound = 0;
  end
end


switch nargout
  case 0
  otherwise
end
end

function vec = get_vector(M)
[~,col] = size(M);
vec = zeros(col,1);
for jj = 1:col
 	extract = M(:,jj);
    vec(jj) = mean(extract(isnan(extract)<1));
end
end



