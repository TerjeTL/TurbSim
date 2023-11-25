function [best_par, Vn_best]=lsoptim(func,parameters,maxSteps,Lin_A,Lin_b,varargin);

% [best_par]=lsoptim(funcName,parameters,maxSteps,A,b,varargin)
%
% Solves an uncostrained non-linear least squares problem using a
% Levenberg-Marquardt like method.
%
% Linear and non-negative constraints is supported through Ax>=b. This is solved using
% log barrier internal point method. Therefore the equality is only approached
% assymptotically.  The intended usage for this feature is to damp the unconstrained
% method so that it doesn't jump into a region with unfeasible parameter values
% during the search.
%
%  funcName     - name of the function to be optimized
%  parameters   - default values for the parameters
%  maxSteps     - maximum number of iterations during the search (optional, default=5)
%  A            - Matrix A for linear inequality constraints Ax>b (optional).
%  b            - Vector b for linear (optional).
%  varargin     - Support for additional argument to be sent to the
%                 criteria function (optional).

% Copyright (C) 1997, 2001, 2003, 2005, 2016, 2017 Lars Eriksson
% Version 1.4
%  + Added check for imaginary numbers.
%  + Added support for linear inequalities with log barrier method.
% Version 1.3
%  + Added support for linear inequalities, simple truncation of step length.
% Version 1.2
%  + Added support for non negative variables.
%  + Added support for varargin and function handles, contributed by Peter Spring
% Version 1.1ima
%  + Changed the name of the function name 'function' to func
%  + Added support for transposed input parameter vector.
%  + Added support for transposed result vector.

global PARAMETER_TRACE
PARAMETER_TRACE=[];

DispFlag = 1 ; % Flag for displaying messages during the optimization

WhosOutput = whos('func') ;
FuncIsHandle = 0;
if strcmp(WhosOutput.class,'function_handle')
  %disp('"func" is a function handle')
  FuncIsHandle = 1 ;
end


% -- provide variable inputs for 'evalstring'
VarParString = [] ;
if nargin > 5
  for i = 6:nargin
    VarParString = [VarParString,',varargin{',num2str(i-3),'}'] ;
  end
end

% -- parameters for the optimization method
if nargin<3||isempty(maxSteps),
  maxSteps=5;
end

%delta  = 1e-3;
delta  = 1e-2;

lambda = 10;

%convergenceTol=1e-6;
convergenceTol=1e-4;


% -- Initialize the data
if FuncIsHandle
  evalstring=['feval(func,temp_par',VarParString,')'] ;
else
  evalstring=[func,'(temp_par',VarParString,')'] ;
end

[m,n]=size(parameters);
if min(m,n)>1,
  error('The parameter vector must be an Nx1 column vector.');
end
if m<n,
  parameters=parameters';
  warning('The input parameter vector should be an Nx1 column vector!');
end
PARAMETER_TRACE=[PARAMETER_TRACE; parameters'];
temp_par=parameters;
F=eval(evalstring);
[m,n]=size(F);
if n>1,
  error('The function ''func'' must return an Nx1 vector.');
end
Vn=F'*F;
best_par=parameters;
mu=2;%sqrt(Vn_best);

% -- initialization of optimization variables and parameters
nvar=length(parameters);
J   = zeros(length(F),nvar); % The Jacobian for f
I   = eye(nvar);       % The Identity matrix for Levenberg-Marquardt method
dk  = ones(nvar,1);
noSteps=0;

% -- Variables for handling linear inequality constraints.
Constrained=0;
if nargin>3,
  % -- Check that both A and b are provided
  if nargin==4,
    error('Only A is given, b must also be given to ''lsoptim''');
  else
    % -- If both A and b are empty nothing is made.
    if isempty(Lin_A)||isempty(Lin_b),
      % -- Test if only one is non-empty.
      if (~isempty(Lin_A))||(~isempty(Lin_b)),
        error('Only one of A or b is non empty');
      end
    else
      % -- Check that the initial guess fulfills the inequality constraint.
      if Lin_A*parameters<=Lin_b,
        error(['The initial condition does not fulfill the ineqality' ...
          ' constraints Ax>b']);
      else
        % -- We are ready to rock'n roll with constraints.
        Constrained=1;
        disp('Problem set up with constraints')
      end
    end
  end
end

ParameterSize=abs(parameters);
ZeroIndex=find(ParameterSize==0);
if ~isempty(ZeroIndex)
  ParameterSize(ZeroIndex)=ones(size(ZeroIndex));
  if DispFlag,
    fprintf(1,'The initial values to the following parameters were zero: ');
    disp(ZeroIndex)
    fprintf(1,'Nominal values used in the numeric differentiation are set to 1.')
  end
end

Vn_old=inf;
Vn_best=Vn;
Vn_prev=Vn;

if DispFlag,
  disp('Optimization procedure started.');
  fprintf(1,'Initial Error: %f,   ',Vn_best);
end
%HRgui('plot');

while ((noSteps<maxSteps) & (abs((Vn_old-Vn)/Vn)>convergenceTol)) %&(dk'*dk>delta*parameters'*parameters*delta)),
    % ------ Calculate the gradient for all optimization parameters
    for i=1:nvar,
        temp_par=parameters;
        deltaX=delta*abs(ParameterSize(i));
        %if (deltaX<1e-14),
        %  deltaX=1e-14;
        %  [delta i]
        %end
        temp_par(i)=parameters(i)+deltaX;
        J(:,i)=(eval(evalstring)-F)/deltaX;
        if DispFlag,
            fprintf(1,'\b\b\b\b%3d%%',round(i/nvar*100));
        end
    end
    lambda=min(Vn_best/5000,lambda);
    if Constrained,
        % ---- Run the log barrier penalty function
        H=J'*J+lambda*I+mu*Lin_A'*diag(1./(Lin_A*parameters-Lin_b).^2)*Lin_A;
        JJ=J'*F-mu*Lin_A'*(1./(Lin_A*parameters-Lin_b));
        dk=-H\JJ;
    else
        dk=-(J'*J + lambda*I)\(J'*F);
    end
    
    % ------ Delta for the calculation of the gradient approximation,
    % ------ which is decreased when we are close to the optimum
    delta=min(max(abs(dk./ParameterSize))/100,delta);
    
    % ------ make sure that constraints are fulfilled within the steplength
    t=1;
    if Constrained,
        % ---- Yes we have constraints
        % ---- Calculate the steplength, t, to the constraints to avid stepping
        % ---- into the constraints.
        t=getMaxStep(parameters,dk,Lin_A,Lin_b,DispFlag);
        if (t==1),
            % -- Constraints not hit decreasing log barrier, weight.
            mu=mu/10;
        else
            % -- Constraints hit increasing log barrier, weight.
            mu=mu/2;
        end
    end
    
    Vn_old=Vn;
    alpha=1;
    F_old=F;
    % ------ Try to ensure that it is a descent.
    Vn_stored=inf;
    while (alpha>1e-9),
        % ---- Take one step in the search direction.
        temp_par=parameters+t*dk;
        
        % ---- Evaluate the criterion
        F=eval(evalstring);
        if ~isreal(F),
            disp(' ')
            disp('The parameters that give imaginary numbers.')
            disp([PARAMETER_TRACE(end,:)'])
            error('Imaginary numbers encountered, consider adding constraints to the problem.');
        end
        
        % ---- Calculate the loss function and test if the current parameter
        % ---- values improves the criterion.
        Vn=F'*F;
        if (Vn>Vn_prev),
            % -- It is not yet a descent.
            % -- Try to pull the step towards the descent direction.
            if (Vn>Vn_stored),
                % Decent attempt did not improve the fit, revert to original Newton step and
                t=t_stored;
                dk=dk_stored;
                Vn=Vn_stored;
                alpha=0;
                if DispFlag,
                    fprintf(1,' descent attempt failed.');
                end
            else
                Vn_stored=Vn; % store to check improvement.
                t_stored=t;
                dk_stored=dk;
                alpha=alpha/8;
                if DispFlag,
                    fprintf(1,' %f',alpha);
                end
                % -- Calculate a new search direction
                if Constrained,
                    H=J'*J+lambda/alpha*I+mu*Lin_A'*diag(1./(Lin_A*parameters-Lin_b).^2)*Lin_A;
                    JJ=J'*F-mu*Lin_A'*(1./(Lin_A*parameters-Lin_b));
                    dk=-H\JJ;
                else
                    dk=-(J'*J + lambda/alpha*I)\(J'*F_old);
                end
                t=1;
                if Constrained,
                    % --- Check the constraints
                    t=getMaxStep(parameters,dk,Lin_A,Lin_b,DispFlag);
                end
            end
        else
            alpha=0;
        end
    end
    Vn_Prev=Vn;
    if DispFlag,
        if (Vn>Vn_best),
            fprintf(1,' no descent');
        end
    end
    parameters=temp_par;
    PARAMETER_TRACE=[PARAMETER_TRACE; parameters'];
    noSteps=noSteps+1;
    if DispFlag,
        fprintf('\nIteration # %3d,  Error: %f,    ',noSteps,Vn);
        %drawnow('update')
    end
    if isnan(Vn),
        disp([])
        disp('The parameter trace')
        disp([PARAMETER_TRACE])
        error('NaN encountered during the optimization, optimization exited.')
    end
    if Vn<Vn_best,
        %	HRgui('plot');
        Vn_best=Vn;
        %	best_par=[best_par parameters];
        best_par=parameters;
    end
end
if DispFlag,
    if (abs((Vn_old-Vn)/Vn)<=convergenceTol),
        fprintf(1,'\nOptimization converged.\n');
    else
        fprintf(1,'\nOptimization exited, maximum number of iterations reached.\n');
    end
end



function t=getMaxStep(parameters,dk,Lin_A,Lin_b,DispFlag);

% ----- Remove the components that have dk(i)==0
index=1:length(dk);

% ----- What is the longest possible steplength?
% ----- Solve all rows in A(x_o+t*d_k)=b for t.
t_vec=(Lin_b-Lin_A(:,index)*parameters(index))./(Lin_A(:,index)*dk(index));

% ----- look for the smallest positive value
t_limit = min(t_vec(find(t_vec>0)));

if isempty(t_limit) | (t_limit>1),
    t=1;
else
    t=t_limit*.99;
    if DispFlag,
        fprintf(1,[' Constraint hit t=' num2str(t)])
    end
end
