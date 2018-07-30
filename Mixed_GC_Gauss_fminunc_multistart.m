% Demo codes for the Copula-based Granger causality for mixed data (e.g., LFP and Spike)
%  with multistart option
%
% Meng Hu, Mingyao Li, Wu Li and Hualou Liang, Joint Analysis of Spikes 
%   and Local Field Potentials using Copula, NeuroImage, 133: 457 ? 467, 2016
%
% Meng Hu @ Liang's lab at Drexel University
%

function [gc12 gc21 para Lik]=Mixed_GC_Gauss_fminunc_multistart(Y1,Y2,porder,options, msflag)

% Establish copula-based model to fit data and estimate Granger
% causality between mixed time series signals
%
% Input:
%   Y1 - cont
%   Y2 - binary
%   porder - model order
%   options - parameter for optimization
%   msfalg - multi-start flag, 0 (false, default), 1 (true) of 10 starts

% Output:
%   gc12 - GC from 1 to 2
%   gc21 - GC from 2 to 1
%   para - parameter of full model
%   Lik - Likelihood of full model

% Additional:
%   Lik1 - Likelihood of reduced model eliminating Y2->Y1 (binary->cont) causal direction 
%   Lik2 - Likelihood of reduced model eliminating Y1->Y2 (cont->binary) causal direction

%%
if nargin < 5,
  msflag = 0;
end
runs=10; % number of random starts

%% whole modeling
gamma_init=corrcoef(Y1,Y2);
gamma_init=log(log(2/(gamma_init(1,2)+1)));
para_init=[rand(1,4*porder+2) gamma_init]; 
if msflag % multistart flag
    % Set options for a optimization solver
    options = optimset('GradObj','on','Display','notify','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);

    % 1. Set up the PROBLEM structure
    problem = createOptimProblem('fminunc','objective',@(para) mixedgc_obj_Gauss(Y1,Y2,para,porder),'x0',para_init,'options',options);

    %2. Construct the MultiStart solver
    ms = MultiStart;

    % 3. Run the solver from 'runs' random start points
    csps = CustomStartPointSet([rand(runs,10) gamma_init*ones(runs,1)]);
    [x,fval] = run(ms,problem,csps);
else
    [x,fval]=fminunc(@(para) mixedgc_obj_Gauss(Y1,Y2,para,porder),para_init,options);
end
para=x;
Lik=-fval;  


%% eliminate a causal direction (Y2->Y1) (binary->cont)
N=length(Y1);
para_init1=para_init([1:1+porder 2+2*porder:end]);
if msflag % multistart flag
    % Set options for a optimization solver
    options = optimset('GradObj','on','Display','notify','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);

    % 1. Set up the PROBLEM structure
    problem = createOptimProblem('fminunc','objective',@(para) mixedgc_obj_Gauss2(Y1,Y2,para,porder),'x0',para_init1,'options',options);

    %2. Construct the MultiStart solver
    ms = MultiStart;

    % 3. Run the solver from 5 random start points
    csps = CustomStartPointSet([rand(runs,8) gamma_init*ones(runs,1)]);
    [x,fval] = run(ms,problem,csps);
else
    [x,fval]=fminunc(@(para) mixedgc_obj_Gauss2(Y1,Y2,para,porder),para_init1,options);
end
para1=x;
Lik1=-fval;      

%% eliminate a causal direction (Y1->Y2) (cont->binary)
N=length(Y1);
para_init2=para_init([1:2*porder+2 3*porder+3:end]);
if msflag % multistart flag
    % Set options for a optimization solver
    options = optimset('GradObj','on','Display','notify','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);

    % 1. Set up the PROBLEM structure
    problem = createOptimProblem('fminunc','objective',@(para) mixedgc_obj_Gauss1(Y1,Y2,para,porder),'x0',para_init2,'options',options);

    %2. Construct the MultiStart solver
    ms = MultiStart;

    % 3. Run the solver from 5 random start points
    csps = CustomStartPointSet([rand(runs,8) gamma_init*ones(runs,1)]);
    [x,fval] = run(ms,problem,csps);
else
    [x,fval]=fminunc(@(para) mixedgc_obj_Gauss1(Y1,Y2,para,porder),para_init2,options);
end
para2=x;
Lik2=-fval;      
  

%% ratio of likelihood

% Y2->Y1
gc21=real(Lik)-real(Lik1); % Spike -> LFP

% Y1->Y2
gc12=real(Lik)-real(Lik2); % LFP -> Spike


end