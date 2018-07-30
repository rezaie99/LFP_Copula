% Demo codes for the Copula-based Granger causality for mixed data (e.g., LFP and Spike)
%
% Meng Hu, Mingyao Li, Wu Li and Hualou Liang, Joint Analysis of Spikes 
%   and Local Field Potentials using Copula, NeuroImage, 133: 457 ? 467, 2016
%
% Meng Hu @ Liang's lab at Drexel University, 2015
%

function [gc12 gc21 para Lik]=Mixed_GC_Gauss_fminunc(Y1,Y2,porder,options)

% Establish copula-based model to fit data and estimate Granger
% causality between mixed time series signals
%
% Input:
%   Y1 - cont
%   Y2 - binary
%   porder - model order
%   options - parameter for optimization

% Output:
%   gc12 - GC from 1 to 2
%   gc21 - GC from 2 to 1
%   para - parameter of full model
%   Lik - Likelihood of full model

% Additional:
%   Lik1 - Likelihood of reduced model eliminating Y2->Y1 (binary->cont) causal direction 
%   Lik2 - Likelihood of reduced model eliminating Y1->Y2 (cont->binary) causal direction

%% whole modeling
gamma_init=corrcoef(Y1,Y2);
gamma_init=log(log(2/(gamma_init(1,2)+1)));
para_init=([rand(1,4*porder+2) gamma_init]); 
[x,fval]=fminunc(@(para) mixedgc_obj_Gauss(Y1,Y2,para,porder),para_init,options);
para=x;
Lik=-fval;  

%% eliminate a causal direction (Y2->Y1) (binary->cont)
N=length(Y1);
para_init1=para_init([1:1+porder 2+2*porder:end]);
[x,fval]=fminunc(@(para) mixedgc_obj_Gauss2(Y1,Y2,para,porder),para_init1,options);
para1=x;
Lik1=-fval;      

%% eliminate a causal direction (Y1->Y2) (cont->binary)
N=length(Y1);
para_init2=para_init([1:2*porder+2 3*porder+3:end]);
[x,fval]=fminunc(@(para) mixedgc_obj_Gauss1(Y1,Y2,para,porder),para_init2,options);
para2=x;
Lik2=-fval;      
  

%% ratio of likelihood

% Y2->Y1
gc21=real(Lik)-real(Lik1); % Spike -> LFP

% Y1->Y2
gc12=real(Lik)-real(Lik2); % LFP -> Spike


end