% Demo codes for the Copula-based Granger causality for mixed data (e.g., LFP and Spike)
% 
% This demo is for Simulation 1 by doing the following:
% 1. Apply the copula-based model to fit multiple-trial mixed data and
% 2. Estimate Granger causality between mixed time series.
% 3. Plot parameter estimation of model
% 4. Plot estimated Granger causality
% 5. Permutation for significance test (warning: it takes a while to run permutation)
%
% Meng Hu, Mingyao Li, Wu Li and Hualou Liang, Joint Analysis of Spikes 
%   and Local Field Potentials using Copula, NeuroImage, 133: 457 ? 467, 2016
%
% Meng Hu @ Liang's lab at Drexel University, 2015
%

%%
clear

M=100; %% trial
N=1000; %% data length
% rho=0.3; %% copula rho
porder=2; %% model order
% qq=-1.5; %% parameter to adjust spike rate
% parameter for optimization 
options = optimset('GradObj','on','Display','notify','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);

% load sample mixed data
load Simulation1MixedData  % [trial x time x channel], The 1st channel is LFP, whereas the 2nd channel is Spike

gc12_all = [];
gc21_all = [];
para_all = [];
for n=1:M
    
    Y1=squeeze(dat(n,:,1)); % LFP
    Y2=squeeze(dat(n,:,2)); % Spike

% Fit data and calculate Granger causality
    [gc12 gc21 para Lik]=Mixed_GC_Gauss_fminunc(Y1,Y2,porder,options);
%     msflag = 0; % 1 - multistart to mitigate the possible local minima, eg. 10); 0(no multistart)
%     [gc12 gc21 para Lik]=Mixed_GC_Gauss_fminunc_multistart(Y1,Y2,porder,options,msflag);
    gc12_all = [gc12_all,gc12]; % GC from 1->2
    gc21_all = [gc21_all,gc21]; % GC from 2->1
    para_all = [para_all;para]; % model parameter

n       
end


%% plot parameter estimation of model

para_all(:,end) = -1+2*exp(-exp(para_all(:,end)));
tv = [0.2 -0.5 0.3 0 1 -1.5 0 0 0.6 0 0.3]; %% true value
evm = mean(para_all,1);%% estimated value mean
ever = std(para_all,[],1)/sqrt(M);%% estimated value se
nerr=2; %% times of se

figure;
b=bar([tv;evm]');
set(b(1),'FaceColor',[1 1 1]*0.85);
set(b(2),'FaceColor','black');
legend('True','Estimated')
hold on
for ii=1:11
   plot([ii,ii]+0.15,[evm(ii)-nerr*ever(ii),evm(ii)+nerr*ever(ii)],'-k','LineWidth',2)
end
set(gca,'FontSize',16,'fontWeight','bold')
xtl={{'\beta_1'} {'\beta_2'} {'\beta_3'} {'\beta_4'} {'\beta_5'} {'\beta_6'} {'\beta_7'} {'\beta_8'} {'\beta_9'} {'\beta_1_0'} {'r'}};
my_xticklabels(gca,1:11,xtl,'FontSize',16);
ax=axis;
axis([ax(1) ax(2) -1.6 1.1])
ylabel('Parameter estimation')


%% plot estimated Granger causality

gc_all =[gc12_all; gc21_all]';
gcm = mean(gc_all,1);%% estimated value mean
gcer = std(gc_all,[],1)/sqrt(M);%% estimated value se
nerr=2; %% times of se

figure;
b=bar([gcm]', 0.5);
set(b(1),'FaceColor',[1 1 1]*0.85);
hold on
for ii=1:2
   plot([ii,ii],[gcm(ii)-nerr*gcer(ii),gcm(ii)+nerr*gcer(ii)],'-k','LineWidth',2)
end
set(gca,'FontSize',16,'fontWeight','bold')
xtl={{'LFP->Spike'} {'Spike->LFP'}};
my_xticklabels(gca,1:2,xtl,'FontSize',16);
ax=axis;
axis([ax(1) ax(2) 0 100])
ylabel('Granger causality')
title('True direction: Spike->LFP')


%% Permutation for significance test

gc12_all_p = [];
gc21_all_p = [];
nn=0;
for k=1:10 % permutations

    dat_p = dat;
    dat_p(:,:,1) = dat_p(randperm(M),:,1);
    dat_p(:,:,2) = dat_p(randperm(M),:,2);

    for n=1:100

        Y1_p=squeeze(dat_p(n,:,1)); % LFP
        Y2_p=squeeze(dat_p(n,:,2)); % Spike

% Fit data and calculate Granger causality
        [gc12_p gc21_p para_p Lik_p]=Mixed_GC_Gauss_fminunc(Y1_p,Y2_p,porder,options);
        gc12_all_p = [gc12_all_p,gc12_p]; % GC from 1->2
        gc21_all_p = [gc21_all_p,gc21_p]; % GC from 2->1

        nn=nn+1     
    end

end

%% plot estimated Granger causality + significance test by permutation

% calculate threshold for significance test (p=0.01)
thresh12 = sort(gc12_all_p);
thresh12 = thresh12(fix(length(thresh12)*0.99));
thresh21 = sort(gc21_all_p);
thresh21 = thresh21(fix(length(thresh21)*0.99));


gc_all =[gc12_all; gc21_all]';
gcm = mean(gc_all,1);%% estimated value mean
gcer = std(gc_all,[],1)/sqrt(M);%% estimated value se
nerr=2; %% times of se

figure;
b=bar([gcm]', 0.5);
set(b(1),'FaceColor',[1 1 1]*0.85);
hold on
for ii=1:2
   plot([ii,ii],[gcm(ii)-nerr*gcer(ii),gcm(ii)+nerr*gcer(ii)],'-k','LineWidth',2)
end
set(gca,'FontSize',16,'fontWeight','bold')
xtl={{'LFP->Spike'} {'Spike->LFP'}};
my_xticklabels(gca,1:2,xtl,'FontSize',16);
ax=axis;
axis([ax(1) ax(2) 0 100])
ylabel('Granger causality')
title('True direction: Spike->LFP')
hold on
plot([1,2],[thresh12,thresh21],'x','color','red')