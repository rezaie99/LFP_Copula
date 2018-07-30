function [gcm, gcer]= demo(LFP, Spikes,index)


M=size(LFP,2); %% trial
N=size(LFP,1); %% data length
% rho=0.3; %% copula rho
porder=2; %% model order
% qq=-1.5; %% parameter to adjust spike rate
% parameter for optimization 
options = optimset('GradObj','on','Display','notify','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);

% load sample mixed data
%load Simulation1MixedData  % [trial x time x channel], The 1st channel is LFP, whereas the 2nd channel is Spike

gc12_all = [];
gc21_all = [];
para_all = [];
for n=1:M
     fprintf('Iteration %d out of %d\n',n, M);
%    Y1=squeeze(dat(n,:,1)); % LFP
%    Y2=squeeze(dat(n,:,2)); % Spike
    Y1 = LFP(:,n);
    Y2 = Spikes(:,n);
    
    Y1 = (Y1-mean(Y1))/std(Y1);
    windowWidth = 20; % Whatever you want.
    kernel = ones(windowWidth,1) ;
    Y1 = Y1 - filter(kernel, 20, Y1);
    
    windowWidth = 50; % Whatever you want.
    kernel = ones(windowWidth,1) ;
    out = filter(kernel, 50, abs(Y1));
    Y1=Y1./out;
    Y1 = Y1(100:end);
    Y2 = Y2(100:end);

    %plot(Y1)
    %hold on
    %plot(out)
    %plot(Y2)

% Fit data and calculate Granger causality
    try
        [gc12, gc21, para, Lik]=Mixed_GC_Gauss_fminunc(Y1,Y2,porder,options);
    %     msflag = 0; % 1 - multistart to mitigate the possible local minima, eg. 10); 0(no multistart)
    %     [gc12 gc21 para Lik]=Mixed_GC_Gauss_fminunc_multistart(Y1,Y2,porder,options,msflag);
        gc12_all = [gc12_all,gc12]; % GC from 1->2
        gc21_all = [gc21_all,gc21]; % GC from 2->1
        para_all = [para_all;para]; % model parameter  
    catch
        continue
    end
    
end


%% plot parameter estimation of model

para_all(:,end) = -1+2*exp(-exp(para_all(:,end)));
tv = [0.2 -0.5 0.3 0 1 -1.5 0 0 0.6 0 0.3]; %% true value
evm = mean(para_all,1);%% estimated value mean
ever = std(para_all,[],1)/sqrt(M);%% estimated value se
nerr=2; %% times of se

figure;
b=bar([evm;evm]');
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
title(sprintf('parameters %d',index))

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
axis([ax(1) ax(2) 0 10])
ylabel('Granger causality')
title(sprintf('True direction: @ %d', index))


%% Permutation for significance test

% gc12_all_p = [];
% gc21_all_p = [];
% nn=0;
% for k=1:1 % permutations
% 
%     dat_p1= LFP(:,randperm(M));
%     dat_p2= Spikes(:,randperm(M));
% 
%     for n=1:M
% 
%         Y1_p=dat_p1(:,n); % LFP
%         Y2_p=dat_p2(:,n); % Spike
%             
%         Y1_p = (Y1_p-mean(Y1_p))/std(Y1_p);
%         windowWidth = 20; % Whatever you want.
%         kernel = ones(windowWidth,1) ;
%         Y1_p = Y1_p - filter(kernel, 20, Y1_p);
%     
%         windowWidth = 50; % Whatever you want.
%         kernel = ones(windowWidth,1) ;
%         out = filter(kernel, 50, abs(Y1_p));
%         Y1_p=Y1_p./out;
%         Y1_p = Y1_p(100:end);
%         Y2_p = Y2_p(100:end);
% 
% 
% 
% % Fit data and calculate Granger causality
%         [gc12_p, gc21_p, para_p, Lik_p]=Mixed_GC_Gauss_fminunc(Y1_p,Y2_p,porder,options);
%         gc12_all_p = [gc12_all_p,gc12_p]; % GC from 1->2
%         gc21_all_p = [gc21_all_p,gc21_p]; % GC from 2->1
% 
%         nn=nn+1;
%         disp(nn)
%     end
% 
% end
% 
% %% plot estimated Granger causality + significance test by permutation
% 
% % calculate threshold for significance test (p=0.01)
% thresh12 = sort(gc12_all_p);
% thresh12 = thresh12(fix(length(thresh12)*0.99));
% thresh21 = sort(gc21_all_p);
% thresh21 = thresh21(fix(length(thresh21)*0.99));
% 
% 
% gc_all =[gc12_all; gc21_all]';
% gcm = mean(gc_all,1);%% estimated value mean
% gcer = std(gc_all,[],1)/sqrt(M);%% estimated value se
% nerr=2; %% times of se
% 
% figure;
% b=bar([gcm]', 0.5);
% set(b(1),'FaceColor',[1 1 1]*0.85);
% hold on
% for ii=1:2
%    plot([ii,ii],[gcm(ii)-nerr*gcer(ii),gcm(ii)+nerr*gcer(ii)],'-k','LineWidth',2)
% end
% set(gca,'FontSize',16,'fontWeight','bold')
% xtl={{'LFP->Spike'} {'Spike->LFP'}};
% my_xticklabels(gca,1:2,xtl,'FontSize',16);
% ax=axis;
% axis([ax(1) ax(2) 0 10])
% ylabel('Granger causality')
% title(sprintf('True direction: @%d',index))
% hold on
% plot([1,2],[thresh12,thresh21],'x','color','red')
end
