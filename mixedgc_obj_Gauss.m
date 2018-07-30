% Demo codes for the Copula-based Granger causality for mixed data (e.g., LFP and Spike)
% 
% Meng Hu, Mingyao Li, Wu Li and Hualou Liang, Joint Analysis of Spikes 
%   and Local Field Potentials using Copula, NeuroImage, 133: 457 ? 467, 2016
%
% Meng Hu @ Liang's lab at Drexel University
%

function [lk der]=mixedgc_obj_Gauss(Y1,Y2,x,porder)
%
% Cost function for optimization for full model
%
% Input:
%     Y1 - Continuous data
%     Y2 - Binary data
%     x - initial parameter
%     porder - model order
%     
% Output:
%     lk - likelihood 
%     der -  derivation


N=length(Y1);
v=length(x);

win=1;
y1(:,1)=Y1(1+porder*win:end);
y2(:,1)=Y2(1+porder*win:end);
ll=length(y1);

for n=1:porder
    for m=1:ll
        y1(m,n+1)=sum(Y1(porder*win+m-1-(n-1)*win:-1:porder*win+m-n*win));
        y2(m,n+1)=sum(Y2(porder*win+m-1-(n-1)*win:-1:porder*win+m-n*win));
    end
end
yhist=[y1(:,2:end) y2(:,2:end)];

mu1=x(1)+yhist*x(2:2*porder+1)';
tmpp=x(2*porder+2)+yhist*x(2*porder+3:v-1)';
mu2=exp(tmpp)./(1+exp(tmpp));
gamma=-1+2*exp(-exp(x(v)));

a0=normpdf(y1(:,1),mu1,1);
a1 = norminv(mu2);
a2 = gamma * (y1(:,1) - mu1);
a3 = sqrt(1 - gamma * gamma);
a4 = normcdf((a1 + a2) / a3,0,1);

LLK=[];
%%
indx21=(y2(:,1)==1);
indxtmp=find(indx21==1);
LLK(indxtmp) = a0(indxtmp) .* a4(indxtmp);

%%
indx20=(y2(:,1)==0);
indxtmp=find(indx20==1);
LLK(indxtmp) = a0(indxtmp) .* (1 - a4(indxtmp));

%% LLK
LLK(isinf(log(LLK))) = .1;
lk=sum(log(LLK));
lk=-lk;


%% derivation

a1=norminv(mu2);
a2=gamma*(y1(:,1)-mu1);
a3=sqrt(1-gamma*gamma);
a4=normpdf((a1+a2)/a3,0,1);
a5=normcdf((a1+a2)/a3,0,1);
a6=normpdf(a1,0,1);


%% Y2=1
indx21=(y2(:,1)==1);
indxtmp=find(indx21==1);

       fvec(indxtmp,1)=(y1(indxtmp,1) - mu1(indxtmp)) -  a4(indxtmp)./a5(indxtmp)*gamma/a3;
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,2*porder+2)=a4(indxtmp)./(a5(indxtmp)*a3.*a6(indxtmp)).*mu2(indxtmp).*(1-mu2(indxtmp));
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-2*porder-2);
       end       
       fvec(indxtmp,v)=a4(indxtmp)./a5(indxtmp) .* ((y1(indxtmp,1) - mu1(indxtmp)) + gamma * a1(indxtmp))/(a3*a3*a3)*(-2*exp(-exp(x(v)))*exp(x(v)));

%% Y2=0       
indx20=(y2(:,1)==0);
indxtmp=find(indx20==1);

       fvec(indxtmp,1)=(y1(indxtmp,1) - mu1(indxtmp)) +  a4(indxtmp)./(1-a5(indxtmp))*gamma/a3;    
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end 
       fvec(indxtmp,2*porder+2)=-a4(indxtmp)./((1-a5(indxtmp))*a3.*a6(indxtmp)).*mu2(indxtmp).*(1-mu2(indxtmp));
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-2*porder-2);
       end      
       fvec(indxtmp,v)=-a4(indxtmp)./(1-a5(indxtmp)) .* ((y1(indxtmp,1) - mu1(indxtmp)) + gamma * a1(indxtmp))/(a3*a3*a3)*(-2*exp(-exp(x(v)))*exp(x(v)));

%%       
fvec_new=fvec;       

if nargout > 1
    der=-sum(fvec_new,1); 
end


end