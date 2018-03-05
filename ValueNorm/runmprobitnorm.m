 addpath(genpath('../../routines'))
    %%%%run mprobit on set size data
clear



K=0; %Trial Specific Regressor / Alt specific parameter
Q=1; %Alternative Specific regressor 
D=4; %Denominator Parameters


%Use normalized model or not
normalized=1;
range=1; %1 if range normalization
%Calculate ses from (standalone) finite difference hessian, or use hessian
%from optimization routine
calcSEsInd=0;

relative='none';
%relative='mean';
%relative='max';



%type='setsize';
 type='trinary';
%type='simul';
%type='simulss';
%
% Which choice sets in analysis?
N=[2,4,6,8,10,12];
%N=12;


    opts.getP=0;
    opts.numInit=1; %1: run using gradient- method first at theta 0. If >1, use random starting points with gradient free method.
    names={'kappa','G0','omega','a','b','w2'};
    opts.R=5000;
    opts.ses=calcSEsInd;
    opts.i0=3;
    opts.Display='final';
    opts.Sinz=0;
    opts.weights=0;
    opts.range=range; %if 1, range normalization. If 0, div norm
    opts.indep=0;


if strcmp(type,'setsize')

    path='';
    S=30;
    load('SetSizeData','SetSizeData');
    load('SetSizeBDMData');

    for s=1:S
        ind=find(SetSizeData.data(:,17)==s);
        T=max(SetSizeData.data(ind,16));

    if length(N)==1
        tind=find(SetSizeData.data(:,13)==N);
        opts.indep=0;
    else
        tind=1:8100;
        opts.indep=1;
    end
    J=max(N); %Number of options
    
        V{s}=num2cell(SetSizeData.data(ind,1:J)',1);
        for t=1:T %Strip out the NaNs and put in cell
           temp=V{s}{:,t};
           V{s}{:,t}=temp(~isnan(V{s}{:,t}));
        end
        
        X{s,1}=zeros(T,K); 
        
        %because choice data for alts 3...N not collected, have to match
        %with chosen value instead
        [~,idx]=sort(SetSizeData.data(ind,1:J)==repmat(SetSizeData.data(ind,14),1,J),2,'descend');
        y{s,1}=idx(:,1);
        y{s}(SetSizeData.data(ind,15)==1)=1; %correct for the matching algorithm by using the hard data
    	y{s}(SetSizeData.data(ind,15)==0)=2; %correct for the matching algorithm by using the hard data

        if strcmp(relative,'none')
            normV{s}=V{s}; %normalization done in earlier step
        elseif strcmp(relative,'mean')  
            meanfun=@(x) (x./mean(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
        elseif strcmp(relative,'max')
            maxfun=@(x) (x./max(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(maxfun,V{s},'uniformoutput',false);
        end
        
    end
    opts.SS=1;  
    
     i0=2;

elseif strcmp(type,'trinary')
    opts.SS=0;
    J=3; %Number of options
    path='';
    S=40;
    load('datatrinaryconverted.mat');
    load('alldata.mat' ,'alldata','CDpop');
    tind=1:10000;
    
    
    for s=1:S
        T=max(data{s}(:,1));
        for t= 1:T
           y{s}(t,1)=find(data{s}(data{s}(:,1)==t,3));
        end
       X{s,1}=zeros(T,K); 
       V{s}=num2cell(reshape(data{s}(:,4+K:4+K+Q-1),J,T),1);
       
       meanV(s)=mean(CDpop.BDM.mean(:,s));
       maxV(s)=max(CDpop.BDM.mean(:,s));
       minV(s)=min(CDpop.BDM.mean(:,s));
       
       if strcmp(relative,'none')
        normV{s}=V{s}; %normalization done in earlier step
%        elseif strcmp(relative,'mean') 
%        meanfun=@(x) (x./meanV(s));
%        normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
%        end
       elseif strcmp(relative,'mean')  
            meanfun=@(x) (x./mean(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
        elseif strcmp(relative,'max')
            maxfun=@(x) (x./max(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(maxfun,V{s},'uniformoutput',false);
        end
       scalefun=@(x) (x-minV(s))./(maxV(s)-minV(s));
       %normV{s}=cellfun(scalefun,V{s},'uniformoutput',false);

    end
    i0=3; %normalize w.r.t. option #

    
elseif strcmp(type,'simul')
    opts.SS=0;
    J=3; %Number of options
    %%%Run simulated data
     path='../Simulated/simuldata/';
     load([path 'simulout.mat']);

     T=length(choice{1});
     S=size(choice,1);
     %%%%
        y=choice;
        
    VV=V;
    clear V
    for s=1:S
        T=length(choice{s});
       X{s,1}=zeros(T,K); 

       V{s,1}=VV(:,:,s);
       normV{s,1}=V{s};
       
    end
    clear VV
 elseif strcmp(type,'simulss')
    J=12; %Number of options
    %%%Run simulated data
    load('SetSizeData','SetSizeData');
     path='../Simulated/simuldata/setsize/sigma0/';
     load([path 'simulout.mat']);

%      T=length(choice{1});
%      S=size(choice,1);
     T=270;
     S=30;
     %%%%
        

    VV=permute(reshape(V',[12 T S]),[2,1,3]);
    clear V    
    for s=1:S
        for t=1:T %Strip out the NaNs and put in cell
           V{s}{:,t}=VV(t,~isnan(VV(t,:,s)),s)';
        end
    end
    
    temp=reshape(cell2mat(choice),T,S)
    for s=1:S
        y{s}=temp(:,s);
        X{s,1}=zeros(T,K); 
       normV{s,1}=V{s};       
    end
    clear VV
    opts.SS=1; 
else
    error('Wrong Data')
end

% True Pars and/or hypothesis testing. 

M1=[-1*ones(J-1,1) eye(J-1)];
L=chol(M1*0.5*eye(J)*M1')';
L2=L(L~=0);
cind=L2(2:end)';
    
%[1 
% c1 c3
% c2 c4 c5]% 

%%%%Params
a=1;
b=1;
% g=b;
w2=1;

kappa0=1;
G0=1;
omega0=1;

par0=[kappa0 G0 omega0 a b w2 cind ];
%par0=[1 0.5793 0.0172 1 3.1147 1 cind];
%par0=[1 0.00001 0.4319 1 17.46 0.0592 cind];
%par0=[1 0.00001 .4423 1 18.8453 cind];
%par0=[1 0.9236 0 1 1 0.01 cind];
%par0=[1 0.00001 .5399 1.17 26.86 cind];
%par0=[1 0.013    .410 1 1 cind];
%par0=[1 1.4965    0 1 1 0 0.7907 0.4353];
%par0=[1 0.00001 .4423 1 511 cind];

% load('../simuldata/PooledEstimatesw1G1b1.mat','par');
% par0=par;

%2*(LLfun([0.0130    0.4100  206.42])-LLfun([0.0130    0.4100  511])) to
%get standard error of beta

%Restrictions
%sigma
GLB=0;
GUB=inf;

%omega
wLB=0;
wUB=inf;

%kappa
kappaLB=1;
kappaUB=1;

aLB=1;
aUB=1;

bLB=1;
bUB=1;

 w2LB=0;
 w2UB=0; 
 
if opts.indep==1  
    cLB=cind;
    cUB=cind;
else

    temp=-inf(J-1);
    temp=tril(temp,-1);
    temp=temp(L~=0)';
%     
     %cLB=[-inf, 0];
     cLB=temp(2:end);
    cUB=inf(1,length(cind));
end

LB=[-inf(1,K) kappaLB GLB wLB aLB bLB w2LB cLB];
UB=[inf(1,K) kappaUB GUB wUB aUB bUB w2UB cUB];


    
% %Use weighted linear LLfun
% opts.weights=1;    
% %sigma
% GLB=0;
% GUB=inf;    
% %omega
% wLB=[0 0 0];
% wUB=[inf inf inf];    
% 
% par0=[0 1 1 1];
% 
% LB=[GLB wLB];
% UB=[GUB wUB];
    
    

% LB=[-inf(1,K) kappaLB GLB wLB aLB bLB gLB cLB];
% UB=[inf(1,K) kappaUB GUB wUB aUB bUB gUB cUB];
%
%LIN=[0 0 0 0 1 -1 0 0]; %Linear constaints on free parameters. Doesn't
%work yet.

%For testing likelihood ratios, restrict to estimates?
% load('../simuldata/Estimatesw1G1b1Restwb.mat','LL','Gh');
% GhRwb=Gh;
% load('../simuldata/Estimatesw1G1b1RestG.mat','LL','Gh');
% GhRG=Gh;

% load('mpestimates.mat','par');
% par0=par;

% parfor s=1:S 
%     s   
% %     par0s=[par0; kappa0 GhRwb(s) omega0 a b cind; kappa0 GhRG(s) omega0 a b cind];
%     par0s=par0;
%     
%     if normalized==1
%      out{s}=mprobitnormalized(y{s},X{s},V{s},par0s,LB,UB,opts);
%      %out{s}=mprobitnormalized(y{s},X{s},V{s},par0s,LB,UB,LIN,indep,i0,J,R,calcSEsInd,names,0);
%     else
%      out{s}=mprobit(y{s},X{s},V{s},[kappa0 cind],indep,i0,J,R,caclSEsInd);
%     end
%     %save([path 'out.mat'],'out') %Save after each s in case of crash
% end
% save([path 'out.mat'],'out')

% for n=1:100
%     ind=randsample(S,30);
%     out{n}=mprobitnormalized(cell2mat(y(ind)'),cell2mat(X(ind)),cell2mat(normV(ind)),par0,LB,UB,indep,i0,J,R,calcSEsInd,names,0);
%     save([path 'out.mat'],'out') %Save after each s in case of crash
% end

% Must normalize by mean bid, need it from Kenway's data.
yc=cell2mat(y(:));
Vc=[normV{:}];

%% Run pooled or Subject Specific

pool=1;

for s=1:S
    SubjData{s}.Xs=normV{s};
    SubjData{s}.ChoiceList=y{s};
end

if pool==1
%outall=mprobitnormalized(yc(tind),cell2mat(X),Vc(tind),par0,LB,UB,opts);
EstimationSMC( SubjData, param, backup_file )
displayResults(outall,names,par0,J,i0);
save([path 'outall.mat'],'outall')
else


parfor s=1:S 
%     s   

     par0s=par0;
     

      out(s)=mprobitnormalized(y{s},[],normV{s},par0s,LB,UB,opts);%     %save([path 'out.mat'],'out') %Save after each s in case of crash
 end
save([path 'out.mat'],'out')
end
 
 
% load sigOmegers
% outsig=mprobitnormalized(cell2mat(y(sig)'),cell2mat(X),cell2mat(normV(sig)),par0,LB,UB,opts);
% save([path 'outsig.mat'],'outsig')

 %%
clear X
if N==2
    X=-diff([Vc{tind}])';
    sumV=sum([Vc{tind}])';
GLM=fitglm(X,2-yc(tind),'distribution','binomial','link','probit','intercept',false)

F=GLM.Fitted.Probability;
Xb=GLM.Fitted.LinearPredictor;
f=normpdf(Xb);
Var=(1./sqrt(F.*(1-F)));
BRNR=fitglm([Var.*f.*X,Var.*f.*Xb.*sumV],Var.*(2-yc(tind) - F),'intercept',false)
display('Null of homoskedasticity (p-val):'); 1-chi2cdf(BRNR.SSR,1)
end