clear
load('simulout3.mat')
for s=1:40
    for t=1:250
        data(s).X{t}=V(t,:,s)';
    end
    data(s).y=choice{s};
    data(s).J=3;
    data(s).K=1;
end
save ExampleData3 'data' 'par'



clear
load('simuloutSS.mat')

for s=1:30
    for t=1:270
        X=V((s-1)*270+t,:)';
        X(isnan(X)) = [];
        data(s).X{t}=X;
        data(s).J(t)=length(X);
        data(s).K=1;
    end
    data(s).y=choice{1}((s-1)*270+1:s*270);  
end

load trueparSS.mat
par.scale=.5;
par.omega=omega; %Normalized?
par.sigma=sigma; %Gain control
par.a=a;
par.b=b; %power inside sum of denominator
par.kappa=1;

save ExampleDataSS 'data' 'par'

clear
N=[2,4,6,8,10,12];
relative='none';
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
       
        
        %because choice data for alts 3...N not collected, have to match
        %with chosen value instead
        [~,idx]=sort(SetSizeData.data(ind,1:J)==repmat(SetSizeData.data(ind,14),1,J),2,'descend');
        y{s,1}=idx(:,1);
        y{s}(SetSizeData.data(ind,15)==1)=1; %correct for the matching algorithm by using the hard data
    	y{s}(SetSizeData.data(ind,15)==0)=2; %correct for the matching algorithm by using the hard data

        data(s).y=y{s};
        if strcmp(relative,'none')
            data(s).X=V{s}; %normalization done in earlier step
        elseif strcmp(relative,'mean')  
            meanfun=@(x) (x./mean(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
        elseif strcmp(relative,'max')
            maxfun=@(x) (x./max(BDMdata(s).data(:,6)));
            normV{s,1}=cellfun(maxfun,V{s},'uniformoutput',false);
%         elseif strcmp(relative,'kenway')
%             maxfun=@(x) (x./max(BDMdata(s).data(:,6)));
%             normV{s,1}=cellfun(maxfun,V{s},'uniformoutput',false);
        end
        data(s).K=1;
        for t=1:length(data(s).X)
            data(s).J(t)=length(data(s).X{t});
        end
    end

    save ExpDataSS 'data'
    
    clear
    opts.SS=0;
    relative='none';
    J=3; %Number of options
    path='';
    S=40;
    load('datatrinaryconverted.mat');
    load('alldata.mat' ,'alldata','CDpop');
    tind=1:10000;
    olddata=data;
    clear data
    for s=1:S
        T=max(olddata{s}(:,1));
        for t= 1:T
           y{s}(t,1)=find(olddata{s}(olddata{s}(:,1)==t,3));
        end

       V{s}=num2cell(reshape(olddata{s}(:,4),J,T),1);
      
       data(s).y=y{s};
       
       meanV(s)=mean(CDpop.BDM.mean(:,s));
       maxV(s)=max(CDpop.BDM.mean(:,s));
       minV(s)=min(CDpop.BDM.mean(:,s));
       
       if strcmp(relative,'none')
        data(s).X=V{s}; %normalization done in earlier step
%        elseif strcmp(relative,'mean') 
%        meanfun=@(x) (x./meanV(s));
%        normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
%        end
       elseif strcmp(relative,'mean')  
            meanfun=@(x) (x./mean(CDpop.BDM.mean(:,s)));
            normV{s,1}=cellfun(meanfun,V{s},'uniformoutput',false);
        elseif strcmp(relative,'max')
            maxfun=@(x) (x./max(CDpop.BDM.mean(:,s)));
            normV{s,1}=cellfun(maxfun,V{s},'uniformoutput',false);
        end
       scalefun=@(x) (x-minV(s))./(maxV(s)-minV(s));
       %normV{s}=cellfun(scalefun,V{s},'uniformoutput',false);

    end
    i0=3; %normalize w.r.t. option #
      data(s).K=1;
        for t=1:length(data(s).X)
            data(s).J(t)=length(data(s).X{t});
        end

    save ExpData3 'data'
