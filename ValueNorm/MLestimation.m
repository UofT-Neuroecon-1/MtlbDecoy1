function out=MLestimation(data,par0,opts)

    %the base altenative, i0, matters only for the beta par vector. The chol
    %decomp of the cov matrix is still parameterized in terms of the 1st item.
    
    %If J is not constant on each trial, cannot let the covariance matrix
    %be free. It must be independent.
    

   
    if isfield(data,'W')
        W=data.W;
    else W=[];
    end
    
    model=opts.Models{1};
    
    options.MaxFunEvals = 2000;
    
%     %Set defaults
%     opts.getP=0;
%     opts.numInit=100;
%     opts.names=sprintf('Var %d',1:length(par0));
%     opts.indep=1;
%     opts.R=5000;
    
%     if indep==0
%         alg='fminunc';  
%     else %and homoskedastic
%         %alg='fmincon';
%         alg='fminsearchbnd';
%     end
   % alg='fminsearchbnd';
   alg='fmincon';
   opts.objfun=@LLfun; 

 %% Choice set Properties
   
    T=size(data.X,2);
    Jt=cellfun(@length,data.X,'uniformoutput',true);
    J=max(Jt);
    Q=1;
    
%% Pre-calculate Image Matrices for Chosen Alternative
    temp=eye(J-1); 
    for i=1:J
        M{i}=[temp(:,1:i-1) -1*ones(J-1,1) temp(:,i:J-1)];
    end
    for t=1:T
       Mi{t}=M{data.y(t)}(1:Jt(t)-1,1:Jt(t)); 
    end  
   
    data.Mi=Mi;   
%% Set Restrictions 
opts=setRestrictions(model,Jt,opts);
LB=opts.LB;
UB=opts.UB;
    
 %% Set starting values 
    if isempty(par0)
        display('No Initial parameters specified, starting point is random');

        par0=randn(1, length(W)*(J-1) +Q+ ((J-1)*J/2)-1);
    end
    theta0=par0(:,LB~=UB);

    
if opts.getP %just getting Choice Probs, call LL and exit now
    [~, out]=LLfun(theta0);
return
end
   
 %% Start Estimation   
    disp('Start estimation');
    fprintf('# of Observations: %d \n', length(data.y));
    
    fprintf('Size of initial sample for starting values: %d \n',opts.numInit);
    tic;
    

    

    if strcmp(alg,'optimize') 
        display('Using Simplex Method');

        %Use optimize toolbox
        %options.Display='plot';
        options.MaxIter=10000;
        options=optimset(options,'Display',opts.Display);

%         optfun=@optimize;
%         varin={[],[],[],[],[],[],[],[],options,alg};

         time=('00:10:00');
        getHess=0;
    elseif strcmp(alg,'fminsearch')
         display('Using Simplex Method');
         options=optimset(options,'Display',opts.Display);
 %        optfun=@fminsearch;
 %        varin={options};
          
    elseif strcmp(alg,'fminsearchbnd')
         display('Using Bounded Simplex Method');
         options=optimset(options,'Display',opts.Display);
         
          optfun=@fminsearchbnd;
 %       varin={LB(LB~=UB),UB(LB~=UB),options};
        
    elseif strcmp(alg,'fminsearchcon')
         display('Using (Constrained) Simplex Method');
         options=optimset(options,'Display',opts.Display);

%          optfun=@fminsearchcon;
%        varin={LB(LB~=UB),UB(LB~=UB),LIN(LB~=UB),0,[],options};

    elseif strcmp(alg,'fminunc')
         display('Using Unconstrained Quasi-Newton');
         options=optimset(options,'Largescale','off','Display',opts.Display);
 %        optfun=@fminunc;
 %        varin={options};

         getHess=1;

    elseif strcmp(alg,'fmincon')
         display('Using Constrained Quasi-Newton');
         options=optimset(options,'Algorithm','interior-point','Display',opts.Display);
        optfun=@fmincon;
%          varin={[],[],LIN(LB~=UB),0,LB(LB~=UB),UB(LB~=UB),[],options};
%        varin={[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options};

         getHess=1;

    end

    display(sprintf('TolX: %g   TolFun: %g',options.TolX,options.TolFun));

%     if isempty(theta0)
%         -LLfun(LB)
%     end
    %matlabpool open
    if opts.numInit==1
        %[thetah, maxLL, exitflags]=feval(optfun,opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        %[thetah, maxLL, exitflags]=fminsearchbnd(opts.objfun,theta0,LB(LB~=UB),UB(LB~=UB),options);
        %[thetah, maxLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,thetah,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        
        [thetah, maxLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        i=1;

        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-maxLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
    else
        [thetah, nLL, exitflags, xstart]=rmsearch(opts.objfun,'fminsearchbnd',theta0,LB(LB~=UB),UB(LB~=UB),'options',options,'InitialSample',opts.numInit);
            [maxLL,i]=max(-nLL);

        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-nLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
        save 'mpnormEst1stStage.mat' 'thetah' 'nLL' 'exitflags'

        for n=1:opts.numInit
            [thetah(n,:), nLL(n), exitflags(n),~,~,grad(n,:),hess]=fmincon(opts.objfun,thetah(n,:),[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        end
        [maxLL,i]=max(-nLL);
 
        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-nLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
        save 'mpnormEst2ndStage.mat' 'thetah' 'nLL' 'exitflags' 'grad'
    end
                       
%         if strcmp(alg,'fminsearch') || strcmp(alg,'fminsearchcon') || strcmp(alg,'fminsearchbnd') || strcmp(alg,'optimize')
%             [thetah(s,:), nLL(s,1), exitflags(s,1)]=optfun(objfun,theta0(s,:),varin{:});
%             grad=[]
%             hess=[];
%         else
%                 if strcmp(alg,'interior-point')
%                     
%                     [xfinal,ffinal,exitflag,xstart] = rmsearch(fun,optname,x0,LB(LB~=UB),UB(LB~=UB),options)
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),xstart]=rmsearch(objfun,optfun,theta0(s,:),LB(LB~=UB),UB(LB~=UB),options);
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),~,~,grad,H(:,:,s)]=optfun(objfun,theta0(s,:),varin{:});
%                 else
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),~,grad,H(:,:,s)]=optfun(objfun,theta0(s,:),varin{:});
%                 end
% 
%         end


    %matlabpool close
   

    %Get likelihood at estimates for Vuong test	
    [~,P]=opts.objfun(thetah(i,:));
    
    %Get gradient for convergence check 
    %Or just report gradient from optimizer
    %grad=gradest(@(x) -LLfun(x),thetah(i,:));    
    
    checkConvergence;


    % if PREDICT == 1;
    %     disp(' ');
    %     disp('Predict shares at estimated coefficients.');
    %     probs=pred(paramhat);
    % end;
    i0=opts.i0;      
    save 'mpestimates.mat' 'data' 'grad' 'Q' 'J' 'i0'

    display('Estimates saved to disk');
    
    if opts.ses==1
    disp('Calculating finite-difference hessian and taking inverse for standard errors.');
    H=hessian(@(x) -LLfun(x),thetah(i,:));
    covh=inv(-H);
    se=sqrt(diag(covh))';

    ses(LB~=UB)=se;
    ses(LB==UB)=0;
    else
        H=0;
        covh=inv(-H);
       ses=zeros(1,length(LB));
    end
    
    save 'mpestimates.mat' 'par' 'data' 'H' 'grad' 'Q' 'J' 'ses' 'i0'
    

    out.exitflag=exitflag;
    out.theta0=theta0;
    out.max=i;
    out.parh=par;
    out.LL=maxLL;
    out.P=P;
    out.se=ses;
    out.cov=covh;
    out.LB=LB;
    out.UB=UB;
    out.model=model;
    out.Prob=opts.Prob;
        
    %%%%%% Nested Functions %%%%%%%   
    function [nLL, Pi]=LLfun(theta)
        
     
        particle.theta=theta;
        Pi=ProbaChoice(data, particle, opts );

        nLL=-sum(log(Pi));
    end
            
    function checkConvergence
        disp(' ');

        if exitflag == 1
          disp('Convergence achieved.');

        elseif exitflag == 2
          disp('Convergence achieved by criterion based on change in parameters.');
          if size(options.TolX,1)>0
             disp(['Parameters changed less than PARAMTOL= ' num2str(options.TolX)]);
          else
             disp('Parameters changed less than PARAMTOL=0.000001, set by default.');
          end
          disp('You might want to check whether this is actually convergence.');
          disp('The gradient vector is');
          disp(grad)
        elseif exitflag == 3
          disp('Convergence achieved by criterion based on change in log-likelihood value.');
          if size(options.TolFun,1)>0
             disp(['Log-likelihood value changed less than LLTOL= ' num2str(options.TolFun)]);
          else
             disp('Log-likelihood changed less than LLTOL=0.000001, set by default.');
          end
         disp('You might want to check whether this is actually convergence.');
         disp('The gradient vector is');
         disp(grad)

        elseif exitflag == 5
            disp('Predicted decrease in the objective function was less than the TolFun tolerance.');
            disp('Convergence likely achieved, but the first-order optimality measure was above the function tolerance.');

        else
            disp('Convergence not achieved.');
            disp('Results are not printed because no convergence.');
            return
        end


    end
end