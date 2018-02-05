function [ theta , accept ] = Mutate(Xs,ChoiceList, model, theta,chol_cov_theta, attrSign, param)
%% define mvn random generator from cholesky decomposition
theta_size = size(chol_cov_theta,1);
mvnrnd_chol = @(chol_cov_theta) (chol_cov_theta * randn(theta_size,1))';
%% Mutate particle
accept=0;
K = size(Xs{1,1},2);
SeqIndex = length(Xs);
if(SeqIndex<10)
    MSteps = 20*param.Msteps;
elseif(SeqIndex<30)
    MSteps = 10*param.Msteps;
else
    MSteps = 5*param.Msteps;
end
if  strcmp(model,'Logit')
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : MSteps
        % Prior: [betarnd(3,1) gamrnd(2,2,1,K)];
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta);
        if all(propTheta > 0) && propTheta(1) <=1
            logPriorRatio = 2 * log (propTheta(1)/theta(1)) ...
                + sum(theta(2:end) - propTheta(2:end),2) / 2;
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'PDN') || strcmp(model,'RemiStandardized') || strcmp(model,'RemiStandardized2')
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : MSteps
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta);
        if all(propTheta > 0) && propTheta(1) <=1
            logPriorRatio = 2 * log (propTheta(1)/theta(1)) ... %9*log(propTheta(1)/theta(1)) +  (theta(1) - propTheta(1))/0.08 ...
                + (theta(2) - propTheta(2))*2 ...
                + sum(theta(3:end) - propTheta(3:end),2) / 2 ;
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end
elseif strcmp(model,'RemiStandardizedL')
    % Prior constants
        mean_att_val = [15, 0.25, 15, 1, 1, 1, 50];
        sd_att_val = [5, 0.10, 5, 0.5, 0.7, 0.7, 45];
    prior_param = [ mean_att_val sd_att_val];
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : MSteps
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta) / 2;
        %theta = [betarnd(3,1) gamrnd(2,2,1,K) gamrnd(mean_att_val,1) gamrnd(sd_att_val,1)];
        if all(propTheta > 0) && propTheta(1) <=1
            logPriorRatio =  log (propTheta(1)/theta(1)) ... %9*log(propTheta(1)/theta(1)) +  (theta(1) - propTheta(1))/0.08 ...
                + sum(theta(2:K+1) - propTheta(2:K+1),2) / 2 ...
                + sum( (prior_param -1) .*  log(propTheta(K+2:end)./theta(K+2:end) ) + theta(K+2:end) - propTheta(K+2:end)  );
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'PDNNew')
    % Prior: theta = [betarnd(3,1) gamrnd(1,0.5,1,1) gamrnd(1,1,1,K)];
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : MSteps
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta);
        if all(propTheta > 0) && propTheta(1) <=1
            logPriorRatio = 2 * log (propTheta(1)/theta(1)) ... %9*log(propTheta(1)/theta(1)) +  (theta(1) - propTheta(1))/0.08 ...
                + (theta(2) - propTheta(2))*2 ...
                + sum(theta(3:end) - propTheta(3:end),2);
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'PDNUnitIndep')
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : MSteps
        % Prior: [betarnd(3,1) gamrnd(2,4,1,K+1)];
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta);
        if all(propTheta > 0) && propTheta(1) <=1
            logPriorRatio = 2 * log (propTheta(1)/theta(1)) ...
                + (theta(2:2+K-1) - propTheta(2:2+K-1)) *2 ...
                + sum(theta(3:end) - propTheta(3:end),2) / 2 ...
                + sum(log(propTheta(3:end)./theta(3:end)));
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'RemiProbit')
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : (2*MSteps)
        % Prior:[gamrnd(1,0.5,1,K) gamrnd(2,3,1,K) gamrnd(1,1)];
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta) / 4;
        if all(propTheta > 0) && propTheta(end) < 1
            logPriorRatio = (theta(1:K) - propTheta(1:K)) *2 ...
                + sum(theta(K+1:2*K) - propTheta(K+1:2*K),2) / 2 ...
                + sum(log(propTheta(K+1:2*K)./theta(K+1:2*K))) ...
                + 3 * log(propTheta(end) / theta(end)); % alpha ~ Beta(4,1)
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/(2*MSteps);
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'RemiProbitNorm') ||strcmp(model,'RemiProbitStandardized') 
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    for m = 1 : (2*MSteps)
        % Prior: [gamrnd(1,0.5,1,K) gamrnd(1,2,1,K) Gamma(1,1) gamrnd(1,0.5,1,1)];
        %% joint resampling
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta) / 4;
        if all(propTheta > 0)
            logPriorRatio = (theta(1:K) - propTheta(1:K)) *2 ...
                + sum(theta(K+1:2*K) - propTheta(K+1:2*K),2) / 2 ...
                + theta(end-1) - propTheta(end-1) + 2 * (theta(end)-propTheta(end));
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/(2*MSteps);
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end
elseif strcmp(model,'MLBA')
    std_theta_g = (diag(cov_param).^0.5)';
    % Theta: [I0 m lambda beta]
    % Prior: Gam(2,4) Beta(3,1) Gam(2,4) Gam(2,4)
    accept=0;
    logLikTheta = LogLikelihood( Xs, ChoiceList , model , theta, attrSign, param );
    %% joint resampling
    for m=1:(param.Msteps + 4 * (SeqIndex < 10) + 4 * (SeqIndex < 5))
        propTheta = theta;
        propTheta = propTheta + mvnrnd_chol(chol_cov_theta) ./ 2;
        if all(propTheta > 0) && propTheta(2) <=1
            logPriorRatio = propTheta(2) - theta(2) ...
                + sum(theta(3:end)-propTheta(3:end)) ...
                + 2 * log (propTheta(2)/theta(2));
            logLikProp = LogLikelihood( Xs, ChoiceList , model , propTheta, attrSign, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=1;
                theta = propTheta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'MLBAvsPDN')
    accept = 0;
    if theta(1)
        [ subtheta , accept] = Mutate(Xs,ChoiceList, 'PDN', theta(2:end),cov_param{1}, attrSign);
    else
        [ subtheta , accept] = Mutate(Xs,ChoiceList, 'MLBA', theta(2:end),cov_param{2}, attrSign);
    end
    theta(2:end) = subtheta;

else

end


end

