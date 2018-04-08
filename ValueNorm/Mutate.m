function [ particle , accept ] = Mutate(CurSubjData, obs, model, particle,HyperParams, param)
%% Number of Msteps varies accrding to observation numbers
% More Msteps are done for the first few observations
if obs < 10
    MSteps = 10 * param.Msteps;
else
    MSteps = param.Msteps;
end
%% Mutate particle
accept=0;
if  strcmp(model,'Logit')
    logLikTheta = particle.log_like;
    for m = 1 : MSteps
        %% joint resampling
        propTheta = particle;
        % Beta proposal
        step = gamrnd(100,0.01,1,param.K);
        propTheta.beta = propTheta.beta .* step;
        logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1./step));
        % r proposal
        step = gamrnd(100,0.01);
        propTheta.r = propTheta.r .* step;
        logQRatio = logQRatio + sum(- 198 .* log(step) + 100 .* (step - 1./step));
        % Compute Prior Ratio
        propTheta.logprior = logPrior( propTheta, HyperParams, model, param );
        logPriorRatio = propTheta.logprior.total - particle.logprior.total;
        logLikProp = LogLikelihood( CurSubjData, obs, model , propTheta, param );
        %accept-reject
        if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
            accept=accept + 1/MSteps;
            particle = propTheta;
            logLikTheta = logLikProp;
        end
    end
    particle.log_like = logLikTheta;

elseif  strcmp(model,'DN')
    logLikTheta = particle.log_like;
    for m = 1 : MSteps
        %% joint resampling
        propTheta = particle;
        % Theta proposal
        step = gamrnd(100,0.01,1,param.size_theta);
        propTheta.theta = propTheta.theta .* step;
        logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1./step));
        % Compute Prior Ratio
        propTheta.logprior = logPrior( propTheta, HyperParams, model, param );
        logPriorRatio = propTheta.logprior.total - particle.logprior.total;
        logLikProp = LogLikelihood( CurSubjData, obs, model , propTheta, param );
        %accept-reject
        if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
            accept=accept + 1/MSteps;
            particle = propTheta;
            logLikTheta = logLikProp;
        end
    end
    particle.log_like = logLikTheta;
else
    error('Mutate : unknown model');
end


end

