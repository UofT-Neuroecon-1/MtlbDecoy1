function [ logprior ] = logPrior( particle, HyperParams, model , param )
%LOGPRIOR Returns the log prior by approximately integrating out the hyper parameters
    logprior = struct;
    if strcmp(model,'Logit')
        logprior.r = log(sum( ...
            gampdf(particle.r,HyperParams.ha_r(particle.clust,:,:),1./HyperParams.hb_r(particle.clust,:,:)) ...
            ));
        logprior.beta = zeros(1,param.K);
        for k = 1:param.K
            logprior.beta(k) = log(sum( ...
                gampdf(particle.beta(k),HyperParams.ha_beta(particle.clust,k,:),1./HyperParams.hb_beta(particle.clust,k,:)) ...
                )); 
        end
        logprior.total = sum(logprior.beta) + logprior.r;
    elseif strcmp(model,'DN')
        logprior.theta = zeros(1,numel(particle.theta));
        for k = 1:numel(particle.theta)
            logprior.theta(k) = log(sum( ...
                gampdf(particle.theta(k),HyperParams.ha_theta(particle.clust,k,:),1./HyperParams.hb_theta(particle.clust,k,:)) ...
                ));
        end
        logprior.total = sum(logprior.theta);
    else
        error('logPrior : unknown model');
    end
end

