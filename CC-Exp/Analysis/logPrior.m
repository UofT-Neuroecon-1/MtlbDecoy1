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
    elseif strcmp(model,'HBC-ADN')
        %% Compute log marginal prior : \int P(theta|H) P(H) dH
        %% Prior for gamma distributed params
        logprior.gam = zeros(1,param.size_gam);
        for k = 1:param.size_gam
            logprior.gam(k) = log(sum( ...
                gampdf(particle.gam(k),HyperParams.h_gam(k,1,particle.clust,:),1./HyperParams.h_gam(k,2,particle.clust,:)) ...
                ));
        end
        %% Prior for mvn distributed params
        logprior.mvn = sum( log( ...
            mvnpdf(particle.mvn,squeeze(HyperParams.h_mvn_m(:,particle.clust,:))',squeeze(HyperParams.h_mvn_s(:,:,particle.clust,:)))...
        ));
        %% Prior for beta distributed params
        logprior.bet = zeros(1,param.size_bet);
        for k = 1:param.size_bet
            logprior.bet(k) = sum( log( ...
                gampdf(particle.bet(1,k,1),HyperParams.h_bet(k,1,particle.clust,:),ones(1,1,1,param.P)) + ...
                gampdf(particle.bet(1,k,2),HyperParams.h_bet(k,2,particle.clust,:),ones(1,1,1,param.P)) ...
                ));
        end
        %% Total log prior
        logprior.total = sum(logprior.gam) + logprior.mvn + sum(logprior.bet);
            
    elseif strcmp(model,'DN3')
        logprior.theta = zeros(1,param.size_theta);
        for k = 1:param.size_theta
            logprior.theta(k) = log(sum( ...
                gampdf(particle.theta(k),HyperParams.h_theta(k,1,particle.clust,:),1./HyperParams.h_theta(k,2,particle.clust,:)) ...
                ));
        end
        logprior.total = sum(logprior.theta);
    else
        error('logPrior : unknown model');
    end
end

