function [ Particles ] = Assimilate( Particles, SubParticles, subj, param )
%ASSIMILATE Assimilates SubParticles into Particles
for m = 1:numel(param.Models)
    model = param.Models{m};
    %% Compute weights
    for g = 1:param.G
        for sp = 1:param.P
            supParticle = Particles{m}.particle{g,sp};
            log_weights = zeros(1,param.P);
            for p = 1:param.P
                particle = SubParticles{m}(g,p);
                %% Compute log_conditional
                logcondit = struct;
                if strcmp(model,'Logit')
                    logcondit.r = log(gampdf( particle.r, ...
                        supParticle.ha_r(particle.clust,:),1./ supParticle.hb_r(particle.clust,:) ...
                        ));
                    logcondit.beta = zeros(1,param.K);
                    for k = 1:param.K
                        logcondit.beta(k) = log(gampdf( particle.beta(k), ...
                            supParticle.ha_beta(particle.clust,k),1./supParticle.hb_beta(particle.clust,k) ...
                            )); 
                    end
                    logcondit.total = sum(logcondit.beta) + logcondit.r;
                elseif strcmp(model,'DN')
                    logcondit.theta = zeros(1,param.size_theta);
                    for k = 1:param.size_theta
                        logcondit.theta(k) = log(gampdf( particle.theta(k), ...
                            supParticle.ha_theta(particle.clust,k),1./supParticle.hb_theta(particle.clust,k) ...
                            ));
                    end
                    logcondit.total = sum(logcondit.theta);
                else
                    error('Assimilate : unknown model');
                end
                %% Compute weights
                log_weights(p) = logcondit.total-particle.logprior.total;
            end
            %% resample
            draws = drawidx(1,exp(log_weights));
            Particles{m}.particle{g,sp}.theta(subj) = SubParticles{m}(g,draws);
            Particles{m}.particle{g,sp}.theta(subj).logprior.total = Particles{m}.particle{g,sp}.theta(subj).logprior.total + log_weights(draws);
        end
    end

end

