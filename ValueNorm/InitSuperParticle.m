function [ particle ] = InitSuperParticle( m , param)
%INITPARTICLE Returns a Theta particle drawn from prior
% Get model info
model = param.Models{m};
K = param.K;
N = param.num_subj;
particle = struct;
%model : 'logit', 'MLBA, 'PDN' or 'PDNUnitIndep'
%K : number of attributes
    if  strcmp(model,'Logit')
        particle.ha_r = nan(param.num_clust,1);
        particle.hb_r = gamrnd(2,1,param.num_clust,1);
        particle.ha_beta = nan(param.num_clust,K);
        particle.hb_beta = gamrnd(2,1,param.num_clust,K);
        % Sample ha_beta via importance resampling
        % Hyperparam: a=1, b=2, c=2
        for c = 1:param.num_clust
            draws_a = gamrnd(2,2,50,K);
            log_weights_a = draws_a/2 - log(draws_a) ...
                + 2 .* draws_a .* log(particle.hb_beta(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_beta(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        % Sample ha_r via importance resampling
        % Hyperparam: a=1, b=2, c=2
        for c = 1:param.num_clust
            draws_a = gamrnd(2,2,50,1);
            log_weights_a = draws_a/2 - log(draws_a) ...
                + 2 .* draws_a .* log(particle.hb_r(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_r(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        
        % Init subparticle container
        particle.theta = struct;
        for subj = 1:N
            particle.theta(subj).clust = nan(1,param.num_clust);
            particle.theta(subj).beta = nan(1,K);
            particle.theta(subj).r = nan;
            particle.theta(subj).log_like = 0;
            particle.theta(subj).logprior = 0;
        end
        
    elseif strcmp(model,'HBC-PNE')
        particle.ha_r = nan(param.num_clust,1);
        particle.hb_r = gamrnd(2,1,param.num_clust,1);
        particle.ha_omega = nan(param.num_clust,K);
        particle.hb_omega = gamrnd(2,1,param.num_clust,K);
        particle.ha_sig = nan(param.num_clust,K);
        particle.hb_sig = gamrnd(2,1,param.num_clust,K);
        % Sample ha_omega via importance resampling
        % Hyperparam: a=0.1, b=2, c=2
        for c = 1:param.num_clust
            draws_a = gamrnd(1,2,50,K);
            log_weights_a = draws_a/2 ...
                +(draws_a - 1) .* log(0.1) + 2 .* draws_a .* log(particle.hb_omega(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_omega(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        % Sample ha_sig via importance resampling
        for c = 1:param.num_clust
            draws_a = gamrnd(1,2,50,K);
            log_weights_a = draws_a/2 ...
                +(draws_a - 1) .* log(0.1) + 2 .* draws_a .* log(particle.hb_sig(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_sig(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        % Sample ha_r via importance resampling
        % Hyperparam: a=1, b=2, c=2
        for c = 1:param.num_clust
            draws_a = gamrnd(2,2,50,1);
            log_weights_a = draws_a/2 - log(draws_a) ...
                +(draws_a - 1) .* log(0.1) + 2 .* draws_a .* log(particle.hb_r(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_r(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        
        % Init subparticle container
        particle.theta = struct;
        for subj = 1:N
            particle.theta(subj).clust = nan(1,param.num_clust);
            particle.theta(subj).omega = nan(1,K);
            particle.theta(subj).sig = nan(1,K);
            particle.theta(subj).r = nan;
            particle.theta(subj).log_like = 0;
            particle.theta(subj).logprior = 0;
        end
    elseif strcmp(model,'DN')
        particle.ha_theta = nan(param.num_clust,param.size_theta);
        particle.hb_theta = gamrnd(2,1,param.num_clust,param.size_theta);
        % Sample ha_omega via importance resampling
        % Hyperparam: a=0.1, b=2, c=2
        % Sample ha_theta via importance resampling
        for c = 1:param.num_clust
            draws_a = gamrnd(1,2,50,2);
            log_weights_a = draws_a/2 ...
                +(draws_a - 1) .* log(0.1) + 2 .* draws_a .* log(particle.hb_theta(c,:)) - 2.*log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            particle.ha_theta(c,:) = draws_a(mnrnd(1,weights_a')'==1)';
        end
        
        % Init subparticle container
        particle.theta = struct;
        for subj = 1:N
            particle.theta(subj).clust = nan(1,param.num_clust);
            particle.theta(subj).theta = nan(1,param.size_theta);
        end
    elseif strcmp(model,'DN3')
        log_p = log(0.5);
        q = 1.5;
        r = 0.1;
        s = 0.4;
        n_draws = 1000;
        particle.h_theta = nan(param.size_theta,2,param.num_clust);
        for c = 1:param.num_clust
            % Sample from prior via importance resampling + MH
            draws = gamrnd(2,3,n_draws,2);     
            %weighting draws
            log_weights = (draws(:,1)-1) .* log_p - draws(:,2) .* q + (draws(:,1).*s) .* log(draws(:,2)) - r .* gammaln(draws(:,1));
            log_weights = log_weights - max(log_weights);
            weights = exp(log_weights) ./ sum(exp(log_weights));
            weights = weights ./ prod(gampdf(draws,2,3),2);
            % sampling
            sample_idx = randsample(1:n_draws,param.size_theta,true,weights);
            prior_sample = draws(sample_idx,:);
            
            % Metropolis Hastings step
            for m = 1:100
                draw_step = gamrnd(100,0.01,param.size_theta,2);
                prop_x = particle.h_theta(:,:,c) .* draw_step;
                log_kernel_ratio = sum( -198 .* log(draw_step) + 100 .* (draw_step - 1./ draw_step) , 2);
                log_target_ratio = log_p .* (prop_x(:,1) - prior_sample(:,1)) - (prop_x(:,2) - prior_sample(:,2)) .* q ...
                    + s .* (prop_x(:,1) .* log(prop_x(:,2)) -  prior_sample(:,1) .* log(prior_sample(:,2))  ) ...
                    - r .* (gammaln(prop_x(:,1))-gammaln(prior_sample(:,1) ));
                bool_accept = log(rand(param.size_theta,1)) < (log_kernel_ratio + log_target_ratio);
                prior_sample(bool_accept,:) = prop_x(bool_accept,:);
            end
            particle.h_theta(:,:,c) = prior_sample;
        end
        
        % Init subparticle container
        particle.theta = struct;
        for subj = 1:N
            particle.theta(subj).clust = nan(1,param.num_clust);
            particle.theta(subj).theta = nan(1,param.size_theta);
        end
    else
        error(sprintf('InitParticle: unknown model used (%s)',model));
    end
    % Create likelihood and prior fields
    for subj = 1:N
        particle.theta(subj).log_like = 0;
        particle.theta(subj).logprior = 0;
    end

end

