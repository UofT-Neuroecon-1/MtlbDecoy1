function [ supParticle ] = UpdateSuperParticle( supParticle, SubjData, subj, model, param )
%UPDATESUPERPARTICLE Summary of this function goes here
%   Detailed explanation goes here
if strcmp(model,'Logit')
    %% Updade Hyperparams
    N = numel(supParticle.theta);
    clust_list = reshape([supParticle.theta.clust],param.num_clust,N)'==1;
    r = reshape([supParticle.theta.r],1,N)';
    beta = reshape([supParticle.theta.beta],param.K,N)';
    for c = 1:param.num_clust
        same_clust = clust_list(:,c);
        num_members = sum(same_clust(1:subj),1);
        % Use 5 steps of Gibbs sampling for each Gamma's hyper param
        % Hyperparam beta
        a_0 = supParticle.ha_beta(c,:);
        log_param_p = sum(log(beta(same_clust(1:subj), : )),1);
        sum_theta = sum( beta(same_clust(1:subj), : ),1);
        for i=1:5
            b_0 = gamrnd( 2+(num_members)*a_0 , 1./(1+sum_theta) );
            % Get a_0|b_0 using importance sampling
            draws_a = gamrnd(2,2,50,size(b_0,2));
            log_weights_a = draws_a/2 - log(draws_a) + (draws_a-1) .* log_param_p  ...
                + draws_a .* (2+num_members) .* log(b_0) - (2+num_members) .* log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            a_0 = draws_a(mnrnd(1,weights_a')'==1)';
        end
        supParticle.ha_beta(c,:) = a_0;
        supParticle.hb_beta(c,:) = b_0;
        % Use 5 steps of Gibbs sampling for each Gamma's hyper param
        % Hyperparam R
        a_0 = supParticle.ha_r(c,:);
        log_param_p = sum(log(r(same_clust(1:subj), : )),1);
        sum_theta = sum( r(same_clust(1:subj), : ),1);
        for i=1:5
            b_0 = gamrnd( 2+(num_members)*a_0 , 1./(1+sum_theta) );
            % Get a_0|b_0 using importance sampling
            draws_a = gamrnd(2,2,50,size(b_0,2));
            log_weights_a = draws_a/2 - log(draws_a) + (draws_a-1) .* log_param_p  ...
                + draws_a .* (2+num_members) .* log(b_0) - (2+num_members) .* log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            a_0 = draws_a(mnrnd(1,weights_a')'==1)';
        end
        supParticle.ha_r(c,:) = a_0;
        supParticle.hb_r(c,:) = b_0;
    end
    %% Updade params
    % resample up to 5 subjects' param at random
    mutate_subj = 1:subj;
    if subj > 5
        mutate_subj = randsample(subj,5);
    end
    for s = 1:numel(mutate_subj)
        ss = mutate_subj(s);
        logLikTheta = supParticle.theta(ss).log_like;
        c = supParticle.theta(ss).clust;
        for m = 1 : param.Msteps
            %% joint resampling
            propTheta = supParticle.theta(ss);
            % Beta proposal
            step = gamrnd(100,0.01,1,param.K);
            propTheta.beta = propTheta.beta .* step;
            logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1./step));
            logPriorRatio = (supParticle.ha_beta(c,:)-1) .* log(propTheta.beta/supParticle.theta(ss).beta) ...
                + supParticle.hb_beta(c,:) .* (supParticle.theta(ss).beta - propTheta.beta);
            % r proposal
            step = gamrnd(100,0.01);
            propTheta.r = propTheta.r .* step;
            logQRatio = logQRatio + sum(- 198 .* log(step) + 100 .* (step - 1./step));
            logPriorRatio = logPriorRatio + (supParticle.ha_r(c,:)-1) .* log(propTheta.r/supParticle.theta(ss).r) ...
                + supParticle.hb_r(c,:) .* (supParticle.theta(ss).r - propTheta.r);
            % Compute Prior Ratio
            logLikProp = LogLikelihood( SubjData{ss}, numel(SubjData{ss}.Ys), model , propTheta, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
                supParticle.theta(ss) = propTheta;
                logLikTheta = logLikProp;
            end
        end
        supParticle.theta(ss).log_like = logLikTheta;
    end
elseif strcmp(model,'DN')
    %% Updade Hyperparams
    N = numel(supParticle.theta);
    clust_list = reshape([supParticle.theta.clust],param.num_clust,N)'==1;
    theta = reshape([supParticle.theta.theta],param.size_theta,N)';
    for c = 1:param.num_clust
        same_clust = clust_list(:,c);
        num_members = sum(same_clust(1:subj),1);
        % Use 5 steps of Gibbs sampling for each Gamma's hyper param
        % Hyperparam theta
        a_0 = supParticle.ha_theta(c,:);
        log_param_p = sum(log(theta(same_clust(1:subj), : )),1);
        sum_theta = sum( theta(same_clust(1:subj), : ),1);
        for i=1:5
            b_0 = gamrnd( 2+(num_members)*a_0 , 1./(1+sum_theta) );
            % Get a_0|b_0 using importance sampling
            draws_a = gamrnd(1,2,50,size(b_0,2));
            log_weights_a = draws_a/2 + (draws_a-1) .* (log_param_p+log(0.1))  ...
                + draws_a .* (2+num_members) .* log(b_0) - (2+num_members) .* log(gamma(draws_a));
            log_weights_a = log_weights_a  - max(log_weights_a,[],1 );
            weights_a = exp(log_weights_a)./sum(exp(log_weights_a),1);
            a_0 = draws_a(mnrnd(1,weights_a')'==1)';
        end
        supParticle.ha_theta(c,:) = a_0;
        supParticle.hb_theta(c,:) = b_0;
    end
    %% Updade params
    % resample up to 5 subjects' param at random
    mutate_subj = 1:subj;
    if subj > 5
        mutate_subj = randsample(subj,5);
    end
    for s = 1:numel(mutate_subj)
        ss = mutate_subj(s);
        logLikTheta = supParticle.theta(ss).log_like;
        c = supParticle.theta(ss).clust;
        for m = 1 : param.Msteps
            %% joint resampling
            propTheta = supParticle.theta(ss);
            % theta proposal
            step = gamrnd(100,0.01,1,param.K);
            propTheta.theta = propTheta.theta .* step;
            logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1./step));
            logPriorRatio = (supParticle.ha_theta(c,:)-1) .* log(propTheta.theta/supParticle.theta(ss).theta) ...
                + supParticle.hb_theta(c,:) .* (supParticle.theta(ss).theta - propTheta.theta);
            % Compute Likelihood Ratio
            logLikProp = LogLikelihood( SubjData{ss}, numel(SubjData{ss}.Ys), model , propTheta, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
                supParticle.theta(ss) = propTheta;
                logLikTheta = logLikProp;
            end
        end
        supParticle.theta(ss).log_like = logLikTheta;
    end
elseif strcmp(model,'DN3')
    %% Updade Hyperparams
    clust_list = vertcat(supParticle.theta.clust)==1;
    theta = vertcat(supParticle.theta.theta);
    for c = 1:param.num_clust
        same_clust = clust_list(:,c);
        num_members = sum(same_clust(1:subj),1);
        % Use MH to update hyper param
        log_p_post = -0.6931471806 + sum(log(theta(same_clust(1:subj), : )),1)';
        q_post = 1.5 + sum(theta(same_clust(1:subj), : ),1)';
        r_post = 0.1 + num_members;
        s_post = 0.4 + num_members;
        
        for m = 1:100
            draw_step = gamrnd(100,0.01,param.size_theta,2);
            prop_x = supParticle.h_theta(:,:,c) .* draw_step;
            log_kernel_ratio = sum( -198 .* log(draw_step) + 100 .* (draw_step - 1./ draw_step) , 2);
            log_target_ratio = log_p_post .* (prop_x(:,1) - supParticle.h_theta(:,1,c)) - (prop_x(:,2) - supParticle.h_theta(:,2,c)) .* q_post ...
                + s_post .* (prop_x(:,1) .* log(prop_x(:,2)) -  supParticle.h_theta(:,1,c) .* log(supParticle.h_theta(:,2,c))  ) ...
                - r_post .* (gammaln(prop_x(:,1))-gammaln(supParticle.h_theta(:,1,c) ));
            bool_accept = log(rand(param.size_theta,1)) < (log_kernel_ratio + log_target_ratio);
            supParticle.h_theta(bool_accept,:,c) = prop_x(bool_accept,:);
        end
    end
    %% Updade params
    % resample up to 5 subjects' param at random
    mutate_subj = 1:subj;
    if subj > 5
        mutate_subj = randsample(subj,5);
    end
    for s = 1:numel(mutate_subj)
        ss = mutate_subj(s);
        logLikTheta = supParticle.theta(ss).log_like;
        c = supParticle.theta(ss).clust;
        for m = 1 : param.Msteps
            %% joint resampling
            propTheta = supParticle.theta(ss);
            % theta proposal
            step = gamrnd(100,0.01,1,param.size_theta);
            propTheta.theta = propTheta.theta .* step;
            logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1./step));
            logPriorRatio = sum( (supParticle.h_theta(:,1,c)-1)' .* log(propTheta.theta/supParticle.theta(ss).theta) ...
                + supParticle.h_theta(:,2,c)' .* (supParticle.theta(ss).theta - propTheta.theta));
            % Compute Likelihood Ratio
            logLikProp = LogLikelihood( SubjData{ss}, numel(SubjData{ss}.Ys), model , propTheta, param );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
                supParticle.theta(ss) = propTheta;
                logLikTheta = logLikProp;
            end
        end
        supParticle.theta(ss).log_like = logLikTheta;
    end
else
    error('UpdateSuperParticle : unknown model');
end

end

