function [ particle , accept ] = Mutate(data, subj, obs, model, particle,algoVars, opts)
%% define mvn random generator from cholesky decomposition
chol_cov_theta= algoVars.chol_cov_theta;
theta_size = size(chol_cov_theta{subj},1);
mvnrnd_chol = @(chol_cov_theta) (chol_cov_theta * randn(theta_size,1))';
%% Number of Msteps varies accrding to observation numbers
% More Msteps are done for the first few observations
if obs < 10
    MSteps = 10 * opts.Msteps;
else
    MSteps = opts.Msteps;
end
%% Mutate particle
accept=0;
if  strcmp(model,'Logit')
    logLikTheta = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), subj , model , particle, opts );
    for m = 1 : MSteps
        % Prior: [betarnd(3,1) gamrnd(2,2,1,K)];
        %% joint resampling
        propTheta = struct;
        propTheta.theta = particle.theta(subj,:);
        step = gamrnd(100,0.01,size(propTheta.theta));
        propTheta.theta = propTheta.theta .* step;
        if all(propTheta.theta > 0) && propTheta.theta(1) <=1
            logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1 ./step));
            logPriorRatio = 2 * log (propTheta.theta(1)/particle.theta(subj,1)) ...
                + sum(particle.theta(subj,2:end) - propTheta.theta(2:end),2) / 2 ...
                + sum(log(propTheta.theta(2:end)./particle.theta(subj,2:end)));
            logLikProp = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), 1 , model , propTheta, opts );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
                accept=accept + 1/MSteps;
                particle.theta(subj,:) = propTheta.theta;
                logLikTheta = logLikProp;
            end
        end
    end
    
elseif  strcmp(model,'DN1')
    logLikTheta = particle.log_like_subj(subj);
    for m = 1 : MSteps
        % Prior: [betarnd(3,1) gamrnd(2,2,1,K)];
        %% joint resampling
        propTheta = struct;
        propTheta.theta = particle.theta(subj,:);
        step = gamrnd(100,0.01,size(propTheta.theta));
        propTheta.theta = propTheta.theta .* step;
        if all(propTheta.theta > 0)
            logQRatio = sum(- 198 .* log(step) + 100 .* (step - 1 ./step));
            logPriorRatio = sum(particle.theta(subj,:) - propTheta.theta,2);
            logLikProp = LogLikelihood( data(subj).X(1:obs), data(subj).y(1:obs), 1 , model , propTheta, opts, algoVars );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta + logQRatio
                accept=accept + 1/MSteps;
                particle.theta(subj,:) = propTheta.theta;
                logLikTheta = logLikProp;
                particle.log_like_subj(subj) = logLikProp;
            end
        end
    end
    
elseif strcmp(model,'PDNNew') || strcmp(model,'RemiStand')
    % Prior: theta = [betarnd(3,1) gamrnd(1,0.5,1,1) gamrnd(1,1,1,K)];
    logLikTheta = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), subj , model , particle, opts );
    for m = 1 : MSteps
        %% joint resampling
        propTheta = struct;
        propTheta.theta = particle.theta(subj,:);
        propTheta.theta = propTheta.theta + mvnrnd_chol(chol_cov_theta{subj}/2);
        if all(propTheta.theta > 0) && propTheta.theta(1) <=1
            logPriorRatio = 2 * log (propTheta.theta(1)/particle.theta(subj,1)) ...
                + (particle.theta(subj,2) - propTheta.theta(2))*2 ...
                + sum(particle.theta(subj,3:end) - propTheta.theta(3:end),2);
            logLikProp = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), 1 , model , propTheta, opts );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                particle.theta(subj,:) = propTheta.theta;
                logLikTheta = logLikProp;
            end
        end
    end

elseif strcmp(model,'HierarchicalProbit')
    %% If it is a new subject,update all previous subjects' parameter
    if obs == 1 && subj > 1
        for ss = 1:subj-1
            num_obs = numel(SubjData{ss}.ChoiceList);
            logLikTheta = LogLikelihood( SubjData{ss}.Xs(1:num_obs), SubjData{ss}.ChoiceList, ss , model , particle, opts );
            %Get sufficient statistics for other subject's parameters
            other_subj_list = (1:size(particle.theta,1))~=ss & (1:size(particle.theta,1)) <= subj;
            q = sum(particle.theta(other_subj_list,:),1);
            q(2) = q(2) + particle.hypertheta(2);
            q(3:end) = q(3:end) + particle.hypertheta(4);
            a_prime = particle.hypertheta([1 3])+subj-1;
            for m = 1 : 2
                %% joint resampling
                propTheta = struct;
                propTheta.theta = particle.theta(ss,:);
                propTheta.theta = propTheta.theta + mvnrnd_chol(chol_cov_theta{ss}/8);
                if all(propTheta.theta > 0) && propTheta.theta(1) <=1
                    logPriorRatio = 2 * log (propTheta.theta(1)/particle.theta(ss,1)) ...
                        - (1+a_prime(1)) * log((q(2)+propTheta.theta(2)) / (q(2)+particle.theta(ss,2)))  ...
                        - sum( (1+a_prime(2)) * log((q(3:end)+propTheta.theta(3:end)) ./ (q(3:end)+particle.theta(ss,3:end)))  );
                    logLikProp = LogLikelihood( SubjData{ss}.Xs(1:num_obs), SubjData{ss}.ChoiceList, 1 , model , propTheta, opts );
                    %accept-reject
                    if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                        particle.theta(ss,:) = propTheta.theta;
                        logLikTheta = logLikProp;
                    end
                end
            end
        end
    end
    %% Update current subject's parameters
    logLikTheta = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), subj , model , particle, opts );
    %Get sufficient statistics for other subject's parameters
    q = sum(particle.theta(1:subj-1,:),1);
    q(2) = q(2) + particle.hypertheta(2);
    q(3:end) = q(3:end) + particle.hypertheta(4);
    a_prime = particle.hypertheta([1 3])+subj-1;
    for m = 1 : MSteps
        %% joint resampling
        propTheta = struct;
        propTheta.theta = particle.theta(subj,:);
        propTheta.theta = propTheta.theta + mvnrnd_chol(chol_cov_theta{subj}/8);
        if all(propTheta.theta > 0) && propTheta.theta(1) <=1
            logPriorRatio = 2 * log (propTheta.theta(1)/particle.theta(subj,1)) ...
                - (1+a_prime(1)) * log((q(2)+propTheta.theta(2)) / (q(2)+particle.theta(subj,2)))  ...
                - sum( (1+a_prime(2)) * log((q(3:end)+propTheta.theta(3:end)) ./ (q(3:end)+particle.theta(subj,3:end)))  );
            logLikProp = LogLikelihood( SubjData{subj}.Xs(1:obs), SubjData{subj}.ChoiceList(1:obs), 1 , model , propTheta, opts );
            %accept-reject
            if log(rand()) <= logPriorRatio + logLikProp - logLikTheta
                accept=accept + 1/MSteps;
                particle.theta(subj,:) = propTheta.theta;
                logLikTheta = logLikProp;
            end
        end
    end
        
else
    error('Mutate : unknown model');
end


end

