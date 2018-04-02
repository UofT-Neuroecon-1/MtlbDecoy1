function [ Particles ] = UpdateParticles( Particles, data, subj, obs, model,opts)
%UPDATEPARTICLES Update particles with a resample-move(SMC) algorithm


%% Get particles info

    %% If new subject, draw subject specific particles values
    if obs == 1
        for g=1:opts.G
            for p=1:opts.P
                Particles.particle{g,p} = UpdateNewSubject( Particles.particle{g,p} , model, subj ,opts );
            end
        end
    end
    %% C Phase
    % Reweight each particle to take into account new observation
    logweights = Particles.logweights;
    dataobs.X = data(subj).X(obs);
    dataobs.Mi = data(subj).Mi(obs);
    [dataobs.J, dataobs.K] = size(dataobs.X{1});
    
    par(Particles.LB==Particles.UB)=Particles.r0;%Set the restricted variables.
    index=Particles.LB~=Particles.UB;
    for g = 1:opts.G
        for p = 1:opts.P          

            if any(index)
                par(index)=Particles.particle{g,p}.theta(subj,:);%Set the unrestricted variables to be those passed to the function.
            end
            particles.theta=par;
            proba_choice = ProbaChoice( dataobs,particles, Particles.model,opts );
            logweights(g,p) = logweights(g,p) + log(proba_choice);
            Particles.particle{g,p}.log_like_subj(subj) = Particles.particle{g,p}.log_like_subj(subj) + log(proba_choice);
        end
    end
    % Compute relative ESS
    ress = sum(sum(exp(logweights)))^2 / (opts.G*opts.P*sum(sum(exp(2*logweights))));
    
    % Save marginal likelihood for current observation
    log_w_bar = log( sum(exp(logweights),2) ./ sum(exp(Particles.logweights),2) );%(param.G*param.P));
    Particles.log_marg_like(subj,:) = Particles.log_marg_like(subj,:) + log_w_bar';
    Particles.log_marg_like_total = Particles.log_marg_like_total + log_w_bar;
    Particles.logweights = logweights;
    Particles.logweights = Particles.logweights - max(Particles.logweights,[],2);
    fprintf('End C: Model %s ; RESS = %.5f ; %d - %d subject logML = %.5f\n',model,ress,subj,obs,mean(Particles.log_marg_like(subj,:)) );
    
    %% S Phase : Importance Resampling, within particle groups
    if ~opts.Adaptive || ress < opts.ress_threshold || obs == numel(data(subj).y) || obs < 5
        for g=1:opts.G
            %Resample if weights are not all the same
            if length(unique(logweights(g,:))) > 1
                resample = drawidx(opts.P,exp(logweights(g,:)));
                % Make a copy of the drawn particles
                NewTheta = cell(1,opts.P);
                for p=1:opts.P
                   NewTheta{1,p}=Particles.particle{g,resample(p)};
                end
                % Copy back the drawn particle to the correct position
                for p=1:opts.P
                   Particles.particle{g,p} = NewTheta{1,p};
                end
            end
        end
        Particles.logweights = zeros(opts.G,opts.P);
        clear temp_Particles;
    end
    
    %% M phase
    if ~opts.Adaptive || ress < opts.ress_threshold || obs == numel(data(subj).y) || obs < 5
        accept_count = zeros(opts.G,opts.P);
        for g = 1:opts.G
            algoVars.LB = Particles.LB;
            algoVars.UB = Particles.UB;
            algoVars.r0 = Particles.r0;
            algoVars.model = Particles.model;
            % Get the Cholesky decomposition of the Covariance matrix for each
            % subject's parameters;
            algoVars.chol_cov_theta = CholCovTheta( Particles.particle(g,:), opts );
            
            % Copy particles for parallel looping
            temp_Particles = Particles.particle(g,:);
            parfor p = 1:opts.P
                [temp_Particles{p},accept_count(g,p)] = Mutate(data, subj, obs, model, temp_Particles{p},algoVars,opts);
            end
            Particles.particle(g,:) = temp_Particles;
        end
        fprintf('End M: Model %s ; accept. ratio = %.5f\n',model,mean(squeeze(sum(accept_count,2)) ./ opts.P));
    end
end


