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
    weights = Particles.weights;
    dataobs.X = data(subj).X(obs);
    dataobs.Mi = data(subj).Mi(obs);
    [dataobs.J dataobs.K] = size(dataobs.X{1});
    choice = data(subj).y(obs);
    
    par(Particles.LB==Particles.UB)=Particles.r0;%Set the restricted variables.
    index=Particles.LB~=Particles.UB;
    for g = 1:opts.G
        for p = 1:opts.P          

            if any(index)
                par(index)=Particles.particle{g,p}.theta(subj,:);%Set the unrestricted variables to be those passed to the function.
            end
            particles.theta=par;
            proba_choice = ProbaChoice( dataobs,particles, Particles.model,opts );
            weights(g,p) = weights(g,p) * proba_choice(choice);
        end
    end
    % Compute relative ESS
    ress = sum(sum(weights))^2 / (opts.G*opts.P*sum(sum(weights.^2)));
    
    % Save marginal likelihood for current observation
    log_w_bar = log( sum(weights,2) ./ sum(Particles.weights,2) );%(param.G*param.P));
    Particles.log_marg_like(subj,:) = Particles.log_marg_like(subj,:) + log_w_bar';
    Particles.log_marg_like_total = Particles.log_marg_like_total + log_w_bar;
    Particles.weights = weights;
    fprintf('End C: Model %d ; RESS = %.5f ; %d - %d subject avg logML = %.5f\n',model,ress,subj,obs,mean(Particles.log_marg_like(subj,:)) );
    
    %% S Phase : Importance Resampling, within particle groups
    if ~opts.Adaptive || ress < opts.ress_threshold || obs == numel(data(subj).y) || obs < 5
        for g=1:opts.G
            %Resample if weights are not all the same
            if length(unique(weights(g,:))) > 1
                resample = drawidx(opts.P,weights(g,:));
                % Make a copy of the drawn particles
                temp_Particles = cell(1,opts.P);
                for p=1:opts.P
                   NewTheta{1,p}=Particles.particle{g,resample(p)};
                end
                % Copy back the drawn particle to the correct position
                for p=1:opts.P
                   Particles.particle{g,p} = NewTheta{1,p};
                end
            end
        end
        Particles.weights = ones(opts.G,opts.P);
        clear temp_Particles;
    end
    
    %% M phase
    if ~opts.Adaptive || ress < opts.ress_threshold || obs == numel(data(subj).y) || obs < 5
        accept_count = zeros(opts.G,opts.P);
        for g = 1:opts.G
            % Get the Cholesky decomposition of the Covariance matrix for each
            % subject's parameters;
            chol_cov_theta = CholCovTheta( Particles.particle(g,:), opts );
            % Copy particles for parallel looping
            temp_Particles = Particles.particle(g,:);
            parfor p = 1:opts.P
                [temp_Particles{p},accept_count(g,p)] = Mutate(data, subj, obs, model, temp_Particles{p},chol_cov_theta,opts);
            end
            Particles.particle(g,:) = temp_Particles;
        end
        fprintf('End M: Model %d ; accept. ratio = %.5f\n',m,mean(squeeze(sum(accept_count,2)) ./ opts.P));
    end
end


