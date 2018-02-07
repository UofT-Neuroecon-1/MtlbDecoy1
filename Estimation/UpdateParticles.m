function [ Particles ] = UpdateParticles( Particles, SubjData, subj, obs, param)
%UPDATEPARTICLES Update particles with a resample-move(SMC) algorithm


%% Get particles info
M = numel(param.Models);
for m=1:M
    model =  param.Models{m};
    %% If new subject, draw subject specific particles values
    if obs == 1
        for g=1:param.G
            for p=1:param.P
                Particles{m}.particle{g,p} = UpdateNewSubject( Particles{m}.particle{g,p} , model, subj ,param );
            end
        end
    end
    %% C Phase
    % Reweight each particle to take into account new observation
    weights = ones(param.G,param.P);
    ChoiceSet = SubjData{subj}.Xs{obs};
    choice = SubjData{subj}.ChoiceList(obs);
    for g = 1:param.G
        parfor p = 1:param.P
            proba_choice = ProbaChoice( ChoiceSet, subj, model , Particles{m}.particle{g,p}, param );
            weights(g,p) = proba_choice(choice);
        end
    end
    
    % Save marginal likelihood for current observation
    log_w_bar = log(sum(sum(weights))/(param.G*param.P));
    Particles{m}.log_marg_like(subj) = Particles{m}.log_marg_like(subj) + log_w_bar;
    Particles{m}.log_marg_like_total = Particles{m}.log_marg_like_total + log_w_bar;
    fprintf('End C: Model %d ; %d - %d log(P(y|M)) = %.5f\n',m,subj,obs,Particles{m}.log_marg_like(subj));
    
    %% S Phase : Importance Resampling, within particle groups
    for g=1:param.G
        %Resample if weights are not all the same
        if length(unique(weights(g,:))) > 1
            resample = drawidx(param.P,weights(g,:));
            % Make a copy of the drawn particles
            temp_Particles = cell(1,param.P);
            for p=1:param.P
               NewTheta{1,p}=Particles{m}.particle{g,resample(p)};
            end
            % Copy back the drawn particle to the correct position
            for p=1:param.P
               Particles{m}.particle{g,p} = NewTheta{1,p};
            end
        end
    end
    clear temp_Particles;
    
    %% M phase
    accept_count = zeros(param.G,param.P);
    for g = 1:param.G
        % Get the Cholesky decomposition of the Covariance matrix for each
        % sumbject's parameters;
        chol_cov_theta = CholCovTheta( Particles{m}.particle(g,:), param );
        % Copy particles for parallel looping
        temp_Particles = Particles{m}.particle(g,:);
        parfor p = 1:param.P
            [temp_Particles{p},accept_count(g,p)] = Mutate(SubjData, subj, obs, model, temp_Particles{p},chol_cov_theta,param);
        end
        Particles{m}.particle(g,:) = temp_Particles;
    end
    fprintf('End M: Model %d ; accept. ratio = %.5f\n',m,mean(squeeze(sum(accept_count,2)) ./ param.P));
end

end

