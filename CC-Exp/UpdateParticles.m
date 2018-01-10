function [ Particles ] = UpdateParticles( Particles, Xs, ChoiceList)
%UPDATEPARTICLES Update particles with a resample-move algorithm
% Xs: NumObs x 1 Cells containing X
% X: J x K attribute values
% ChoiceList: NumObs x 1 selected options

%% Get particles info
M = size(Particles,1);
SeqIndex = size(Xs,1);

for m=1:M
    param = Particles{m}.param;
    model = Particles{m}.model;
    attrSign = Particles{m}.attrSign;
    OptimTheta = Particles{m}.OptimTheta;
    NewTheta = coder.nullcopy(OptimTheta);
    %% C Phase
    weights = ones(param.G,param.P);
    X = Xs{end};
    choice = ChoiceList(end);
    for g = 1:param.G
        parfor p = 1:param.P
            proba_choice = ProbaChoice( X , model , OptimTheta{g,p}, attrSign, param );
            weights(g,p) = proba_choice(choice);
        end
    end
    
    % Save marginal likelihood for current observation
    if size(Particles{m}.c_phase_transitions,1) == 0
        Particles{m}.LogMargLik_hist_group(:,SeqIndex) = log((1/param.P) * sum(weights,2)) - sum(Particles{m}.LogMargLik_hist_group(:,1:SeqIndex-1),2);
        Particles{m}.LogMargLik_hist(SeqIndex) = log(sum(sum(weights,2)) / (param.P*param.G)) - sum(Particles{m}.LogMargLik_hist(1,1:SeqIndex-1),2);
    else
        Particles{m}.LogMargLik_hist_group(:,SeqIndex) = log((1/param.P) * sum(weights,2)) - sum(Particles{m}.LogMargLik_hist_group(:,Particles{m}.c_phase_transitions(end)+1:SeqIndex-1),2);
        Particles{m}.LogMargLik_hist(SeqIndex) = log(sum(sum(weights,2)) / (param.P*param.G)) - sum(Particles{m}.LogMargLik_hist(1,Particles{m}.c_phase_transitions(end)+1:SeqIndex-1),2);
    end
    Particles{m}.log_w_bar =[Particles{m}.log_w_bar; log(sum(sum(weights))/(param.G*param.P)) ];
    Particles{m}.c_phase_transitions = [Particles{m}.c_phase_transitions; SeqIndex];
    Particles{m}.log_marg_like = Particles{m}.log_marg_like + Particles{m}.log_w_bar(end);
    fprintf('End: %.0f log(P(y|M)) = %.5f\n',SeqIndex,Particles{m}.log_marg_like);
    fprintf('Average MargLik: %.4f\n',exp(Particles{m}.log_marg_like/SeqIndex));
    
    %% S Phase
    %Resample if weights are not all the same
    resample = 1;
    coder.varsize('resample',[1 inf],[0 1]);
    for g=1:param.G
        if length(unique(weights(g,:))) > 1
            resample = drawidx(param.P,weights(g,:));
            for n=1:param.P
               NewTheta{g,n}= OptimTheta{g,resample(n)};
            end
            for n=1:param.P
               OptimTheta{g,n} = NewTheta{g,n};
            end
        end
    end
    
    %% M phase
    accept = zeros(param.G,param.P);
    vecTheta = vectorizeTheta( OptimTheta );
    cov_theta = zeros(size(vecTheta,3),size(vecTheta,3));
    for g = 1:param.G
        cov_theta = cov(squeeze(vecTheta(g,:,:)))./2;
        %set min step variance
        cov_theta(eye(size(vecTheta,3))==1) = max(diag(cov_theta),0.0001);
        parfor n = 1:param.P
            [OptimTheta{g,n},accept(g,n)] = Mutate(Xs,ChoiceList, model, OptimTheta{g,n},cov_theta,attrSign,param);
        end
    end
    mean(squeeze(sum(accept,2)) ./ param.P)
    vecTheta = vectorizeTheta( OptimTheta );
    vectGroupedTheta = reshape(permute(vecTheta,[2 1 3]),[],4);
    postmeans = squeeze(mean(vecTheta,2))
    % save results
    Particles{m}.OptimTheta = OptimTheta;
    Particles{m}.postmeans = postmeans;
end

end

