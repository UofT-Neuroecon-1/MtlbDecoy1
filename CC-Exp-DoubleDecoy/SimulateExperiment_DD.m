% clear all;
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
global attrVals attrNames attrSign
attrVals = cell(2,1);
attrNames = cell(2,1);
attrVals{1}= (5:50);          attrNames{1}="Annual fee ($)";
attrVals{2}= 0.20:0.1:5;    attrNames{2}="Reward (%)";
attrSign = [-1, 1];
num_option_list = [2*ones(10,1);3*ones(20,1)];
opt_num_quest = numel(num_option_list);
num_double_decoy = 25;
list_double_decoy = opt_num_quest + num_double_decoy + randsample(num_double_decoy,num_double_decoy);
num_non_check = opt_num_quest + 2*num_double_decoy;
num_consistency_check = 20;
list_const_check = datasample(opt_num_quest+1:num_non_check,num_consistency_check,'Replace',false);
debug_mode = false;

%% Simulation param
trueModel = 'PDN';
trueTheta = [1, 0, 8, 15];
% trueModel = 'RemiProbitNorm';
% trueTheta = [1.11,1.06,0.88,0.70,0.85,0.64,0.62,2.28,19.47,1.43,1.82,1.17,1.14,0.49,0.55,0.93]; %subject 1
% trueTheta = [0.83,1.27,1.08,1.14,1.07,1.12,1.71,2.83,7.03,1.37,7.73,1.86,2.84,1.26,0.99,0.67]; %subject 2
% trueTheta = [0.84 1.01 0.97 0.59 0.59 0.72 0.42 1.76 20.76 2.16 1.22 2.04 0.68 0.59 0.51 0.95]; %subject 3

%% Particles Initialization
param = struct;
param.G = 2;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 3;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
Models = {'PDN','Logit'};%;'Logit'
Particles = cell(numel(Models),1);
initTheta = cell(param.G,param.P);
for m=1:numel(Models)
    Particles{m} = struct;
    Particles{m}.model = Models{m};
    Particles{m}.OptimTheta = coder.nullcopy(initTheta);
    Particles{m}.NewTheta = coder.nullcopy(initTheta);
    Particles{m}.param = param;
    Particles{m}.attrSign = attrSign;
    %tracking
    Particles{m}.log_w_bar =[];
    Particles{m}.c_phase_transitions = [];
    Particles{m}.LogMargLik_hist_group = zeros(param.G,opt_num_quest+num_consistency_check);
    Particles{m}.LogMargLik_hist = zeros(1,opt_num_quest+num_consistency_check);
    Particles{m}.log_marg_like = 0;
    for g=1:param.G
        for p=1:param.P
            Particles{m}.OptimTheta{g,p} = InitParticle(Particles{m}.model,param.K);
        end
    end
end
%% Phase 1
% X_Doptim = DOptimDesign( opt_num_quest,3 , attrVals );
Xs = cell(num_non_check+num_consistency_check,1);
XIndifs = cell(num_non_check+num_consistency_check,1);
TargetAndAltX = nan(num_non_check+num_consistency_check,2);
ChoiceList = nan(num_non_check+num_consistency_check,1);
ConsistencyCheck = false(num_consistency_check,1);
time = zeros(opt_num_quest+num_consistency_check+1,3);
timer = tic;
list_true_proba = zeros(opt_num_quest+num_consistency_check,1);
for obs = 1:num_non_check+num_consistency_check
    %% Optimal Design
    fprintf('\n==Iteration %.0f==\n',obs);
%     Xs{obs,1} = reshape( X_Doptim(obs,:),[param.K size(X_Doptim(obs,:),2)/param.K] )';
    if obs <= opt_num_quest
        %% Optimal Questions
        X_optim = OptimDesign( Particles, num_option_list(obs), attrVals, attrSign, 'Shannon' );
        Xs{obs,1} = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    elseif obs <= opt_num_quest + num_double_decoy
        %% Compute Decoys
        X_indif = FindIndif( Particles, attrVals, attrSign );
        X_indif = reshape( X_indif,[param.K size(X_indif,2)/param.K] )';
        XIndifs{obs,1} = X_indif;
        [ X_1decoy, X_2decoy , TargetAndAltX1, TargetAndAltX2 ] = AddDecoy(X_indif, attrVals, attrSign );
        Xs{obs,1} = X_1decoy; %Single decoy
        Xs{list_double_decoy(obs-opt_num_quest),1} = X_2decoy; %Double decoy
        TargetAndAltX(obs,:) = TargetAndAltX1;
        TargetAndAltX(list_double_decoy(obs-opt_num_quest),:) = TargetAndAltX2;
        J = size(Xs{obs,1},1);
        num_option_list = [num_option_list;J];
        %% debug
%         list_double_decoy(obs-opt_num_quest)
%         Xs{list_double_decoy(obs-opt_num_quest),1}
%         TargetAndAltX(list_double_decoy(obs-opt_num_quest),:)
%         TargetAndAltX2
    elseif obs <= opt_num_quest + 2*num_double_decoy
        J = size(Xs{obs,1},1);
        num_option_list = [num_option_list;J];
    else
        %% Consistency check
        redo_obs = list_const_check(obs-num_non_check);
        J = size(Xs{redo_obs,1},1);
        num_option_list = [num_option_list;J];
        option_permuted = randperm(J);
        Xs{obs,1} =  Xs{redo_obs,1}(option_permuted,:);
        expected_choice = sum((1:J) .* (ChoiceList(redo_obs)==option_permuted));
        TargetAndAltX(obs,1) = find(TargetAndAltX(redo_obs,1)==option_permuted);
        TargetAndAltX(obs,2) = find(TargetAndAltX(redo_obs,2)==option_permuted);
    end
    %% Observe
    [choice,list_true_proba(obs)] = SimChoice(Xs{obs},trueModel,trueTheta,attrSign, param);
    fprintf("Choice: %d\n",choice);
    ChoiceList(obs) = choice;
    time(obs+1,2) = toc(timer)-sum(time(:,1));
    % Consistency check
    if obs > num_non_check
        ConsistencyCheck(obs-opt_num_quest) = expected_choice == choice;
    end
    %% Update Particles
    [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
%     BF = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like);
    if numel(Models) > 2
        BF31 = exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like);
        BF32 = exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like);
        fprintf('BF 2/1: %.3g\nBF 3/1: %.3g\nBF 3/2: %.3g\n',BF,BF31,BF32);
    end
    time(obs+1,3) = toc(timer)-sum(time(:,1))-time(obs+1,2);
    time(obs+1,1) = toc(timer)-sum(time(:,1));
    [time(obs+1,:);mean(time(2:obs+1,:),1)]
    NumDecoychoice = sum(ChoiceList(1:obs) == TargetAndAltX(1:obs,:))
    %% Plot results
    LML = [];
    progLML = [];
    progBFvsUnif = [];
    for m=1:numel(Models)
        LML = [LML;real(Particles{m}.LogMargLik_hist(1:obs)) ];
        progLML = [progLML;cumsum(real(Particles{m}.LogMargLik_hist(1:obs)))./(1:obs) ];
        progBFvsUnif = [progBFvsUnif;cumsum(Particles{m}.LogMargLik_hist(1:obs))-cumsum(log(1./num_option_list(1:obs)))' ];
    end
    progBFvsUnifShift20 = [zeros(numel(Models),20) progBFvsUnif(:,1:max(0,obs-20))];
    progBFvsUnifLast20 = progBFvsUnif- progBFvsUnifShift20(:,1:obs) ;
    ML = [exp(LML);list_true_proba(1:obs)' ];
    subplot(2,2,1, 'YScale', 'log');
    plot(ML')
    title('ML')
    subplot(2,2,2, 'YScale', 'log');
    plot(exp(progLML)')
    title('average progLML')
    subplot(2,2,3, 'YScale', 'log');
    semilogy(exp(progBFvsUnif)')
    title('BF vs Unif')
    subplot(2,2,4, 'YScale', 'log');
    semilogy(exp(progBFvsUnifLast20)')
    title('BF vs Unif (Last 20)')
    drawnow;
end

%% Plot with posterior probas
postmarglikelihood = nan(numel(Xs),1);
list_ML = zeros(4,numel(Xs));
for obs = 1:numel(Xs)
    P = numel(Particles{1}.OptimTheta);
    likelihood = zeros(num_option_list(obs),P);
    for p=1:P
        likelihood(:,p) = ProbaChoice( Xs{obs} , Particles{1}.model, Particles{1}.OptimTheta{p}, attrSign, param );
    end
    marginal_likelihood = sum(likelihood,2)/P;
    list_ML(1:numel(marginal_likelihood),obs)=marginal_likelihood;
    postmarglikelihood(obs) = marginal_likelihood(ChoiceList(obs));
end
subplot(2,2,1, 'YScale', 'log');
plot([ML' postmarglikelihood]);
title('ML');

%% Prop Decoy and Alt chosen
list_decoy_quest = find(~isnan(TargetAndAltX(:,1)));
PropChoice = NumDecoychoice ./ numel(list_decoy_quest)
ExpectedPropChoice = [0 0];
TrueBinaryChoice = [0 0];
TruePropChoice = [0 0];
test = [];
proba_target =zeros(numel(list_decoy_quest),1);
proba_target_binary =zeros(numel(list_decoy_quest),1);
for i = 1:numel(list_decoy_quest)
    obs = list_decoy_quest(i);
    ExpectedPropChoice = ExpectedPropChoice +  list_ML(TargetAndAltX(obs,:),obs)';
    likelihood = ProbaChoice( Xs{obs} , Particles{1}.model, trueTheta, attrSign, param );
    likelihoodBinary = ProbaChoice( Xs{obs}(TargetAndAltX(obs,:),:) , Particles{1}.model, trueTheta, attrSign, param );
    proba_target(i) = likelihood(TargetAndAltX(obs,1)) / likelihood(TargetAndAltX(obs,2));
    proba_target_binary(i) = likelihoodBinary(1) / likelihoodBinary(2);
    TruePropChoice = TruePropChoice + likelihood(TargetAndAltX(obs,:))';
    TrueBinaryChoice = TrueBinaryChoice + likelihoodBinary';
end
subplot(2,2,4);
plot([proba_target proba_target_binary]);
title('P(T|Decoy)/P(Alt|Decoy) and P(T|Bin)/P(Alt|Bin)')
ExpectedPropChoice = ExpectedPropChoice / numel(list_decoy_quest)
TruePropChoice = TruePropChoice / numel(list_decoy_quest);
TrueBinaryChoice = TrueBinaryChoice / numel(list_decoy_quest);
Ratios = [PropChoice(1)/PropChoice(2) TruePropChoice(1)/TruePropChoice(2),TrueBinaryChoice(1)/TrueBinaryChoice(2)] 
%% Results
% logML = [Particles{2}.LogMargLik_hist;Particles{1}.LogMargLik_hist];
% proglnML = [exp(cumsum(Particles{2}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2))) ; exp(cumsum(Particles{1}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2)))];
% progRatio = proglnML(1,:) ./ proglnML(2,:);
% logBF_list = cumsum(Particles{2}.LogMargLik_hist) - cumsum(Particles{1}.LogMargLik_hist);
% plot(logBF_list)
