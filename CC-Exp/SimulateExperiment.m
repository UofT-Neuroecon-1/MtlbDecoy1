% clear all;
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
Debug = false;
global attrVals attrNames attrSign
attrVals = cell(7,1);
attrNames = cell(7,1);
attrVals{1}= (5:25);          attrNames{1}="Annual fee ($)";
attrVals{2}= (0:0.25:1.5);    attrNames{2}="Transaction fee ($)";
attrVals{3}= (5:1:25);        attrNames{3}="Interest charged (%)";
attrVals{4}= 0:0.25:2.5;    attrNames{4}="Reward (%)";
attrVals{5}= (0:0.25:2.5);    attrNames{5}="ATM surcharge ($)";
attrVals{6}= (0:0.25:2.5);    attrNames{6}="Currency conversion fee (%)";
attrVals{7}= [0,15,30,60,90,120];    attrNames{7}="Warranty on purchase (days)";
attrSign = [-1, -1, -1, 1, -1, -1, 1];
num_option_list = [2*ones(1,1);3*ones(2,1);4*ones(3,1)];
opt_num_quest = numel(num_option_list);

%% Simulation param
trueModel = 'PDN';
trueTheta = [0.8, 0.1, 4, 6, 4, 6, 4, 6, 0];
trueModel = 'RemiProbitNorm';
trueTheta = [0.83,1.27,1.08,1.14,1.07,1.12,1.71,2.83,7.03,1.37,7.73,1.86,2.84,1.26,0.99,0.67];

%% Particles Initialization
param = struct;
param.G = 2;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 1;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
Models = {'Logit';'PDN';'RemiProbitNorm'};
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
    Particles{m}.LogMargLik_hist_group = zeros(param.G,opt_num_quest);
    Particles{m}.LogMargLik_hist = zeros(1,opt_num_quest);
    Particles{m}.log_marg_like = 0;
    for g=1:param.G
        for p=1:param.P
            Particles{m}.OptimTheta{g,p} = InitParticle(Particles{m}.model,param.K);
        end
    end
end
%% Phase 1
% X_Doptim = DOptimDesign( opt_num_quest,3 , attrVals );
Xs = {};
ChoiceList = [];
time = zeros(opt_num_quest+1,3);
timer = tic;
list_true_proba = zeros(opt_num_quest,1);
for obs = 1:opt_num_quest
    %% Optimal Design
    fprintf('\n==Iteration %.0f==\n',obs);
%     Xs{obs,1} = reshape( X_Doptim(obs,:),[param.K size(X_Doptim(obs,:),2)/param.K] )';
    X_optim = OptimDesign( Particles, num_option_list(obs), attrVals, attrSign, 'ADO' );
    Xs{obs,1} = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    %% Observe
    [choice,list_true_proba(obs)] = SimChoice(Xs{obs,1},trueModel,trueTheta,attrSign, param);
    fprintf("Choice: %d\n",choice);
    ChoiceList = [ChoiceList; choice];
    time(obs+1,2) = toc(timer)-sum(time(:,1));
    %% Update Particles
    [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
    BF = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like);
    if numel(Models) > 2
        BF31 = exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like);
        BF32 = exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like);
        fprintf('BF 2/1: %.3g\nBF 3/1: %.3g\nBF 3/2: %.3g\n',BF,BF31,BF32);
    end
    time(obs+1,3) = toc(timer)-sum(time(:,1))-time(obs+1,2);
    time(obs+1,1) = toc(timer)-sum(time(:,1));
    [time(obs+1,:);mean(time(2:obs+1,:),1)]
    %%
    LML = [];
    progLML = [];
    progBFvsUnif = [];
    for m=1:numel(Models)
        LML = [LML;Particles{m}.LogMargLik_hist(1:obs) ];
        progLML = [progLML;cumsum(Particles{m}.LogMargLik_hist(1:obs))./(1:obs) ];
        progBFvsUnif = [progBFvsUnif;cumsum(Particles{m}.LogMargLik_hist(1:obs))-cumsum(log(1./num_option_list(1:obs)))' ];
    end
    ML = [exp(LML);list_true_proba(1:obs)' ];
    subplot(3,1,1, 'YScale', 'log');
    plot(ML')
    title('ML')
    subplot(3,1,2, 'YScale', 'log');
    plot(exp(progLML)')
    title('average progLML')
    subplot(3,1,3, 'YScale', 'log');
    semilogy(exp(progBFvsUnif)')
    title('BF vs Unif')
    drawnow;
end

%% Results
% logML = [Particles{2}.LogMargLik_hist;Particles{1}.LogMargLik_hist];
% proglnML = [exp(cumsum(Particles{2}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2))) ; exp(cumsum(Particles{1}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2)))];
% progRatio = proglnML(1,:) ./ proglnML(2,:);
% logBF_list = cumsum(Particles{2}.LogMargLik_hist) - cumsum(Particles{1}.LogMargLik_hist);
% plot(logBF_list)
