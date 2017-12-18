clear all;
if(exist('data')==0)
   mkdir('data') 
end
%% start parallel pool
% poolobj = gcp;

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
num_option_list = [2*ones(15,1);3*ones(15,1)];
opt_num_quest = numel(num_option_list);

%% Particles Initialization
param = struct;
param.G = 4;
param.P = 64;
param.K = size(attrVals,1);
param.Msteps = 1;
Models = {'Logit';'PDN'};
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

%% Simulation param
trueModel = 'PDN';
trueTheta = [0.8, 0.1, 4, 6, 4, 6, 4, 6, 0];

%% Phase 1
X_Doptim = DOptimDesign( opt_num_quest,3 , attrVals );
Xs = {};
ChoiceList = [];
time = 0;
timer = tic;
for obs = 1:opt_num_quest
    %% Optimal Design
    fprintf('\n==Iteration %.0f==\n',obs);
%     Xs{obs,1} = reshape( X_Doptim(obs,:),[param.K size(X_Doptim(obs,:),2)/param.K] )';
    X_optim = OptimDesign( Particles, num_option_list(obs), attrVals, attrSign, 'ADO' );
    Xs{obs,1} = reshape( X_optim,[param.K size(X_optim,2)/param.K] )';
    %% Observe
    choice = SimChoice(Xs{obs,1},trueModel,trueTheta,attrSign);
    fprintf("Choice: %d\n",choice);
    ChoiceList = [ChoiceList; choice];
    %% Update Particles
    [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
    BF = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like)
    time = [time;toc(timer)-sum(time)];
end

%% Results
logML = [Particles{2}.LogMargLik_hist;Particles{1}.LogMargLik_hist];
proglnML = [exp(cumsum(Particles{2}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2))) ; exp(cumsum(Particles{1}.LogMargLik_hist)./(1:size(Particles{1}.LogMargLik_hist,2)))];
progRatio = proglnML(1,:) ./ proglnML(2,:);
logBF_list = cumsum(Particles{2}.LogMargLik_hist) - cumsum(Particles{1}.LogMargLik_hist);
plot(logBF_list)