clear all;
if(exist('data')==0)
   mkdir('data') 
end
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

%% Param Initialization
param = struct;
param.G = 2;
param.P = 256;
param.K = size(attrVals,1);
param.Msteps = 6;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
Models = {'Logit';'PDN';'RemiProbitNorm'};

%% Load files list
fileslist = dir(['data' filesep 'Opt*.mat']);
num_subj = size(fileslist,1);
Data = cell(num_subj,1);
for file = 1:num_subj
    Data{file} = load([fileslist(file).folder filesep fileslist(file).name]);
end
savefilename = ['data' filesep 'Aggreg-' sprintf('%.0fx%.0f-M%.0f-',param.G,param.P,param.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];

%% Estimation
SubjData = cell(num_subj,1);
SubjBF = zeros(num_subj,1);
SubjLogBF_list = zeros(40,num_subj);
for subj = 1:num_subj
    SubjData{subj} = struct;
    opt_num_quest = size(Data{subj}.ChoiceList,1);
    % Particle initialization
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
        Particles{m}.postmeans = [];
        for g=1:param.G
            for p=1:param.P
                Particles{m}.OptimTheta{g,p} = InitParticle(Particles{m}.model,param.K);
            end
        end
    end
    % Inference
    Xs = Data{subj}.Xs;
    ChoiceList = Data{subj}.ChoiceList;
    TimerAnalysis = tic;
    AnalysisTime = zeros(opt_num_quest,1);
    for obs = 1:opt_num_quest
        %% Optimal Design
        fprintf('\n== Subject %.0f - Iteration %.0f ==\n',subj,obs);
        %% Update Particles
        [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
        SubjBF(subj) = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like)
        SubjData{subj}.BF =SubjBF(subj);
        if numel(Models) > 2
           SubjData{subj}.BF32 =exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like); 
           SubjData{subj}.BF31 =exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like); 
        end
        %% vs RemiProbit
        [exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like) ; ...
            exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like)]
        AnalysisTime(obs) = toc(TimerAnalysis);
        fprintf('last / average analysis time: %.2f %.2f\n',AnalysisTime(obs),AnalysisTime(obs)/obs);
    end
    SubjData{subj}.AnalysisTime = AnalysisTime;
    SubjData{subj}.Particles = Particles;
    SubjData{subj}.logML = [];
    SubjData{subj}.progAvgML = [];
    for m=1:numel(Models)
        SubjData{subj}.logML = [SubjData{subj}.logML; Particles{m}.LogMargLik_hist];
        SubjData{subj}.progAvgML = [SubjData{subj}.progAvgML ; exp(cumsum(Particles{m}.LogMargLik_hist)./(1:numel(Particles{m}.LogMargLik_hist)))];
    end
    SubjData{subj}.progRatio21 = SubjData{subj}.progAvgML(2,:) ./ SubjData{subj}.progAvgML(1,:);
    SubjData{subj}.progRatio32 = SubjData{subj}.progAvgML(3,:) ./ SubjData{subj}.progAvgML(2,:);
    SubjData{subj}.progRatio31 = SubjData{subj}.progAvgML(3,:) ./ SubjData{subj}.progAvgML(1,:);
    SubjData{subj}.logBF_list = cumsum(Particles{2}.LogMargLik_hist) - cumsum(Particles{1}.LogMargLik_hist);
    SubjLogBF_list(:,subj) = SubjData{subj}.logBF_list';
    save(savefilename,'SubjData','SubjBF','SubjLogBF_list');
end

%% Misc

for subj = 1:num_subj
    LML = [];
    for m=1:numel(Models)
        LML = [LML;SubjData{subj}.Particles{m}.LogMargLik_hist ];
    end
    subplot(2,num_subj,subj, 'YScale', 'log');
    plot(exp(LML)')
    subplot(2,num_subj,3+subj, 'YScale', 'log');
    plot(LML')
end
waitforbuttonpress;
%%
semilogy(exp(SubjLogBF_list))
title('Progression of BF');
subplot(2,3,2, 'XScale', 'log', 'YScale', 'log')
plot([SubjData{3}.Particles{1}.LogMargLik_hist;SubjData{3}.Particles{2}.LogMargLik_hist;SubjData{3}.Particles{3}.LogMargLik_hist]')
