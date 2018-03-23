clear all;
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
Debug = false;
global attrVals attrNames attrSign
attrVals = cell(7,1);
attrNames = cell(7,1);
attrVals{1}= (5:25);          attrNames{1}="Annual fee ($)";
attrVals{2}= (0:0.1:0.5);    attrNames{2}="Transaction fee ($)";
attrVals{3}= (5:1:25);        attrNames{3}="Interest charged (%)";
attrVals{4}= 0:0.25:2.5;    attrNames{4}="Reward (%)";
attrVals{5}= (0:0.25:2.5);    attrNames{5}="ATM surcharge ($)";
attrVals{6}= (0:0.25:2.5);    attrNames{6}="Currency conversion fee (%)";
attrVals{7}= [0,15,30,60,90,120];    attrNames{7}="Warranty on purchase (days)";
attrSign = [-1, -1, -1, 1, -1, -1, 1];

attrMax = zeros(1,numel(attrVals));
for k=1:numel(attrVals)
    attrMax(k) = max(attrVals{k});
end

%% Param Initialization
param = struct;
param.G = 2;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 4;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
param.attrVals = attrVals;
param.attrSign = attrSign;
param.attrMax = attrMax;
Models = {'Logit';'RemiProbitStandardized'};%'RemiProbitStandardized'};
Tag = 'NewVsStandMu2';

%% Load files list
fileslist = dir(['Analysis' filesep 'LabData' filesep 'Optim-BRLAB*.mat']);
num_subj = size(fileslist,1);
Data = cell(num_subj,1);
MaxNumQuest = 0;
subjList = [];
for file = 1:num_subj
    Data{file} = load([fileslist(file).folder filesep fileslist(file).name]);
    MaxNumQuest = max(MaxNumQuest,numel(Data{file}.ChoiceList));
    if numel(Data{file}.ChoiceList) > 50
        if isfield(Data{file},'ConsistencyCheck')
            if mean(Data{file}.ConsistencyCheck) > 0.3
                subjList = [subjList;file];
            end
        else
            subjList = [subjList;file];
        end
    end
end
num_subj = numel(subjList);
savefilename = ['Analysis' filesep 'Aggreg' filesep 'Aggr-' Tag sprintf('-%.0fx%.0f-M%.0f-',param.G,param.P,param.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];

%% Load previous analysis
OldAggreg = struct;
listoldfiles = dir(['Analysis' filesep 'Aggreg' filesep 'Aggr-' Tag sprintf('-%.0fx%.0f-M%.0f-*',param.G,param.P,param.Msteps)]);
dtnum_index = 0;
if numel(listoldfiles) > 0
    for f = 1:numel(listoldfiles)
        if datenum(listoldfiles(f).date) > dtnum_index
            OldAggreg = load([listoldfiles(f).folder filesep listoldfiles(f).name]);
            dtnum_index = datenum(listoldfiles(f).date);
        end
    end
end
%% Estimation
SubjData = cell(num_subj,1);
SubjBF = ones(num_subj,1);
SubjLogBF_list = zeros(MaxNumQuest,num_subj);
SubjLML = zeros(num_subj,numel(Models));
SubjLML20 = zeros(num_subj,numel(Models));
SubjLML40 = zeros(num_subj,numel(Models));
for subj = 1:num_subj
    %% Match with former run
    matched_with = 0;
    run_analysis = true;
    if isfield(OldAggreg,'SubjData')
        for oldsubj = 1:numel(OldAggreg.SubjData)
            if ~isempty(OldAggreg.SubjData{oldsubj})
                % Compare ids in old and new list
                if strcmp(OldAggreg.SubjData{oldsubj}.id,fileslist(subjList(subj)).name)
                    matched_with = oldsubj;
                    % Check if LML is not -inf and recompute if so
                    if isfield(OldAggreg.SubjData{matched_with},'num_option_list')
                        if ~any(isinf(sum(OldAggreg.SubjData{matched_with}.LML,2)))
                            SubjData{subj} = OldAggreg.SubjData{matched_with};
                            SubjLML(subj,:) = sum(SubjData{subj}.LML,2)';
                            SubjLML20(subj,:) = sum(SubjData{subj}.LML(:,20:end),2)';
                            SubjLML40(subj,:) = sum(SubjData{subj}.LML(:,40:end),2)';
                            num_option_list = OldAggreg.SubjData{matched_with}.num_option_list;
                            Particles = OldAggreg.SubjData{matched_with}.Particles;
                            run_analysis = false;
                        end
                    end
                end
            end
        end
    end
    %% Run Analysis
    if run_analysis
        Xs = Data{subjList(subj)}.Xs;
        ChoiceList = Data{subjList(subj)}.ChoiceList;
        SubjData{subj} = struct;
        SubjData{subj}.id = fileslist(subjList(subj)).name;
        SubjData{subj}.ChoiceList = ChoiceList;
        SubjData{subj}.Xs = Data{subjList(subj)}.Xs;
        opt_num_quest = size(Data{subjList(subj)}.ChoiceList,1);
        if isfield(Data{subjList(subj)},'ConsistencyCheck')
            SubjData{subj}.ConsistencyCheck = Data{subjList(subj)}.ConsistencyCheck;
            SubjData{subj}.list_const_check = Data{subjList(subj)}.list_const_check;
            SubjData{subj}.Consistency = sum(Data{subjList(subj)}.ConsistencyCheck) / numel(Data{subjList(subj)}.ConsistencyCheck);
        end
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
        TimerAnalysis = tic;
        AnalysisTime = zeros(opt_num_quest,1);
        num_option_list = [];
        for obs = 1:opt_num_quest
            %% Optimal Design
            fprintf('\n== Subject %.0f - Iteration %.0f ==\n',subj,obs);
            num_option_list = [num_option_list;size(Xs{obs},1)];
            %% Update Particles
            [ Particles ] = UpdateParticles( Particles, Xs(1:obs), ChoiceList(1:obs));
            SubjBF(subj) = exp(Particles{2}.log_marg_like - Particles{1}.log_marg_like);
            SubjData{subj}.BF =SubjBF(subj);
            if numel(Models) > 2
               SubjData{subj}.BF32 =exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like); 
               SubjData{subj}.BF31 =exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like); 
            end
            %% vs RemiProbit
            AnalysisTime(obs) = toc(TimerAnalysis);
            fprintf('last / average analysis time: %.2f %.2f\n',AnalysisTime(obs),AnalysisTime(obs)/obs);
            subjobs = sprintf("%d-%d",subj,obs);
            %% results for plots
            SubjData{subj}.LML = [];
            SubjData{subj}.progLML = [];
            SubjData{subj}.progBFvsUnif = [];
            for m=1:numel(Models)
                SubjData{subj}.LML = [SubjData{subj}.LML;Particles{m}.LogMargLik_hist(1:obs) ];
                SubjData{subj}.progLML = [SubjData{subj}.progLML;cumsum(Particles{m}.LogMargLik_hist(1:obs))./(1:obs) ];
                SubjData{subj}.progBFvsUnif = [SubjData{subj}.progBFvsUnif;cumsum(Particles{m}.LogMargLik_hist(1:obs))-cumsum(log(1./num_option_list(1:obs)))' ];
            end
            SubjLML(subj,:) = sum(SubjData{subj}.LML,2)';
            SubjLML20(subj,:) = sum(SubjData{subj}.LML(:,20:obs),2)';
            SubjLML40(subj,:) = sum(SubjData{subj}.LML(:,40:obs),2)';
            LBF = sum(log(SubjBF));
            LML = sum(SubjLML);
            LML20 = sum(SubjLML20);
            LML40 = sum(SubjLML40);
            % Plotting
            subplot(1,5,mod(subj-1,5)+1, 'YScale', 'log');
            semilogy(exp(SubjData{subj}.progBFvsUnif'))
            title(['BF' num2str(subj)])
            fig=gcf;
            fig.Position = [500 200 1400 750];
            drawnow;
            disp([LML;LML20;LML40]);
            %% Save
            if mod(obs,10) == 0
                save([savefilename  '-tmp.mat'],'SubjData','SubjLML','subj','obs','Models','LBF','LML','LML20','LML40');
            end
        end
        SubjData{subj}.AnalysisTime = AnalysisTime;
        SubjData{subj}.Particles = Particles;
        SubjData{subj}.num_option_list = num_option_list;
    end %end of analysis
    LBF = sum(log(SubjBF));
    LML = sum(SubjLML);
    LML20 = sum(SubjLML20);
    LML40 = sum(SubjLML40);
    save(savefilename,'SubjData','SubjLML','subj','Models','LBF','LML','LML20','LML40');
    % Plotting
    subplot(1,5,mod(subj-1,5)+1, 'YScale', 'log');
    semilogy(exp(SubjData{subj}.progBFvsUnif'))
    title(['BF' num2str(subj)])
    fig=gcf;
    fig.Position = [500 200 1400 750];
    drawnow;

end
