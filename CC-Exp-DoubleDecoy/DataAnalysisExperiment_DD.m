clear all;
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
Debug = false;
global attrVals attrNames attrSign
attrVals = cell(2,1);
attrNames = cell(2,1);
attrVals{1}= (5:50);          attrNames{1}="Annual fee ($)";
attrVals{2}= 0:0.25:3;    attrNames{2}="Reward (%)";
attrSign = [-1, 1];

attrMax = zeros(1,numel(attrVals));
for k=1:numel(attrVals)
    attrMax(k) = max(attrVals{k});
end

%% Param Initialization
param = struct;
param.G = 2;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 5;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
param.attrVals = attrVals;
param.attrSign = attrSign;
param.attrMax = attrMax;
Models = {'Logit';'PDNNew';'RemiStandardized';'RemiProbitStandardized'};% 'RemiProbit';'PDNUnitIndep'};
Tag = 'NewVsStandDD';

%% Load files list
fileslist = dir(['Analysis' filesep 'LabData' filesep 'DD-BRLAB*.mat']);
num_subj = size(fileslist,1);
Data = cell(num_subj,1);
MaxNumQuest = 0;
subjList = [];

for file = 1:num_subj
    Data{file} = load([fileslist(file).folder filesep fileslist(file).name]);
    MaxNumQuest = max(MaxNumQuest,numel(Data{file}.ChoiceList));
    if numel(Data{file}.ChoiceList) > 50
        bin_decoy_obs = ~isnan(Data{file}.TargetAndAltX(:,1));
        DecoyOrAlt = sum(mean(Data{file}.TargetAndAltX(bin_decoy_obs,:)== Data{file}.ChoiceList(bin_decoy_obs,:)));
        if mean(Data{file}.ConsistencyCheck) > 0.3 && DecoyOrAlt > 0.7
            subjList = [subjList;file];
        end
    end
end
num_subj = numel(subjList);
savefilename = ['Analysis' filesep 'Aggreg' filesep 'AggDD-' Tag sprintf('-%.0fx%.0f-M%.0f-',param.G,param.P,param.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];


%% Decoy analysis
DecoyModelFree = struct;
DecoyModelFree.NumDecoychoiceList = nan(0,2);
Consistency = nan(num_subj,1);
AggregTargetChoice = nan(0,3);
AggregNumOptions = [];
AggregSubj = [];
single_prop_subj = [];
double_prop_subj = [];
num_sdec_questions = [];
num_ddec_questions = [];
mediantime = [];
for subj = 1:num_subj
    sub_data = Data{subjList(subj)};
    Consistency(subj) = mean(sub_data.ConsistencyCheck);
    bin_decoy_obs = ~isnan(sub_data.TargetAndAltX(:,1));
    DecoyOrAlt = sum(mean(sub_data.TargetAndAltX(bin_decoy_obs,:)== sub_data.ChoiceList(bin_decoy_obs,:)));
    median_rt = median(sub_data.timeRecords.answer-sub_data.timeRecords.show);
    if Consistency(subj) > 0.3 && DecoyOrAlt > 0.7 && median_rt > 2
        mediantime = [mediantime; median_rt];
        DecoyModelFree.NumDecoychoiceList = [DecoyModelFree.NumDecoychoiceList;sub_data.NumDecoychoice];
        for obs =1:numel(sub_data.ChoiceList)
            if bin_decoy_obs(obs)
                AggregNumOptions = [AggregNumOptions; size(sub_data.Xs{obs},1)];
                target_choice = sub_data.TargetAndAltX(obs,:)== sub_data.ChoiceList(obs);
                AggregTargetChoice = [AggregTargetChoice; target_choice 1-sum(target_choice)];
                AggregSubj = [AggregSubj;subj];
            end
        end
        num_sdec_questions = [num_sdec_questions;size(AggregTargetChoice(AggregNumOptions==3 & AggregSubj == subj,:),1)];
        num_ddec_questions = [num_ddec_questions;size(AggregTargetChoice(AggregNumOptions==4 & AggregSubj == subj,:),1)];
        single_prop_subj = [single_prop_subj; mean(AggregTargetChoice(AggregNumOptions==3 & AggregSubj == subj,:))];
        double_prop_subj = [double_prop_subj; mean(AggregTargetChoice(AggregNumOptions==4 & AggregSubj == subj,:))];
    end
end
DecoyModelFree.TotalDecoyChoice = sum(DecoyModelFree.NumDecoychoiceList);
propDecoy = DecoyModelFree.TotalDecoyChoice ./ sum(DecoyModelFree.TotalDecoyChoice)
DecoyModelFree.single_prop_subj = single_prop_subj;
DecoyModelFree.double_prop_subj = double_prop_subj;
DecoyModelFree.single_prop = mean(single_prop_subj(:,1) ./ sum(single_prop_subj(:,1:2),2))
DecoyModelFree.double_prop = mean(double_prop_subj(:,1) ./ sum(double_prop_subj(:,1:2),2))
DecoyModelFree.single_prop_se = std(single_prop_subj(:,1) ./ sum(single_prop_subj(:,1:2),2)) / sqrt(size(single_prop_subj(:,1),1));
DecoyModelFree.double_prop_se = std(double_prop_subj(:,1) ./ sum(double_prop_subj(:,1:2),2)) / sqrt(size(double_prop_subj(:,1),1));
DecoyModelFree.num_sdec_questions = num_sdec_questions;
DecoyModelFree.num_ddec_questions = num_ddec_questions;

%Double decoy prop | single decoy prop = true
list_has_sdec = single_prop_subj(:,1) > single_prop_subj(:,2);
DecoyModelFree.single_prop_subj_c = single_prop_subj(list_has_sdec,:);
DecoyModelFree.double_prop_subj_c = double_prop_subj(list_has_sdec,:);
prop_cond_sdec = mean(single_prop_subj(list_has_sdec,:));
prop_cond_ddec = mean(double_prop_subj(list_has_sdec,:));
DecoyModelFree.singleProp_cond = prop_cond_sdec(1:2) / sum(prop_cond_sdec(1:2));
DecoyModelFree.doubleProp_cond = prop_cond_ddec(1:2) / sum(prop_cond_ddec(1:2));
DecoyModelFree.singleProp_cond_se = std(single_prop_subj(list_has_sdec,1) ./ sum(single_prop_subj(list_has_sdec,1:2),2)) / sqrt(sum(list_has_sdec));
DecoyModelFree.doubleProp_cond_se = std(double_prop_subj(list_has_sdec,1) ./ sum(double_prop_subj(list_has_sdec,1:2),2))  / sqrt(sum(list_has_sdec));

%sort by median response time
[sorted_med_rtime,sorted_med_rtime_i] = sort(mediantime);

ratio_ptarget_single_double = (single_prop_subj(:,1)./sum(single_prop_subj(:,1:2),2)) ./ ...
    (double_prop_subj(:,1)./sum(double_prop_subj(:,1:2),2));
plot(repmat(sorted_med_rtime,1,2), [ratio_ptarget_single_double(sorted_med_rtime_i) ones(size(single_prop_subj,1),1)])
hold on
plot(repmat(sorted_med_rtime,1,2), [single_prop_subj(sorted_med_rtime_i) 0.5*ones(size(single_prop_subj,1),1)])
mean_ratio = mean(ratio_ptarget_single_double);
se_ratio = std(ratio_ptarget_single_double)/numel(ratio_ptarget_single_double)^0.5;
plot(sorted_med_rtime,mean(ratio_ptarget_single_double) * ones(size(single_prop_subj,1),1));
xlabel('median rt')
title({'P(T|S)/P(A|S) / P(T|D)/P(A|D)  per indiv (sorted by median RT)',sprintf("mean: %.4f , se: %.4f",mean_ratio,se_ratio)})
hold off

%% p values for P(Compet)
list_proba =[DecoyModelFree.single_prop_subj(:,2)./sum(DecoyModelFree.single_prop_subj(:,1:2),2) ...
    DecoyModelFree.double_prop_subj(:,2)./sum(DecoyModelFree.double_prop_subj(:,1:2),2)];
fprintf("\n\nProportions of Comp Choice (absolute)\n");
%% p values for P(Compet) | decoy effect
list_proba =[DecoyModelFree.single_prop_subj(list_has_sdec,2)./sum(DecoyModelFree.single_prop_subj(list_has_sdec,1:2),2) ...
    DecoyModelFree.double_prop_subj(list_has_sdec,2)./sum(DecoyModelFree.double_prop_subj(list_has_sdec,1:2),2)];
fprintf("\n\nProportions of Comp Choice (condit)\n");
%% Compute p-val
N = size(list_proba,1);
means = mean(list_proba);
std_err = std(list_proba,[],1) / sqrt(N);
tstat =  (means - 0.5)./ std_err;
p_val = tcdf(tstat,size(single_prop_subj,1)-1);
fprintf("n = %d\n",N);
fprintf("\t p_single \t p_double \nmean: \t%.4f \t%.4f\n",means(1),means(2));
fprintf("s.e.: \t%.4f \t%.4f\n",std_err(1),std_err(2));
fprintf("t(%d): %.3f\t%.3f\n",N-1,tstat(1),tstat(2));
fprintf("p-val: \t%.4f \t%.4f\n---\n",p_val(1),p_val(2));

% paired diff
diff_list = list_proba(:,1) - list_proba(:,2);
mean_diff=mean(diff_list);
sd_diff=std(diff_list);
tstat= sqrt( N ) * mean_diff / sd_diff;
p_val_diff = tcdf(tstat,N-1);
fprintf("mean diff: \t\t%.5f\n",mean_diff)
fprintf("s.e.(paired): \t%.5f\n",sd_diff/sqrt(size(single_prop_subj,1)))
fprintf("t (paired): \t%.5f\n",tstat)
fprintf("p (paired): \t%.5f\n",p_val_diff)


%% p values (Absolute)
tstat = ([DecoyModelFree.single_prop DecoyModelFree.double_prop]- 0.5)./ [DecoyModelFree.single_prop_se DecoyModelFree.double_prop_se];
p_val = tcdf(tstat,size(single_prop_subj,1)-1,'upper');
fprintf("\nAbsolute proportions of Target Choice\n");
fprintf("n = %d\n",size(single_prop_subj,1));
fprintf("\t p_single \t p_double \nmean: \t%.4f \t%.4f\n",DecoyModelFree.single_prop,DecoyModelFree.double_prop);
fprintf("s.e.: \t%.4f \t%.4f\n",DecoyModelFree.single_prop_se,DecoyModelFree.double_prop_se);
fprintf("t(%d): \t%.4f \t%.4f\n",size(single_prop_subj,1)-1,tstat(1),tstat(2));
fprintf("p-val: \t%.4f \t%.4f\n---\n",p_val(1),p_val(2));

% paired diff
mean_diff=mean(single_prop_subj(:,1) ./ sum(single_prop_subj(:,1:2),2) - double_prop_subj(:,1) ./ sum(double_prop_subj(:,1:2),2));
sd_diff=std(single_prop_subj(:,1) ./ sum(single_prop_subj(:,1:2),2) - double_prop_subj(:,1) ./ sum(double_prop_subj(:,1:2),2));
tstat= sqrt( size(single_prop_subj,1) ) * mean_diff / sd_diff;
p_val_diff = tcdf(tstat,size(single_prop_subj,1)-1,'upper');
fprintf("mean diff: \t\t%.5f\n",mean_diff)
fprintf("s.e.(paired): \t%.5f\n",sd_diff/sqrt(size(single_prop_subj,1)))
fprintf("t (paired): \t%.5f\n",tstat)
fprintf("p (paired): \t%.5f\n",p_val_diff)
%% p values (Conditional)
tstat = ([DecoyModelFree.singleProp_cond(:,1) DecoyModelFree.doubleProp_cond(:,1)]- 0.5)./ [DecoyModelFree.singleProp_cond_se DecoyModelFree.doubleProp_cond_se];
p_val = tcdf(tstat,size(single_prop_subj,1)-1,'upper');
fprintf("\n\nProportions of Target Choice Conditional on P(T|Single)>0.5\n");
fprintf("n = %d\n",sum(list_has_sdec,1));
fprintf("\t p_single \t p_double \nmean: \t%.4f \t%.4f\n",DecoyModelFree.singleProp_cond(:,1),DecoyModelFree.doubleProp_cond(:,1));
fprintf("s.e.: \t%.4f \t%.4f\n",DecoyModelFree.singleProp_cond_se,DecoyModelFree.doubleProp_cond_se);
fprintf("t(%d): \t%.4f \t%.4f\n",sum(list_has_sdec,1)-1,tstat(1),tstat(2));
fprintf("p-val: \t%.4f \t%.4f\n---\n",p_val(1),p_val(2));

% paired diff
diff_list = DecoyModelFree.single_prop_subj_c(:,1) ./ sum(DecoyModelFree.single_prop_subj_c(:,1:2),2) ...
    - DecoyModelFree.double_prop_subj_c(:,1) ./ sum(DecoyModelFree.double_prop_subj_c(:,1:2),2);
mean_diff=mean(diff_list);
sd_diff=std(diff_list);
tstat= sqrt( size(diff_list,1) ) * mean_diff / sd_diff;
p_val_diff = tcdf(tstat,size(diff_list,1)-1,'upper');
fprintf("mean diff: \t\t%.5f\n",mean_diff)
fprintf("s.e.(paired): \t%.5f\n",sd_diff/sqrt(size(single_prop_subj,1)))
fprintf("t (paired): \t%.5f\n",tstat)
fprintf("p (paired): \t%.5f\n",p_val_diff)

%% p values (Conditional on No effect)
n_cond_no = sum(~list_has_sdec,1);
single_prop_cond_no = single_prop_subj(~list_has_sdec,:);
double_prop_cond_no = double_prop_subj(~list_has_sdec,:);
se_cond_no = std([single_prop_cond_no(:,1),double_prop_cond_no(:,1)]) / sqrt(n_cond_no);
tstat =([mean(single_prop_cond_no(:,1)) mean(double_prop_cond_no(:,1))]- 0.5)./ se_cond_no;
p_val = tcdf(tstat,n_cond_no-1,'upper');
fprintf("\n\nProportions of Target Choice Conditional on P(T|Single)<=0.5\n");
fprintf("n = %d\n",n_cond_no);
fprintf("\t p_single \t p_double \nmean: \t%.4f\t%.4f\n",mean(single_prop_cond_no(:,1)),mean(double_prop_cond_no(:,1)));
fprintf("s.e.: \t%.4f \t%.4f\n", se_cond_no(1),se_cond_no(2));
fprintf("t(%d): \t%.3f\t%.3f\n",n_cond_no-1,tstat(1),tstat(2));
fprintf("p-val: \t%.4f \t%.4f\n---\n",p_val(1),p_val(2));

% paired diff
diff_list = single_prop_cond_no(:,1) ./ sum(single_prop_cond_no(:,1:2),2) ...
    - double_prop_cond_no(:,1) ./ sum(double_prop_cond_no(:,1:2),2);
mean_diff=mean(diff_list);
sd_diff=std(diff_list);
tstat= sqrt( size(diff_list,1) ) * mean_diff / sd_diff;
p_val_diff = tcdf(tstat,size(diff_list,1)-1,'upper');
fprintf("mean diff: \t\t%.5f\n",mean_diff)
fprintf("s.e.(paired): \t%.5f\n",sd_diff/sqrt(size(single_prop_subj,1)))
fprintf("t (paired): \t%.5f\n",tstat)
fprintf("p (paired): \t%.5f\n",p_val_diff)

%%
histogram(single_prop_subj(:,1)./sum(single_prop_subj(:,1:2),2),0.2:0.05:0.8)
hold on
histogram(double_prop_subj(:,1)./sum(double_prop_subj(:,1:2),2),0.2:0.05:0.8)
title('P(T|S)/P(A|S) (blue) and P(T|D)/P(A|D) (red)  per indiv ')
hold off

%%
subplot(2,1,1)
histogram(single_prop_subj(:,1)./sum(single_prop_subj(:,1:2),2),0.2:0.05:0.8)
subplot(2,1,2)
histogram(double_prop_subj(:,1)./sum(double_prop_subj(:,1:2),2),0.2:0.05:0.8,'FaceColor','red')
%% 
means = [DecoyModelFree.singleProp_cond(1),DecoyModelFree.doubleProp_cond(1)]
serr = [DecoyModelFree.singleProp_cond_se ,DecoyModelFree.doubleProp_cond_se ]
bar(means)
hold on
errorbar(1:2,means,serr,'.')
plot([0 3],[0.5 0.5])
hold off
title('P(target) / (P(target)+P(Alt)) (Conditional)')
set(gca,'xticklabel',{'Single Decoy','Double Decoy'})
%% 
means = [DecoyModelFree.single_prop(1),DecoyModelFree.double_prop(1)]
serr = [DecoyModelFree.single_prop_se ,DecoyModelFree.double_prop_se ]
bar(means)
hold on
errorbar(1:2,means,serr,'.')
plot([0 3],[0.5 0.5])
hold off
title('P(target) / (P(target)+P(Alt))')
set(gca,'xticklabel',{'Single Decoy','Double Decoy'})

%% Load previous analysis
OldAggreg = struct;
listoldfiles = dir(['Analysis' filesep 'Aggreg' filesep 'AggDD-' Tag sprintf('-%.0fx%.0f-M%.0f-*',param.G,param.P,param.Msteps)]);
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
                            SubjBF(subj) = OldAggreg.SubjBF(matched_with);
                            SubjLML(subj,:) = sum(SubjData{subj}.LML,2)';
                            SubjLogBF_list(1:size(OldAggreg.SubjLogBF_list,1),subj) = OldAggreg.SubjLogBF_list(:,matched_with);
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
            [exp(Particles{3}.log_marg_like - Particles{2}.log_marg_like) ; ...
                exp(Particles{3}.log_marg_like - Particles{1}.log_marg_like)]
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
            % Plotting
            subplot(1,5,mod(subj-1,5)+1, 'YScale', 'log');
            semilogy(exp(SubjData{subj}.progBFvsUnif'))
            title(['BF' num2str(subj)])
            fig=gcf;
            fig.Position = [200 200 1450 750];
            drawnow;
            %% Save
            LBF = sum(log(SubjBF));
            LML = sum(SubjLML);
            save([savefilename  '-tmp.mat'],'SubjData','SubjBF','SubjLogBF_list','SubjLML','subj','subjobs','Models','LBF','LML','DecoyModelFree');
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
        SubjLogBF_list(1:obs,subj) = SubjData{subj}.logBF_list';
        SubjData{subj}.num_option_list = num_option_list;
    end %end of analysis
    LBF = sum(log(SubjBF));
    LML = sum(SubjLML);
    save(savefilename,'SubjData','SubjBF','SubjLogBF_list','SubjLML','subj','Models','LBF','LML','DecoyModelFree');
    % Plotting
    subplot(1,5,mod(subj-1,5)+1, 'YScale', 'log');
    semilogy(exp(SubjData{subj}.progBFvsUnif'))
    title(['BF' num2str(subj)])
    fig=gcf;
    fig.Position = [200 200 1450 750];
    drawnow;

end
