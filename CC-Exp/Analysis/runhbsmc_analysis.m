clear all;
%%
addpath('D:\Dropbox\Dev\DiscreteChoiceSMC\HBCSMCMatlab');
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
attrVals = cell(7,1);
attrNames = cell(7,1);
attrVals{1}= (5:25);          attrNames{1}="YrFee";
attrVals{2}= (0:0.1:0.5);    attrNames{2}="TransacFee";
attrVals{3}= (5:1:25);        attrNames{3}="Interest";
attrVals{4}= 0:0.25:2.5;    attrNames{4}="Reward";
attrVals{5}= (0:0.25:2.5);    attrNames{5}="ATMFee";
attrVals{6}= (0:0.25:2.5);    attrNames{6}="CurrExFee";
attrVals{7}= [0,15,30,60,90,120];    attrNames{7}="Warranty";
attrSign = [-1, -1, -1, 1, -1, -1, 1];

attrMax = zeros(1,numel(attrVals));
for k=1:numel(attrVals)
    attrMax(k) = max(attrVals{k});
end

%% Load files list
fileslist = dir(['LabData' filesep 'Optim-BRLAB*.mat']);
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



%% Param Initialization
param = struct;
param.G = 2;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 5;
param.ress_threshold = 0.8;
param.num_clust= 2;
param.Adaptive = true;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
param.attrVals = attrVals;
param.attrSign = attrSign;
param.attrMax = attrMax;
param.Models = {'HBC-PNE'};%'RemiProbitStandardized'};%'RemiProbitStandardized'};
param.Tag = 'PNE_beta';
param.basepath = 'D:\Data\Lab_Analysis\Out';
param.savefile = ['Aggr-' param.Tag sprintf('-%.0fx%.0f-M%.0f-',param.G,param.P,param.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];
%% Build SubjData List
SubjData = {};
for ss = 1:num_subj
    SubjData{ss}.Xs = Data{subjList(ss)}.Xs;
    SubjData{ss}.Ys = Data{subjList(ss)}.ChoiceList;
end


%% Estimation
backup_file = '';
EstimationOutput = HBCSMC( SubjData, param, backup_file );

%% Vectorize param
param = EstimationOutput.param;
N = numel(EstimationOutput.Particles{1}.particle{1,1}.theta);
vect_omega = nan(N,param.K,param.G,param.P);
vect_sig = nan(N,param.K,param.G,param.P);
for g = 1:param.G
    for p = 1:param.P
        for ss = 1:N
            vect_omega(ss,:,g,p)=EstimationOutput.Particles{1}.particle{g,p}.theta(ss).omega;
            vect_sig(ss,:,g,p)=EstimationOutput.Particles{1}.particle{g,p}.theta(ss).sig;
        end
    end
end
%% omega plot
for k=1:param.K
    subplot(3,param.K,k);
    histogram(vect_omega(:,k,:,:),0:0.5:15)
    title([attrNames{k} 'omega']);
    subplot(3,param.K,param.K+k);
    histogram(vect_sig(:,k,:,:),0:0.5:15)
    title('beta');
    subplot(3,param.K,2*param.K+k);
    histogram(vect_sig(:,k,:,:)./vect_omega(:,k,:,:),0:0.25:10)
    title('beta/omega');
end
%% beta plot
vect_beta = 1/vect_sig;
vect_om2 = vect_omega ./ vect_sig;
for k=1:param.K
    vect_beta_k = vect_beta(:,k,:,:);
    vect_om2_k = vect_om2(:,k,:,:);
    vect_beta_k = vect_beta_k(:);
    vect_om2_k = vect_om2_k(:);
    subplot(2,param.K,k);
    histogram(vect_beta_k,0:0.25:6,'Normalization','pdf')
    title([attrNames{k} 'beta']);
    subplot(2,param.K,param.K+k);
    histogram(vect_om2_k,0:0.5:8,'Normalization','pdf')
    title('omega');
end