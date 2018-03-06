%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation / Model-comparison Example File %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

addpath(['..' filesep 'EstimationAdaptive' filesep])

%% Load Backup
% if you want to load a backup from a previous estimation (set to [] if
% not)
backup_file = [];'backup1.mat';
%% Estimation Parameters
% You can create your own parameters that will be passed to the likelihood
% function and the particles initialization functions
param = struct;
param.G = 3; % Number of particles group
param.P = 128; % Number of particles per group
param.Adaptive = true; % Use the adaptive SMC (see Durham, Geweke 2014).
param.Msteps = 10; % Number of mutate steps
param.Tag = 'StndPDN'; % This tag will be added to the output file
param.ress_threshold = 0.8;

param.savefile = ['Analysis' filesep param.Tag sprintf('-%.0fx%.0f-M%.0f-',param.G,param.P,param.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];
if ~exist('Analysis','dir')
    mkdir('Analysis')
end
% Models to use for estimations
% Each model should have a corresponding entry in the following files:
% InitParticle : Returns a draw from prior for one particle
% Mutate : The Metropolis-Hastings Mutation step
% ProbaChoice : The likelihood of one observation given a model and particle
param.Models = {'PDNNew'}; %{'RemiStand';'HierarchicalProbit'};

% Model specific parameters
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000); % Draws for GHK
% param.attrVals = cell(7,1);
% param.attrNames = cell(7,1);
% param.attrVals{1}= (5:25);          param.attrNames{1}="Annual fee ($)";
% param.attrVals{2}= (0:0.1:0.5);    param.attrNames{2}="Transaction fee ($)";
% param.attrVals{3}= (5:1:25);        param.attrNames{3}="Interest charged (%)";
% param.attrVals{4}= 0:0.25:2.5;    param.attrNames{4}="Reward (%)";
% param.attrVals{5}= (0:0.25:2.5);    param.attrNames{5}="ATM surcharge ($)";
% param.attrVals{6}= (0:0.25:2.5);    param.attrNames{6}="Currency conversion fee (%)";
% param.attrVals{7}= [0,15,30,60,90,120];    param.attrNames{7}="Warranty on purchase (days)";
% param.attrSign = [-1, -1, -1, 1, -1, -1, 1];

param.attrVals{1}= (0:4);          param.attrNames{1}="Bid";
param.attrSign = [1];

param.K = numel(param.attrVals); 
param.attrMax = zeros(1,param.K);
for k=1:param.K
    param.attrMax(k) = max(param.attrVals{k});
end

clear k
%% Data format
% Requires a cell array (e.g. SubjData{}) of size N x 1
% Each cell must contain a structure with at least 2 fields:
% SubjData{n}.Xs : cell array of T observations.
% - SubjData{n}.Xs{t}: Matrix of J(t) options x K(t) attributes
% SubjData{n}.ChoiceList : Vector of T x 1 Choices

load('simulout.mat')
for s=1:40
    for t=1:250
        SubjData{s}.Xs{t}=V(t,:,s)';
    end
    SubjData{s}.ChoiceList=choice{s};
end
save ExampleData 'SubjData' 'par'

% load('simulout.mat')
% for s=1:40
%     for t=1:250
%         temp{t}=V(t,:,s)';
%     end
%     SubjData{s}.ChoiceList=choice{s};
% end
% SubjData{s}.Xs{t}
% save ExampleData 'SubjData' 'par'

% % Load files list
% fileslist = dir(['..\CC-Exp\Analysis\LabData' filesep 'Optim-BRLAB*.mat']);
% num_subj = size(fileslist,1);
% Data = cell(num_subj,1);
% subjList = [];
% for file = 1:num_subj
%     Data{file} = load([fileslist(file).folder filesep fileslist(file).name]);
%     % add files to subjList if they satisfy conditions (min # answers,
%     % consistency,...)
%     if numel(Data{file}.ChoiceList) > 50
%         if isfield(Data{file},'ConsistencyCheck')
%             if mean(Data{file}.ConsistencyCheck) > 0.3
%                 subjList = [subjList;file];
%             end
%         else
%             subjList = [subjList;file];
%         end
%     end
% end
% % Keep only valid subjects
% SubjData = Data(subjList);
% 
% clear Data file fileslist num_subj subjList

load ExampleData

%% Estimation
%use adaptive algorthm
% EstimationOutput = EstimationAdaptiveSMC( SubjData, param, backup_file)
% Particles = EstimationOutput.Particles;
% 
% if strcmp(param.Models{1},'PDNNew')
%     EstimationOutput.Particles{1}.postmeans
%     fprintf('Posterior Mean (Across Subjects) \n')
%     fprintf('alpha: %f \n sigma: %f \n omega: %f \n',mean(EstimationOutput.Particles{1}.postmeans))
% end


% % Full posterior (all the subjects superposed)
% prior = [betarnd(3,1,10000,1) gamrnd(1,0.5,10000,1) gamrnd(1,1,10000,1)];
% prior_mean = mean(prior,1);
% % Compute posterior means
% size_NK = size(Particles{1}.particle{1}.theta);
% VectorizedTheta = nan(param.P*param.G,size_NK(1),size_NK(2));
% for p=1:param.G*param.P
%     VectorizedTheta(p,:,:) = Particles{1}.particle{p}.theta;
% end
% post_mean = squeeze(mean(mean(VectorizedTheta,1),2,'omitnan'))
% % Plot posteriors
% for theta=1:3
%     subplot(2,3,theta)
%     histogram(prior(:,theta),20,'Normalization','pdf')
%     title(sprintf('prior mean: %.2f',prior_mean(theta)));
%     subplot(2,3,3+theta)
%     hist_data = VectorizedTheta(:,:,theta);
%     hist_data = hist_data(~isnan(hist_data(:)))
%     histogram(hist_data(:),20,'Normalization','pdf')
%     title(sprintf('post mean: %.2f',post_mean(theta)));
% end

%% Maximum Likelihood (within)
% theta0 = [1,1,1];
% theta_indiv_ml = nan(numel(SubjData),3);
% for subj = 1:numel(SubjData)
%     Target = @(theta) -LogLikelihood( SubjData{subj}.Xs, SubjData{subj}.ChoiceList, 1 , 'PDNProbit' , theta, param );
%     theta_indiv_ml(subj,:) = fminunc(Target,theta0);
% end
%% Maximum Likelihood (pooled)
% pool data
PooledXs = {};
PooledChoiceList = [];
for subj=1:numel(SubjData)
    PooledXs = [PooledXs,SubjData{subj}.Xs];
    PooledChoiceList = [PooledChoiceList;SubjData{subj}.ChoiceList];
end
theta0 = [1,0,1];
options = optimoptions('fminunc','OptimalityTolerance',1.0e-9,'Display','iter-detailed');
Target = @(theta) -LogLikelihood( PooledXs, PooledChoiceList, 1 , 'DNv' , theta, param );
theta_pooled = fminunc(Target,theta0,options)

%% Plot Likelihoods conditional on true values
true_theta = [par.a,par.sigma,par.omega]
gridpoints = linspace(0,3);
yy = nan(numel(gridpoints),3);
for th = 1:3
    theta = true_theta;
    for i = 1:numel(gridpoints)
        theta(th) = gridpoints(i);
        yy(i,th) = LogLikelihood( PooledXs, PooledChoiceList, 1 , 'DNv' , theta, param );
    end
    subplot(3,1,th);
    plot(gridpoints,yy(:,th));
    title(sprintf('LogLik param %d',th));
end
