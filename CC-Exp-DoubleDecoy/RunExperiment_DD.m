clear all;
if(exist('data')==0)
   mkdir('data') 
end
%% start parallel pool
poolobj = gcp;

%% Experiment Parameters
global attrVals attrNames attrSign
attrVals = cell(2,1);
attrNames = cell(2,1);
attrVals{1}= (5:50);          attrNames{1}="Annual fee ($)";
attrVals{2}= 0:0.25:3;    attrNames{2}="Reward (%)";
attrSign = [-1, 1];
num_option_list = [2*ones(10,1);3*ones(20,1)];
opt_num_quest = numel(num_option_list);
num_double_decoy = 25;
list_double_decoy = opt_num_quest + num_double_decoy + randsample(num_double_decoy,num_double_decoy);
num_non_check = opt_num_quest + 2*num_double_decoy;
num_consistency_check = 20;
list_const_check = datasample(opt_num_quest+1:num_non_check,num_consistency_check,'Replace',false);
Debug = false;

%% Particles Initialization
param = struct;
param.G = 3;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 4;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
Models = {'PDN'};
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


%% Subject info
%Prompt to enter subject's details, which will later be used in the save
%file.

prompt = {'age', 'gender'};
defaults = {'0', 'F'};
answer = inputdlg(prompt, 'ChoiceRT', 2, defaults);
[subage, gender] = deal(answer{:}); % all inputs are strings
subid = num2str(floor(rand*1000));
%create data folder if doesnt exist
if ~exist('data', 'dir')
  mkdir('data');
end

%% Initialise Experiment Screen
rand('state', sum(100*clock));
Screen('Preference', 'SkipSyncTests', 1);
KbName('UnifyKeyNames');
LeftKey=KbName('LeftArrow'); UpKey=KbName('UpArrow'); RightKey = KbName('RightArrow'); DownKey = KbName('DownArrow');
LeftNKey=KbName('4'); UpNKey=KbName('8'); RightNKey = KbName('6'); DownNKey = KbName('2');
spaceKey = KbName('space'); escKey = KbName('ESCAPE');


gray = [127 127 127 ]; white = [ 255 255 255]; black = [ 0 0 0];
bgcolor = white; textcolor = black;

%% Run Experiment
time = 0;
timer = tic;

Xs = cell(num_non_check+num_consistency_check,1);
XIndifs = cell(num_non_check+num_consistency_check,1);
TargetAndAltX = nan(num_non_check+num_consistency_check,2);
ChoiceList = nan(num_non_check+num_consistency_check,1);
ConsistencyCheck = false(num_consistency_check,1);

run('Multiattribute_DD_Task.m');

toc(timer)

%% Stop Parallel Pool
delete(poolobj);
