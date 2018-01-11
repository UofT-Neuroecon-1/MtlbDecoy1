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
num_option_list = [2*ones(10,1);3*ones(20,1);4*ones(30,1)];
opt_num_quest = numel(num_option_list);

%% Particles Initialization
param = struct;
param.G = 4;
param.P = 128;
param.K = size(attrVals,1);
param.Msteps = 2;
param.NormDraw = mvnrnd(zeros(4,1),eye(4),1000);
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
spaceKey = KbName('space'); escKey = KbName('ESCAPE');
gray = [127 127 127 ]; white = [ 255 255 255]; black = [ 0 0 0];
bgcolor = white; textcolor = black;

%% Run Experiment
% X_Doptim = DOptimDesign( opt_num_quest,3 , attrVals );
Xs = {};
ChoiceList = [];
time = 0;
timer = tic;

run('Multiattribute_Ternary_Task.m');

toc(timer)

%% Stop Parallel Pool
delete(poolobj);
