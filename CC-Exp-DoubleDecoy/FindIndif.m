function [ X ] = FindIndif( Particles, attrVals, attrSign )
%FINDINDIF Returns a set of 2 options close to indiference

% Get particles info
M = size(Particles,1);
K = size(attrVals,1);

% Init variables:
num_opt = 2;

%% Objective function
function obj_val = ObjectiveIndif( X, subTheta,param, model )
    subP = numel(subTheta);
    X = reshape( X,[param.K size(X,2)/param.K] )';
    likelihood = zeros(num_opt,subP);
    obj_val = 0 + 1000 * all(X(1,:)==X(2,:));
    for p=1:subP
        likelihood(:,p) =ProbaChoice( X , model, subTheta{p}, attrSign, param );
        obj_val = obj_val + (likelihood(1,p)-likelihood(2,p))^2;
    end
end


function marginal_likelihood = MarginalLikelihood( X, subTheta,param, model )
    subP = numel(subTheta);
    X = reshape( X,[param.K size(X,2)/param.K] )';
    likelihood = zeros(num_opt,subP);
    for p=1:subP
        likelihood(:,p) =ProbaChoice( X , model, subTheta{p}, attrSign, param );
    end
    marginal_likelihood = sum(likelihood,2)/subP;
end

%% Algo
OptimTheta = Particles{1}.OptimTheta;
param = Particles{1}.param;
subTheta = OptimTheta(:,1:floor(param.P/2));
objective = @(X) ObjectiveIndif( X, subTheta,param, Particles{1}.model );

%Draw num_start starting points
% num_start = 4;
% X0 = cell(num_start,1);
% target_val = zeros(num_start,1);
% for start=1:num_start
%     X0{start} = zeros(1,num_opt*K);
%     for opt =1:num_opt
%         for k = 1:K
%             X0{start}(1,K*(opt-1)+k) = randsample(attrVals{k},1,true);
%         end
%     end
%     target_val(start) = objective(X0{start});
% end
% best_start = find(target_val == min(target_val),1);
X_1 = zeros(1,2);
X_2 = zeros(1,2);
range1 = numel(attrVals{1});
range2 = numel(attrVals{2});
%find a valid indifference set
valid = false;
while ~valid
    % First attribute
    X_1(1) = randsample(attrVals{1}(floor(range1/8):floor(range1/2)-2),1);
    X_1(2) = randsample(attrVals{1}(ceil(range1/2)+2:ceil(7*range1/8)),1);
    % Second attribute
    X_2(1) = randsample(attrVals{2}(2:floor(range2/2)),1);
    X_2(2) = randsample(attrVals{2}(ceil(range2/2):end),1);

    objective = @(X_2) ObjectiveIndif( [X_1(1) X_2(1) X_1(2) X_2(2)], subTheta,param, Particles{1}.model );
    %% Optimize from best starting point
    try
        X_2 = PathToMin(objective,X_2, {attrVals{2}(2:end)} );
    catch ME
        sca();
        error('Error: Call the supervisor');
    end
    X = [X_1(1) X_2(1) X_1(2) X_2(2)];
    %% Check if they are at least 2 stepsizes away
    A1 = abs([X_1(2)-X_1(1) X_2(2)-X_2(1)]);
    A2 = 4*[attrVals{1}(2)-attrVals{1}(1) attrVals{2}(2)-attrVals{2}(1)];
    valid = all(A1 >= A2);
end


%% Check choice proba
marginal_likelihood = MarginalLikelihood( X, subTheta,param, Particles{1}.model )
if marginal_likelihood(1) < marginal_likelihood(2)
    X = [X_1(2) X_2(2) X_1(1) X_2(1)];
end


end

