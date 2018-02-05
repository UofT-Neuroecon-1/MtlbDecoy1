function [ X ] = OptimDesign( Particles, num_opt, attrVals, attrSign, type )
%OPTIMDESIGN Returns the BOSIE optimal design
% Particles: A cell containing particles group structures (One structure
% per model)
% num_opt: number of alternatives in the question
% attrVals: possible attribute values
% type: of optimal desighn: 'model_select' (ADO), 'Shannon' ...

% Get particles info
M = 2;size(Particles,1);
K = size(attrVals,1);

%% Objective functions

function obj_val = ObjectiveFunctionShannon( X, subTheta,param )
    subP = numel(subTheta);
    X = reshape( X,[param.K size(X,2)/param.K] )';
    likelihood = zeros(num_opt,subP);
    for p=1:subP
        likelihood(:,p) =ProbaChoice( X , param.model, subTheta{p}, attrSign, param );
    end
    marginal_likelihood = sum(likelihood,2)/subP;
    obj_val = 0;
    for j=1:num_opt
        for p=1:subP
            obj_val = obj_val + log(likelihood(j,p)/marginal_likelihood(j)) * likelihood(j,p);
        end
    end
end

function obj_val = ObjectiveFunctionADO( X, Particles)
    X = reshape( X,[K size(X,2)/K] )';
    u_m = zeros(M,1);
    marg_like = zeros(M,1);
    for m=1:M
        param = Particles{m}.param;
        marg_like(m) = exp(Particles{m}.log_marg_like);
        subTheta = Particles{m}.OptimTheta(:,1:floor(param.P/2));
        subP = numel(subTheta);
        likelihood = zeros(num_opt,subP);
        parfor p=1:subP
            likelihood(:,p) =ProbaChoice( X , Particles{m}.model , subTheta{p}, attrSign, param );
        end
        marginal_pred_likelihood = sum(likelihood,2)/subP;
        for j=1:num_opt
            for p=1:subP
                u_m(m) = u_m(m) + log(likelihood(j,p)/marginal_pred_likelihood(j)) * likelihood(j,p);
            end
        end
    end
    obj_val = sum(marg_like.*u_m);
end

function obj_val = ObjectiveFunctionModel( X, Particles )
    X = reshape( X,[K size(X,2)/K] )';
    u_m = zeros(M,1);
    marg_like = zeros(M,1);
    marginal_py = cell(M,1);
    likelihoods = cell(M,1);
    for m=1:M
        param = Particles{m}.param;
        marg_like(m) = exp(Particles{m}.log_marg_like);
        subTheta = Particles{m}.OptimTheta(:,1:floor(param.P/2));
        subP = numel(subTheta);
        likelihoods{m} = zeros(num_opt,subP);
        for p=1:subP
            likelihoods{m}(:,p) =ProbaChoice( X , Particles{m}.model , subTheta{p}, attrSign ,param );
        end
        marginal_py{m} = sum(likelihoods{m},2)/subP;
    end
    for m=1:M
        for j=1:num_opt
            for p=1:subP
                u_m(m) = u_m(m) + log(marginal_py{3}(j)/marginal_py{2}(j))^2 * likelihoods{m}(j,p);
            end
        end
    end
    obj_val = sum(marg_like.*u_m);
end



function obj_val = ObjectiveFunctionModelADO( X, Particles )
    X = reshape( X,[K size(X,2)/K] )';
    u_m = zeros(M,1);
    marg_like = zeros(M,1);
    model_marginal_pred_like = cell(M,1);
    likelihoods = cell(M,1);
    marg_pred_like = 0;
    subP = zeros(M,1);
    for m=1:M
        param = Particles{m}.param;
        marg_like(m) = exp(Particles{m}.log_marg_like);
        subTheta = Particles{m}.OptimTheta(:,1:floor(param.P/2));
        subP(m) = numel(subTheta);
        likelihoods{m} = zeros(num_opt,subP(m));
        for p=1:subP(m)
            likelihoods{m}(:,p) =ProbaChoice( X , Particles{m}.model , subTheta{p}, attrSign, param );
        end
        model_marginal_pred_like{m} = sum(likelihoods{m},2)/subP(m);
        marg_pred_like = marg_pred_like + model_marginal_pred_like{m} * marg_like(m);
    end
    marg_pred_like = marg_pred_like ./ sum(marg_like);
    for m=1:M
        for j=1:num_opt
            for p=1:subP(m)
                if likelihoods{m}(j,p) > 0
                    u_m(m) = u_m(m) + log(model_marginal_pred_like{m}(j)/marg_pred_like(j)) * likelihoods{m}(j,p);
                end
            end
        end
    end
    obj_val = sum(marg_like.*u_m);
end

%% Find optimal question
if strcmp(type,'Shannon')
    for m=1:M
       if strcmp(Particles{m}.model,'PDN')
           OptimTheta = Particles{m}.OptimTheta;
           param = Particles{m}.param;
       end
    end
    subTheta = OptimTheta(:,1:floor(param.P/2));
    objective = @(X) -ObjectiveFunctionShannon( X, subTheta,param );
elseif strcmp(type,'ADO')
    objective = @(X) -ObjectiveFunctionModelADO( X, Particles );
elseif strcmp(type,'BosieModel')
    objective = @(X) -ObjectiveFunctionModel( X, Particles );
end
%Draw 2 starting points
num_start = 4;
X0 = cell(num_start,1);
target_val = zeros(num_start,1);
for start=1:num_start
    X0{start} = zeros(1,num_opt*K);
    for opt =1:num_opt
        for k = 1:K
            X0{start}(1,K*(opt-1)+k) = randsample(attrVals{k},1,true);
        end
    end
    target_val(start) = objective(X0{start});
end
best_start = find(target_val == min(target_val),1);
X = PathToMin(objective,X0{best_start}, attrVals );

%     options = optimoptions(@fmincon,'OptimalityTolerance',1.0e-3,'DiffMinChange',0.25,'Display','iter'); %,'Display','off');
%     [X,fval,exitflag,output] = fmincon(objective, X0(1:A,:),[],[],[],[],param.Xmin+zeros(1,A*param.K),zeros(1,A*param.K)+param.Xmax,[],options);

end
