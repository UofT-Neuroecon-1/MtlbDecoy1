function [ X_optim ] = DOptimDesign( dopt_num_quest, num_opt, attrVals )
%DOPTIMDESIGN Returns D-Optimal num_quest questions, with num_opt choice
%options and using attrVals as an attribute space.
    %% get random start value
    K=size(attrVals,1);
    X0 = zeros(dopt_num_quest,num_opt*K);
    for opt =1:num_opt
        for k = 1:K
            X0(:,(opt-1)*K+k) = randsample(attrVals{k},dopt_num_quest,true);
        end
    end
    %%
    DObjectiveFunction = @(X) -det(X'*X);
    X_optim = PathToMin(DObjectiveFunction,X0, attrVals );    
end

