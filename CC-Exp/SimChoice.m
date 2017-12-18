function [ choice_num ] = SimChoice( X , model , trueTheta,attrSign )
% X: J x K matrix of choice set
% model: the true model to use
% (returns) choice_num : the alternative j chosen by simulation

J = size(X,1);
proba_choice = ProbaChoice( X , model , trueTheta,attrSign );
proba_choice = proba_choice ./ sum(proba_choice );
choice_num = mnrnd(1,proba_choice')*(1:J)';

end

