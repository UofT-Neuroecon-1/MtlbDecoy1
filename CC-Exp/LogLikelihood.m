function [ logLik ] = LogLikelihood( Xs, ChoiceList , model , Theta, attrSign, param )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
logLik = 0;
T = size(Xs,1);
parfor t = 1:T
   proba_choice = ProbaChoice( Xs{t} , model , Theta, attrSign, param );
   logLik = logLik + log(proba_choice(ChoiceList(t)));
end

end

