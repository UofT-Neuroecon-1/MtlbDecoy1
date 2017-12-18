function [ logLik ] = LogLikelihood( Xs, ChoiceList , model , Theta, attrSign )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
logLik = 0;
T = size(Xs,1);
for t = 1:T
   proba_choice = ProbaChoice( Xs{t} , model , Theta, attrSign );
   logLik = logLik + log(proba_choice(ChoiceList(t)));
end

end

