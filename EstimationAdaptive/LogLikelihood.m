function [ logLik ] = LogLikelihood( Xs, ChoiceList, subj , model , theta, param )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
logLik = 0;
T = numel(Xs);
<<<<<<< HEAD
particle.theta=theta;
% for t = 1:T
%    proba_choice = ProbaChoice( Xs{t}, ChoiceList(t), subj , model , particle, param );
%    logLik = logLik + log(proba_choice(ChoiceList(t))); %this selects the prob for the choice only, so why calculate all of them?
% end

   P= ProbaChoice( Xs{t}, ChoiceList(t), subj , model , particle, param );
   logLik = sum(log(P))
=======

proba_c = @(X) ProbaChoice( X, subj , model , particle, param );
probas = cellfun(proba_c,Xs,'UniformOutput',0);
for t = 1:T
   proba_choice = probas{t};
   logLik = logLik + log(proba_choice(ChoiceList(t)));
end
% for t = 1:T
%    proba_choice = ProbaChoice( Xs{t}, subj , model , particle, param );
%    logLik = logLik + log(proba_choice(ChoiceList(t)));
% end
>>>>>>> origin/master

end

