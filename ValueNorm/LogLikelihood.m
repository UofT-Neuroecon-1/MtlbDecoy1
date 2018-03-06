function [ logLik ] = LogLikelihood( Xs, ChoiceList, subj , model , theta, param )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
logLik = 0;
T = numel(Xs);
Jt=cellfun(@length,Xs);

particle.theta=theta;

for t=1:T
    temp=eye(Jt(t)-1); 
    for i=1:Jt(t)
        M{i}=[temp(:,1:i-1) -1*ones(Jt(t)-1,1) temp(:,i:Jt(t)-1)];
    end 
    Mi{t}=M{ChoiceList(t)}(1:Jt(t)-1,1:Jt(t)); 
end

%Vectorized Version
P= ProbaChoice( Xs, Mi, subj , model , particle, param );
logLik = sum(log(P));

% %CellfunVersion
% proba_c = @(X) ProbaChoice( X, subj , model , particle, param );
% probas = cellfun(proba_c,Xs,'UniformOutput',0);
% for t = 1:T
%    proba_choice = probas{t};
%    logLik = logLik + log(proba_choice(ChoiceList(t)));
% end


% for t = 1:T
%    proba_choice = ProbaChoice( Xs{t}, subj , model , particle, param );
%    logLik = logLik + log(proba_choice(ChoiceList(t)));
% end

end

