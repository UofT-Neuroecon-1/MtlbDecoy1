function [ logLik ] = LogLikelihood( Xs, ChoiceList, subj , model , particle, opts, algoVars )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here

T = numel(Xs);
Jt=cellfun(@length,Xs);

dataobs.X = Xs;
dataobs.Mi = {};
for t=1:T
    [dataobs.J(t), dataobs.K(t)] = size(dataobs.X{t});
    temp=eye(Jt(t)-1); 
    for i=1:Jt(t)
        M{i}=[temp(:,1:i-1) -1*ones(Jt(t)-1,1) temp(:,i:Jt(t)-1)];
    end 
    dataobs.Mi{t}=M{ChoiceList(t)}(1:Jt(t)-1,1:Jt(t)); 
end


par(algoVars.LB==algoVars.UB)=algoVars.r0;%Set the restricted variables.
if any(algoVars.LB~=algoVars.UB)
    par(algoVars.LB~=algoVars.UB)=particle.theta;%Set the unrestricted variables
end
expPart.theta = par;
Pi=ProbaChoice(dataobs, expPart, algoVars.model, opts );
logLik = sum(log(Pi));

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

