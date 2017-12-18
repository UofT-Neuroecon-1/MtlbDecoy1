function [ theta ] = InitParticle( model , K )
%INITPARTICLE Returns a Theta particle drawn from prior
%model : 'logit', 'MLBA, 'PDN' or 'PDNUnitIndep'
%K : number of attributes
    if  strcmp(model,'Logit')
        theta = [betarnd(3,1) gamrnd(2,2,1,K)];
    elseif strcmp(model,'PDN')
        theta = [betarnd(4,1) gamrnd(1,0.5,1,1) gamrnd(2,2,1,K)];
    elseif strcmp(model,'PDNUnitIndep')
        theta = [betarnd(3,1) gamrnd(1,0.5,1,K) gamrnd(2,2,1,K)];
    elseif strcmp(model,'RemiProbit')
        theta = [gamrnd(1,0.5,1,K) gamrnd(2,2,1,K)];
    elseif strcmp(model,'RemiProbitNorm')
        theta = [gamrnd(1,0.5,1,K) gamrnd(2,2,1,K) gamrnd(1,1) gamrnd(1,0.5,1,1)];
    elseif strcmp(model,'MLBA')
        theta = [gamrnd(1,1,1,1) betarnd(3,1)  gamrnd(1,1,1,3)];
    elseif strcmp(model,'MLBAvsPDN')
        model = binornd(1,0.5);
        if model
            theta = [1 InitParticle( 'PDN' , K) ];
        else
            theta = [0 InitParticle( 'MLBA' , K)];
        end
    else
        theta = nan;
    end

end

