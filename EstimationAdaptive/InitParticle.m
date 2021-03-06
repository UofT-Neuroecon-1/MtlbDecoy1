function [ Particles] = InitParticle(Particles,m,J,opts)
%INITPARTICLE Returns a Theta particle drawn from prior
% Get model info

model=opts.Models{m};
K = opts.K;
N = opts.num_subj;

for g=1:opts.G
    for p=1:opts.P
%model : 'logit', 'MLBA, 'PDN' or 'PDNUnitIndep'
%K : number of attributes
        if  strcmp(model,'Logit')
            Particles.particle{g,p}.theta = nan(N,K+1);
        elseif strcmp(model,'DN1')
            Particles.particle{g,p}.theta =  nan(N,2);
            Particles.model='DN';
        elseif strcmp(model,'DN2')
            Particles.particle{g,p}.theta =  nan(N,3);
            Particles.model='DN';
        elseif strcmp(model,'PDNNew')
            Particles.particle{g,p}.theta =  nan(N,K+2);
        elseif strcmp(model,'RemiStand')
            Particles.particle{g,p}.theta =  nan(N,K+2);
        elseif strcmp(model,'HierarchicalProbit')
            %% Hyper prior
            % Prior for sigma_0 ~ Gamma(1,1), omega_0k ~ Gamma(1,1)
            Particles.particle{g,p}.hypertheta = [1,1,1,1];
            %% Particles
            % [alpha sigma Omega(1,K)]
            Particles.particle{g,p}.theta = nan(N,K+2);
        else
            error('InitParticle: unknown model used');
        end
    end
end

[Particles.model, Particles.LB, Particles.UB]=setRestrictions(model,J,opts);
Particles.r0=Particles.LB(Particles.LB==Particles.UB);
end
