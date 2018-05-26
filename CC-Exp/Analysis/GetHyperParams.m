function [ VectHyperParams ] = GetHyperParams( Particles,param, model)
%GETHYPERPARAMS Returns verctorized hyper paramaeters
VectHyperParams = cell(param.G,1);
for g = 1:param.G
    VectHyperParams{g} = struct;
    if strcmp(model,'Logit')
        VectHyperParams{g}.ha_r = nan(param.num_clust,1,param.P);
        VectHyperParams{g}.hb_r = nan(param.num_clust,1,param.P);
        VectHyperParams{g}.ha_beta = nan(param.num_clust,param.K,param.P);
        VectHyperParams{g}.hb_beta = nan(param.num_clust,param.K,param.P);
        for p = 1:param.P
            VectHyperParams{g}.ha_r(:,:,p) = Particles.particle{g,p}.ha_r;
            VectHyperParams{g}.hb_r(:,:,p) = Particles.particle{g,p}.hb_r;
            VectHyperParams{g}.ha_beta(:,:,p) = Particles.particle{g,p}.ha_beta;
            VectHyperParams{g}.hb_beta(:,:,p) = Particles.particle{g,p}.hb_beta;
        end
        
    elseif strcmp(model,'HBC-ADN')
        VectHyperParams{g}.h_gam = nan(param.size_gam,2,param.num_clust,param.P);
        VectHyperParams{g}.h_mvn_m = nan(param.size_mvn,param.num_clust,param.P);
        VectHyperParams{g}.h_mvn_s = nan(param.size_mvn,param.size_mvn,param.num_clust,param.P);
        VectHyperParams{g}.h_bet = nan(param.size_bet,2,param.num_clust,param.P);
        for p = 1:param.P
            VectHyperParams{g}.h_gam(:,:,:,p) = Particles.particle{g,p}.h_gam;
            VectHyperParams{g}.h_mvn_m(:,:,p) = Particles.particle{g,p}.h_mvn_m;
            VectHyperParams{g}.h_mvn_s(:,:,:,p) = Particles.particle{g,p}.h_mvn_s;
            VectHyperParams{g}.h_bet(:,:,:,p) = Particles.particle{g,p}.h_bet;
        end
    elseif strcmp(model,'DN3')
        VectHyperParams{g}.h_theta = nan(param.size_theta,2,param.num_clust,param.P);
        for p = 1:param.P
            VectHyperParams{g}.h_theta(:,:,:,p) = Particles.particle{g,p}.h_theta;
        end
        
    else
        error('GetHyperParams : unknown model');
    end
end

end

