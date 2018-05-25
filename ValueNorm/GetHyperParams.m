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
        
    elseif strcmp(model,'DN')
        VectHyperParams{g}.ha_theta = nan(param.num_clust,2,param.P);
        VectHyperParams{g}.hb_theta = nan(param.num_clust,2,param.P);
        for p = 1:param.P
            VectHyperParams{g}.ha_theta(:,:,p) = Particles.particle{g,p}.ha_theta;
            VectHyperParams{g}.hb_theta(:,:,p) = Particles.particle{g,p}.hb_theta;
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

