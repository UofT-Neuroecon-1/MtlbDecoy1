function [ EstimationOutput ] = EstimationAdaptiveSMC( data, opts, backup_file )
% External libraries
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath([filepath filesep 'ProbabilityDistributions'])

%% get data size
opts.num_subj = numel(data); 
M = numel(opts.Models);

%% Pre-calculate Image Matrices for Chosen Alternative
    Jm=max([data.J]);
    temp=eye(Jm-1); 
    for i=1:Jm
        Image{i}=[temp(:,1:i-1) -1*ones(Jm-1,1) temp(:,i:Jm-1)];
    end
    
    for s=1:opts.num_subj
        for t=1:length(data(s).X)
           Mi{t}=Image{data(s).y(t)}(1:data(s).J(t)-1,1:data(s).J(t)); 
        end  
        data(s).Mi=Mi;  
    end

%% Initialize particles
for m=1:M
    Particles(m).particle = cell(opts.G,opts.P);
    Particles(m).logweights = zeros(opts.G,opts.P);
    % Initialize output values
    Particles(m).log_marg_like = zeros(opts.num_subj,opts.G);
    Particles(m).log_marg_like_total = zeros(opts.G,1);
    Particles(m).group_postmeans = cell(opts.G,1);
    Particles(m).postmeans = [];
    Particles(m).group_postsd = cell(opts.G,1);
    Particles(m).model = cell(opts.G,1);
    Particles(m).LB =0;
    Particles(m).UB =0;
    Particles(m).r0 =0;
    Particles(m) = InitParticle(Particles(m),m,[data.J],opts); %Pass vector of choice set sizes, so that max can be determined.
    for g = 1:opts.G
        for p = 1:opts.P
            Particles(m).particle{g,p}.log_like_subj = zeros(opts.num_subj,1);
        end
    end
end

%% check if backup file provided and restore data
start_subj = 1;
if nargin > 2
   if exist(backup_file,'file')
       BackupData = load(backup_file);
       start_subj = BackupData.subj+1;
       Particles = BackupData.Particles;
       % check if subjects were added since backup
       if BackupData.opts.num_subj < opts.num_subj
           % Extend particles
           for m=1:M
               Particles(m).log_marg_like(opts.num_subj,1) = 0;
           end
       end
   end
end
clear BackupData
%% Check if ress threshold provided
if ~isfield(opts,'ress_threshold')
    opts.ress_threshold = 0.8;
end
%% Run SMC
for m=1:M
    for subj = 1:1%start_subj:opts.num_subj
        fprintf('Begin Subject %d\n',subj)
        num_obs = numel(data(subj).y);
        for obs = 1:num_obs
            Particles(m) = UpdateParticles( Particles(m), data, subj, obs,opts.Models{m}, opts); 
        end
        %% Save temporary file in case of crash
        save(['backup1-' opts.Tag  '.mat']);
        save(['backup2-' opts.Tag  '.mat']);
    end
end

%% Descriptive stats - Post Estimation
% save marginal likelihoods
log_marg_like = zeros(opts.G,M);
for m=1:M
    log_marg_like(:,m) = Particles(m).log_marg_like_total;
    %% Compute posterior means
    size_NK = size(Particles(m).particle{1}.theta);
    Particles(m).postmeans = nan(opts.G,size_NK(1),size_NK(2));
    for g=1:opts.G
        VectorizedTheta = nan(opts.P,size_NK(2),size_NK(1));
        for p = 1:opts.P
            VectorizedTheta(p,:,:) = Particles(m).particle{g,p}.theta';
        end
        Particles(m).group_postmeans{g} = squeeze(mean(VectorizedTheta,1))';
        Particles(m).group_postsd{g} = squeeze(std(VectorizedTheta,[],1))';
        Particles(m).postmeans(g,:,:) = Particles(m).group_postmeans{g};
    end
    Particles(m).postmeans = squeeze(mean(Particles(m).postmeans));
end

%% Save output
% create output object
EstimationOutput = struct;
% Save Data
EstimationOutput.param = opts;
EstimationOutput.Particles = Particles;
EstimationOutput.log_marg_like = log_marg_like;
%%
save(opts.savefile,'EstimationOutput','opts');

end

