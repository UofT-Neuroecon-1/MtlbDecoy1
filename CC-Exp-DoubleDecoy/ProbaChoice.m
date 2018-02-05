function [ proba_choice ] = ProbaChoice( X , model , Theta, attrSign, param )
% X: J x K matrix of choice set
% model: the true model to use
% (returns) proba_choice : Jx1 vector of choice probabilities

J = size(X,1);
K = size(X,2);

%% function definition
%MLBA PDF(j,t|drift_rate)
% function pdf_value = PDF_LBA(d,t)
%     pdf_value = zeros(size(d,2),1);
%     cdfmin = normcdf(log(1./t),log(d)-0.5,1);
%     cdfmax = normcdf(log(2./t),log(d)-0.5,1);
%     cdfmin_psig = normcdf(log(1./t),log(d)+0.5,1);
%     cdfmax_psig = normcdf(log(2./t),log(d)+0.5,1);
%     u_LBA = cdfmax_psig - cdfmin_psig ;
%     v_LBA = cdfmax - cdfmin;
%     u_dash_LBA = (normpdf(log(1./t)-log(d)-0.5) - normpdf(log(2./t)-log(d)-0.5) ) ./t  ;
%     v_dash_LBA = (normpdf(log(1./t)-log(d)+0.5) - normpdf(log(2./t)-log(d)+0.5) ) ./t  ;
%     Z_LBA = d .* u_LBA./ v_LBA ;
%     z_LBA = d .* (u_dash_LBA .* v_LBA - u_LBA .* v_dash_LBA) ./ v_LBA.^2; 
%     F_LBA = 1 + (t.*Z_LBA - 2) .* cdfmax + (1 - t.*Z_LBA) .* cdfmin;
%     f_LBA = (Z_LBA + t .* z_LBA ) .* ( cdfmax - cdfmin ) ...
%         - (2./t.^2) .* (t.*Z_LBA-2) .* lognpdf(2./t,log(d) - 0.5,1) ...
%         - (1./t.^2) .* (1-t.*Z_LBA) .* lognpdf(1./t,log(d) - 0.5,1);
%     for j = 1:size(d,2)
%         pdf_value(j) = max([0,f_LBA(j) .* prod(1-F_LBA(1:J ~= j),'omitnan')]);
%     end
% end
%%

if  strcmp(model,'PDN')
    %True params
    alpha = Theta(1);
    sigma = Theta(2);
    Beta = (attrSign .* Theta(3:3+K-1))';
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u =  Beta .* u_x';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        norm_coefs = sum(1 ./ (sigma + u_x(j,:) + u_y),1);%./(J-1);
        v(j) = norm_coefs * unnorm_u(:,j); 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
elseif strcmp(model,'PDNNew')
    %True params
    alpha = Theta(1);
    sigma = Theta(2);
    Beta = (attrSign)';
    Omega = Theta(3:3+K-1);
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u = (attrSign .* u_x)';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        norm_coefs = sum(1 ./ (sigma + Omega .* (u_x(j,:) + u_y)),1);%./(J-1);
        v(j) = norm_coefs * unnorm_u(:,j); 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
elseif strcmp(model,'RemiStandardized')
    %True params
    alpha = Theta(1);
    sigma = Theta(2);
    Beta = (attrSign .* Theta(3:3+K-1))';
    %utility computation
    x_mean = mean(X,1);
    sd_x = std(X,[],1) * alpha + sigma;
    x_standardized = (X - x_mean) ./ sd_x;
    attr_signal = 1 ./ (1+exp(-x_standardized));
    v = attr_signal * Beta;
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
elseif strcmp(model,'RemiStandardized2')
    %True params
    alpha = Theta(1);
    sigma = Theta(2);
    Beta = (attrSign .* Theta(3:3+K-1))';
    %utility computation
    x_mean = mean(X,1);
    sd_x = std(X,[],1) + sigma;
    x_standardized = (X - x_mean) ./ sd_x;
    attr_signal = 1 ./ (1+exp(-x_standardized));
    v = attr_signal * Beta;
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
elseif strcmp(model,'RemiStandardizedL')
    %True params
    alpha = Theta(1);
    Beta = (attrSign .* Theta(2:1+K))';
    Mu0 = Theta(2+K:1+2*K);
    Sigma0 = Theta(2+2*K:1+3*K);
    %utility computation
    x_mean = alpha * mean(X,1) + (1-alpha) * Mu0;
    sd_x = alpha * std(X,[],1) + (1-alpha) * Sigma0;
    x_standardized = (X - x_mean) ./ sd_x;
    attr_signal = 1 ./ (1+exp(-x_standardized));
    v = attr_signal * Beta;
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
elseif strcmp(model,'RemiProbitStandardized')
    %True params
    sigma = Theta(1:K);
    Beta = (attrSign .*Theta(K+1:2*K))';
    alpha = Theta(end-1);
    gamma = Theta(end);
    % standardize
    x_mean = mean(X,1);
    sd_x = alpha * std(X,[],1) + gamma;
    x_standardized = (X - x_mean) ./ sd_x;
    % encode
    attr_signal = 1 ./ (1+exp(-x_standardized));
    %utility computation
    v = attr_signal * Beta;
    Cov_elements = zeros(J,J,K);
    for k=1:K
        Cov_elements(:,:,k)= attr_signal(:,k) * attr_signal(:,k)';
    end
    VarCov = sum(Cov_elements,3) + 0.1 * eye(J,J);
    %Cholesky decomp + check positive def
    [CholeskyUpper,psd] = chol(VarCov);
    while psd
        VarCov = VarCov + 0.00001 * eye(size(VarCov,1));
        [CholeskyUpper,psd] = chol(VarCov);
    end
    sim = repmat(v',[size(param.NormDraw,1) 1]) + param.NormDraw(:,1:J) * CholeskyUpper;
    proba_choice = mean(sim == max(sim,[],2));
elseif strcmp(model,'RemiProbit')
    %True params
    sigma = Theta(1:K);
    Beta = (attrSign .* Theta(K+1:2*K))';
    alpha = Theta(end);
    %utility computation
    X_norm = X.^alpha;
    u_x = sum(Beta .* X_norm')';
    Cov_elements = zeros(J,J,K);
    for k=1:K
        Cov_elements(:,:,k)= sigma(k)^2 * X_norm(:,k) * X_norm(:,k)';
    end
    VarCov = sum(Cov_elements,3) + 0.1 * eye(J,J);
    %Cholesky decomp + check positive def
    [CholeskyUpper,psd] = chol(VarCov);
    while psd
        VarCov = VarCov + 0.00001 * eye(size(VarCov,1));
        [CholeskyUpper,psd] = chol(VarCov);
    end
    sim = repmat(u_x',[size(param.NormDraw,1) 1]) + param.NormDraw(:,1:J) * CholeskyUpper;
    proba_choice = mean(sim == max(sim,[],2));
elseif strcmp(model,'RemiProbitNorm')
    %True params
    sigma = Theta(1:K);
    Beta = (attrSign .*Theta(K+1:2*K))';
    alpha = Theta(end-1);
    gamma = Theta(end);
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    norm_coefs=zeros(J,K);
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        unnorm_u =  Beta .* u_x(j,:)';
        norm_coefs(j,:) = sum(1 ./ (gamma +  u_x(j,:) + u_y),1);
        v(j) = norm_coefs(j,:) * unnorm_u; 
    end
    X_norm = u_x.*norm_coefs;
    Cov_elements = zeros(J,J,K);
    for k=1:K
        Cov_elements(:,:,k)= sigma(k)^2 * X_norm(:,k) * X_norm(:,k)';
    end
    VarCov = sum(Cov_elements,3) + 0.1 * eye(J,J);
    %Cholesky decomp + check positive def
    [CholeskyUpper,psd] = chol(VarCov);
    while psd
        VarCov = VarCov + 0.00001 * eye(size(VarCov,1));
        [CholeskyUpper,psd] = chol(VarCov);
    end
    sim = repmat(v',[size(param.NormDraw,1) 1]) + param.NormDraw(:,1:J) * CholeskyUpper;
%     sim = mvnrnd(v,VarCov,1000);
    proba_choice = mean(sim == max(sim,[],2));
elseif strcmp(model,'PDNUnitIndep')
    %True params
    alpha = Theta(1);
    sigma = Theta(2:K+1);
    Beta = Theta(K+2:end)';
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    for j=1:J %for each alternative
        u_y = u_x;
        u_y(j,:)=[]; %remove current alternative
        unnorm_u = Beta .* u_x(j,:)';
        norm_coefs = sum( 1 ./ (sigma .* (u_x(j,:) + u_y) +1) ,1);
        v(j) = norm_coefs * unnorm_u; 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
% elseif strcmp(model,'MLBA')
%     %% 
%     proba_choice = (1:3)' * 0;
%     % Debug: Theta= [1 0.8 1 1]
%     a_lba = 1;
%     x_lba = 2;
%     s_lba = 1;
%     I0 = Theta(1);
%     m = Theta(2);
%     beta = Theta(3);
%     lambda = Theta(4:5); %[positive negative]
%     vectJ = 1:J;
%     driftMean = 0*(vectJ) + I0;
%     % Corrected model
%     y_intercept = zeros(J,1);
%     x_intercept = zeros(J,1);
%     u_proj = zeros(J,2);
%     x_intercept = X(:,1) + (1/beta) *  X(:,2); %a in Trueblood's paper
%     y_intercept = X(:,1) + beta *  X(:,2); %b in Trueblood's paper
%     tan_th = X(:,2)./X(:,1);
%     u=zeros(2,2);
%     u(:,1) = y_intercept ./ (tan_th.^m + (y_intercept ./x_intercept ).^m ).^(1/m);
%     u(:,2) = x_intercept .* tan_th ./ (tan_th.^m + (y_intercept ./x_intercept ).^m ).^(1/m);
%     %Compute drift means
%     for j1=vectJ
%         vectJ2 = vectJ(vectJ~=j1);
%         for j2=vectJ2
%             sign = 1 + ([u(j1,1)-u(j2,1) ;  u(j1,2)-u(j2,2)] < 0);
%             W1 = exp( -lambda(sign(1)) * abs( u(j1,1) -u(j2,1) ) );
%             W2 = exp( -lambda(sign(2)) * abs( u(j1,2) -u(j2,2) ) );
%             driftMean(j1) = driftMean(j1) + W1 * (u(j1,1) -u(j2,1)) + W2 * ( u(j1,2) -u(j2,2) );
%         end
%     end
%     driftMean(driftMean<0)=0;
%     % uncomment for MLBA (removes C compatibility)
%     proba_choice = integral( @(t) PDF_LBA(driftMean,t) ,0.01,Inf,'ArrayValued',true);
% elseif strcmp(model,'MLBAvsPDN')
%     if Theta(1)
%         proba_choice = ProbaChoice( X , 'PDN' , Theta(2:end) );
%     else
%         proba_choice = ProbaChoice( X , 'MLBA' , Theta(2:end) );
%     end
elseif strcmp(model,'Logit')
    %True params
    alpha = Theta(1);
    Beta = (attrSign .* Theta(2:end))';
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    for j=1:J
        v(j) = sum(Beta .* u_x(j,:)'); 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
else
    proba_choice = zeros(J,1);
end
%mixture 99.9% model and 0.1% unif
proba_choice = 0.99 .* proba_choice + 0.01/J;

end

