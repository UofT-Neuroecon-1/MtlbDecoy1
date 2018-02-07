function parmhat = gasinfit(x)
%GASINFIT Parameter estimates for the Generalized Arcsine Distribution
%   PARMHAT = GASINFIT(X) Returns the maximum likelihood estimates of the 
%   parameters of the Generalized Arcsine Distribution given the data
%   in vector X. 
%
%   PARMHAT are the estimates for ALPHA, A, and B respectively.
%
%   Type: Continuous, bounded
%
%   See also GASINPDF, GASINCDF, GASININV, GASINSTAT, GASINLIKE, GASINRND
%

%   Mike Sheppard
%   Last Modified 3-Dec-2011


if isscalar(x)
    error('gasinfit:Vector Required');
end

% Remove missing values from the data.
x = x(~isnan(x));

aahat = min(x(:));
bbhat = max(x(:));
%normalized
x=(x(:)-aahat)./(bbhat-aahat);


% Cannot fit constant data.
if abs(aahat - bbhat) <= 2*eps(bbhat)
    error('gasinfit:BadData',...
          'Cannot fit a Generalized Arcsine Distribution if all data values are the same.');
end


%MODIFIED CODE FROM BETAFIT WITH RESTRICTION OF BETA(1-alpha,alpha)
%------------------------------
xmin = min(x); xmax = max(x);

% Initial parameter estimates.
n = length(x);
tmp1 = prod((1-x) .^ (1/n));
tmp2 = prod(x .^ (1/n));
tmp3 = (1 - tmp1 - tmp2);
ahat = 0.5*(1-tmp1) / tmp3;
pstart = ahat;

% If all values are strictly within the interval (0,1), use
% maximum likelihood with the usual continuous log-likelihood.
xl = sqrt(realmin(class(x))); % some tolerance above zero
xu = 1 - eps(class(x))/2;
if (xl <= xmin) && (xmax <= xu)
    sumlogx = sum(log(x));
    sumlog1mx = sum(log1p(-x));
    negloglike = @negloglike_cts;
    
% If some values are zero or one, maximize a mixed likelihood that
% includes discrete probabilities for those values.  Note that the asymmetry
% in xl and xu (relative to 0 and 1, respectively) means that when the data
% vector x contains exact zeros or ones, the parameter estimates for x and
% (1-x) are typically not just flipped.  But that's true even without exact
% ones and zeros, because of floating point's differing precision at 0 and 1.
else
    x0 = (x < xl);
    n0 = sum(x0);
    x1 = (x > xu);
    n1 = sum(x1);
    x2 = x(~x0 & ~x1);
    n2 = length(x2);
    sumlogx2 = sum(log(x2));
    sumlog1mx2 = sum(log1p(-x2));
    negloglike = @negloglike_mixed;
end

% Maximize the likelihood using a log transform for the parameters, to ensure
% the parameters are positive.
pstart = log(pstart);
opts = optimset('Display','off','TolX',1e-6,'TolFun',1e-6);
phat = fminsearch(negloglike,pstart,opts);
phat = exp(phat);

parmhat=[phat aahat bbhat];


    % Negative log-likelihood for data with no zeros or ones.
    function nll = negloglike_cts(p)
        p = exp(p); % remove log transform
        nll = n*betaln(1-p,p) - ((1-p)-1)*sumlogx - (p-1)*sumlog1mx;
    end

   % Negative log-likelihood for data with zeros or ones.
    function nll = negloglike_mixed(p)
        p = exp(p); % remove log transform
        
        %Contained within [0 1]
        p=max(0,min(p,1));
        
        nll = n2*betaln(1-p,p) - ((1-p)-1)*sumlogx2 - (p-1)*sumlog1mx2;
        
        % Include F(xl) = Pr(X <= xl) for data that are zeros.
        if n0 > 0
            nll = nll - n0*log(betainc(xl,1-p,p,'lower'));
        end
        
        % Include 1-F(xu) = Pr(X >= xu) for data that are ones.
        if n1 > 0
            nll = nll - n1*log(betainc(xu,1-p,p,'upper'));
        end
    end

end