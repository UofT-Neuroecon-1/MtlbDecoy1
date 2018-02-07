function x = invginv(p,mu,lambda)
%INVGINV Inverse of the Inverse Gaussian cumulative distribution function
%   X = INVGINV(P,MU,LAMBDA) returns the inverse cumulative distribution 
%   function of the Inverse Gaussian distribution with mean MU and scale
%   parameter LAMBDA, evaluated at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        MU, LAMBDA > 0
%
%   Note: The Inverse Gaussian Distribution is also known as the
%   Wald Distribution.
%
%   See also INVGPDF, INVGCDF, INVGSTAT, INVGFIT, INVGLIKE, INVGRND, 
%            INVGSF, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 7-Dec-2011


if nargin < 3
    error('invginv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode p mu lambda] = distchck(3,p,mu,lambda);

if errorcode > 0
    error('invginv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(p,'single') || isa(mu,'single') || isa(lambda,'single')
    x = zeros(size(p),'single');
else
    x = zeros(size(p));
end

% Return NaN for out of range parameters.
p(p<0 | p>1)=NaN;
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;


%From addinvg.m in statistics toolbox
%-----


try
    okParams = (0 < mu) & (0 < lambda & lambda < Inf);
    k = (okParams & (0 < p & p < 1));
catch
    error(message('invginv:InputSizeMismatch'));
end
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
if ~allOK
    if isa(p,'single') || isa(mu,'single') || isa(lambda,'single')
        x = NaN(size(k),'single');
    else
        x = NaN(size(k));
    end
    x(p == 0 & okParams) = 0;
    x(p == 1 & okParams) = Inf;
    
    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(mu) > 1, mu = mu(k); end
        if numel(lambda) > 1, lambda = lambda(k); end
    else
        return;
    end
end

% Newton's Method to find a root of invgcdf(x,mu,lambda) = p
%
% Choose a starting guess for q.  Use quantiles from a lognormal
% distribution with the same mean (==1) and variance (==lambda0) as
% IG(1,lambda0).
lambda0 = lambda./mu;
sigsqLN = log(1./lambda0 + 1);
muLN = -0.5 .* sigsqLN;
q = exp(muLN - sqrt(2.*sigsqLN).*erfcinv(2*p));

% Limit the number of iterations.
maxiter = 500;
reltol = eps(class(q)).^(3/4);

F = invgcdf(q,1,lambda0);
dF = F - p;
for iter = 1:maxiter
    % Compute the Newton step, but limit its size to prevent stepping to
    % negative or infinite values.
    f = invgpdf(q,1,lambda0);
    h = dF ./ f;
    qNew = max(q/10, min(10*q, q - h)); % Avoid taking too large of a step
    
    % Break out of the iteration loop when the relative size of the last step
    % is small for all elements of q.
    done = (abs(h) <= reltol*q);
    if all(done(:))
        q = qNew;
        break
    end
    
    % Check the error, and backstep for those elements whose error did not
    % decrease.  The direction of h is always correct, because f > 0
    dFold = dF;
    for j = 1:25 % If this fails, the outer loop may too
        F = invgcdf(qNew,1,lambda0);
        dF = F - p;
        worse = (abs(dF) > abs(dFold)) & ~done;
        if ~any(worse(:))
            break
        end
        qNew(worse) = (q(worse) + qNew(worse))/2;
    end
    q = qNew;
end

badcdf = (abs(dF./F) > reltol.^(2/3));
if iter>maxiter || any(badcdf(:))   % too many iterations or cdf is too far off
    didnt = find(~done | badcdf, 1, 'first');
    if numel(mu) == 1, mubad = mu; else mubad = mu(didnt); end
    if numel(lambda) == 1, lambdabad = lambda; else lambdabad = lambda(didnt); end
    if numel(p) == 1, pbad = p; else pbad = p(didnt); end
    warning('invginv:NoConvergence',...
        'INVGINV did not converge for mu = %g, lambda = %g, p = %g.',...
        mubad,lambdabad,pbad);
end

% Add in the scale factor, and broadcast the values to the correct place if
% need be.
if allOK
    x = q .* mu;
else
    x(k) = q .* mu;
end

%-----

end



