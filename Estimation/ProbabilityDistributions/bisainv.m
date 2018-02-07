function x = bisainv(p,gamma,beta,mu)
%BISAINV Inverse of the Birnbaum-Saunders cdf
%   X = BISAINV(P,GAMMA,BETA,MU) returns the inverse cdf of the
%   Birnbaum-Saunders distribution with shape parameter GAMMA, scale 
%   parameter BETA, and location parameter MU, evaluated at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for BETA and MU are 1 and 0, respectively.
%
%   Distribution: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%        GAMMA, BETA > 0
%
%   NOTE: The Birnbaum-Saunders Distribution is also known as the Fatigue
%   Life Distribution. Alternative definitions have BETA defined as
%   1/LAMBDA with LAMBDA being the scale parameter. The definition used
%   here is the notation (X-MU)/BETA
%
%   See also BISAPDF, BISACDF, BISASTAT, BISAFIT, 
%            BISALIKE, BISARND, BISASF, BISAHAZ
%

%   Mike Sheppard
%   Last Modified: 17-Dec-2011


if nargin < 2
    error('bisainv:TooFewInputs',...
        'Requires at least two input arguments.');
end
if nargin < 3, beta = 1; end
if nargin < 4, mu = 0; end


[errorcode p gamma beta mu] = distchck(4,p,gamma,beta,mu);

if errorcode > 0
    error('bisainv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize x to zero.
if isa(p,'single') || isa(gamma,'single') || isa(beta,'single') || isa(mu,'single')
    x = zeros(size(p),'single');
else
    x = zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<gamma & gamma<Inf) & (0<beta & beta<Inf) & (-Inf<mu & mu<Inf);
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=mu(okparam & p==0);
x(okparam & p==1)=Inf;

if any(ok)
    p=p(ok); gamma=gamma(ok); beta=beta(ok); mu=mu(ok);
    %Modified from addbisa.m in statistics toolbox
    %----
    gamz = -sqrt(2).*erfcinv(2*p) .* gamma;
    x(ok) = mu + (0.25 .* beta .* (gamz + sqrt(4+gamz.^2)).^2);
end



end