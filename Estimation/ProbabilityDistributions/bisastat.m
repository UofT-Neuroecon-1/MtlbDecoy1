function [m,v] = bisastat(gamma,beta,mu)
%BISASTAT Mean and variance for the Birnbaum-Saunders Distribution.
%   [M,V] = BISASTAT(GAMMA,BETA,MU) returns the mean and variance of the
%   Birnbaum-Saunders Distribution with shape parameter GAMMA, scale
%   parameter BETA, and location parameter MU.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default values for BETA and MU are 1 and 0, respectively.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%        GAMMA, BETA > 0
%
%   Note: The Birnbaum-Saunders Distribution is also known as the Fatigue
%   Life Distribution. Alternative definitions have BETA defined as
%   1/LAMBDA with LAMBDA being the scale parameter. The definition used
%   here is the notation (X-MU)/BETA
%
%   See also BISAPDF, BISACDF, BISAINV, BISAFIT, BISALIKE, BISARND, 
%            BISASF, BISAHAZ
%

%   Mike Sheppard
%   Last Modified: 17-Dec-2011



if nargin < 1
    error('bisastat:TooFewInputs',...
        'Requires at least one input argument.');
end
if nargin < 2, beta = 1; end
if nargin < 3, mu = 0; end


[errorcode gamma beta mu] = distchck(3,gamma,beta,mu);

if errorcode > 0
    error('bisastat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

if isa(gamma,'single') || isa(beta,'single') || isa(mu,'single')
    m = zeros(size(gamma),'single');
else
    m = zeros(size(gamma));
end
v = m;


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<gamma & gamma<Inf) & (0<beta & beta<Inf) & (-Inf<mu & mu<Inf);
k=(okparam);

m(~k)=NaN; v(~k)=NaN;

if any(k)
    %Modified from addbisa.m in statistics toolbox
    %----
    m(k) = mu(k) + (beta(k) .* (0.5 .* gamma(k).^2 + 1));
    v(k) = (beta(k).*gamma(k)).^2 .* (1.25 .* gamma(k).^2 + 1);
    %---
end



end