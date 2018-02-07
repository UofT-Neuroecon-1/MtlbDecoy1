function y = bisapdf(x,gamma,beta,mu)
%BISAPDF Birnbaum-Saunders probability density function
%   Y=BISAPDF(X,GAMMA,BETA,MU) returns the probability density function of
%   the Birnbaum-Saunders Distribution with shape parameter GAMMA,
%   scale parameter BETA, and location parameter MU, evaluated at the
%   values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for BETA and MU are 1 and 0, respectively.
%
%   Distribution: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%        GAMMA, BETA > 0
%
%   Note: The Birnbaum-Saunders Distribution is also known as the Fatigue
%   Life Distribution. Alternative definitions have BETA defined as
%   1/LAMBDA with LAMBDA being the scale parameter. The definition used
%   here is the translation-scale notation (X-MU)/BETA
%
%   See also BISACDF, BISAINV, BISASTAT, BISAFIT, 
%            BISALIKE, BISARND, BISASF, BISAHAZ
%

%   Mike Sheppard
%   Last Modified: 17-Dec-2011


if nargin < 2
    error('bisapdf:TooFewInputs',...
        'Requires at least two input arguments.');
end
if nargin ==2
    beta=1;
    mu=0;
elseif nargin==3
    mu=0;
end

[errorcode x gamma beta mu] = distchck(4,x,gamma,beta,mu);

if errorcode > 0
    error('bisapdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
% if isa(x,'single') || isa(gamma,'single') || isa(beta,'single') || isa(mu,'single')
%     y = zeros(size(x),'single');
% else
%     y = zeros(size(x));
% end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<gamma & gamma<Inf) & (0<beta & beta<Inf) & (-Inf<mu & mu<Inf);
okvar = (mu < x) & (x < Inf);


%Translate x, but keep same variable name
x=(x-mu);

%From addbisa.m in statistics toolbox
%-----
nonpos = (x <= 0);
x(nonpos) = realmin;
z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
w = (sqrt(x./beta) + sqrt(beta./x)) ./ gamma;
ynorm = exp(-0.5 .* z.^2) ./ sqrt(2.*pi);
y = ynorm .* w ./ (2.*x);
% this would happen automatically for x==0, but generates DivideByZero warnings
y(nonpos) = 0;
%-----

y(~okparam)=NaN;
y(okparam & ~okvar)=0;

%Catch round off
y(y<0)=0;

end
