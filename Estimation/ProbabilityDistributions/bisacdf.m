function p = bisacdf(x,gamma,beta,mu)
%BISACDF Birnbaum-Saunders cumulative distribution function
%   P=BISACDF(X,GAMMA,BETA,MU) returns the cumulative distribution function
%   of the Birnbaum-Saunders Distribution with shape parameter GAMMA,
%   scale parameter BETA, and location parameter MU, evaluated at the
%   values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
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
%   here is the notation (X-MU)/BETA
%
%   See also BISAPDF, BISAINV, BISASTAT, BISAFIT, 
%            BISALIKE, BISARND, BISASF, BISAHAZ
%

%   Mike Sheppard
%   Last Modified: 17-Dec-2011


if nargin < 2
    error('bisacdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end
if nargin==2
    beta=1; mu=0;
elseif nargin==3
    mu=0;
end


[errorcode x gamma beta mu] = distchck(4,x,gamma,beta,mu);

if errorcode > 0
    error('bisacdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(gamma,'single') || isa(beta,'single') || isa(mu,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<gamma & gamma<Inf) & (0<beta & beta<Inf) & (-Inf<mu & mu<Inf);
%okvar = (mu < x) & (x < Inf);

%Translate x, but keep same variable name
x=(x-mu);

%From addbisa.m in statistics toolbox
%-----
nonpos = (x <= 0);
x(nonpos) = realmin;
z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
p = 0.5 * erfc(-z ./ sqrt(2));
% this would happen automatically for x==0, but generates DivideByZero warnings
p(nonpos) = 0;
%-----

%Remember, x is now redefined as translated by mu
p(~okparam)=NaN;
p(okparam & (x+mu)<=mu)=0;
p(okparam & (x+mu)==Inf)=1;

%Catch round off
p(p<0)=0; p(p>1)=1;

end
