function y=sknormpdf(x,alpha,mu,sigma)
%SKNORMPDF Skew Normal probability density function
%   Y = sknormpdf(X,ALPHA,MU,SIGMA) returns the probability density
%   function of the Skew Normal Distribution with shape parameter ALPHA,
%   location parameter MU, and scale parameter SIGMA.
%
%   SKNORMPDF(X,ALPHA) is the same as SKNORMPDF(X,ALPHA,0,1)
%   SKNORMPDF(X,ALPHA,MU) is the same as SKNORMPDF(X,ALPHA,MU,1)
%
%   Type: Continuous, unbounded
%   Restrictions:
%        SIGMA > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 5-Jun-2011


if nargin < 2
    error('sknormpdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin==2
    mu=0; sigma=1;
end

if nargin==3
    sigma=1;
end

[errorcode x alpha mu sigma] = distchck(4,x,alpha,mu,sigma);

if errorcode > 0
    error('sknormpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(alpha,'single') || isa(mu,'single') || ...
        isa(sigma,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

k=(sigma>0);
if any(k)
    temp=-alpha(k).*(x(k)-mu(k))./(sqrt(2).*sigma(k));
    term1=erfc(temp);
    term2=normpdf(x(k),mu(k),sigma(k));
    y(k)=term1.*term2;
end


end