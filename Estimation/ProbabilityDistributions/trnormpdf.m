function y = trnormpdf(x,mu,sigma,xmin,xmax)
%TRNORMPDF Truncated Normal Probability Density Function
%   Y = TRNORMPDF(X,MU,SIGMA,XMIN,XMAX) returns the probability density of
%   the Truncated Normal Distribution with (untruncated) Normal
%   Distribution parameters of MU and SIGMA, truncated between XMIN and
%   XMAX, evaluated at the values in X.
%
%   TRNORMPDF(X,MU,SIGMA) returns the untruncated Normal Distribution
%   TRNORMPDF(X,MU,SIGMA,XMIN) returns the Left-Truncated Normal Distribution
%   For a right-truncated distribution use TRNORMPDF(X,MU,SIGMA,-Inf,XMAX)
%
%   Type: Continuous, unbounded / semi-bounded / bounded
%   Restrictions:
%        XMIN <= X <= XMAX
%        SIGMA > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 31-Jul-2011


if nargin < 3
    error('trnormpdf:TooFewInputs',...
        'Requires at least three input arguments.');
end

if nargin==3  %Untruncated
    xmin=-Inf; xmax=Inf;
end

if nargin==4  %Left truncated
    xmax=Inf;
end

[errorcode x mu sigma xmin xmax] = distchck(5,x,mu,sigma,xmin,xmax);

if errorcode > 0
    error('trnormpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(mu,'single') || isa(sigma,'single') || ...
        isa(xmin,'single') || isa(xmax,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end



k=(x>=xmin & x<=xmax & sigma>0);

if any(k)
    z1=(x(k)-mu(k))./sigma(k);
    z2=(xmax(k)-mu(k))./sigma(k);
    z3=(xmin(k)-mu(k))./sigma(k);
    %In terms of the Normal Distribution
    y(k)=(normpdf(z1)./sigma)./(normcdf(z2)-normcdf(z3));
end

%Catch out of bounds
y(x<xmin | x>xmax)=0;

end

