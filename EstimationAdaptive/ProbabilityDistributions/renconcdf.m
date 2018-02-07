function p = renconcdf(x,n)
%RENCONCDF Rencontres Cumulative Distribution Function (CDF)
%   P = RENCONCDF(x,n) returns the Recontres Cumulative Distribution Function
%   of having at most X fixed points in a uniformly distributed random
%   permutation of { 1, ..., N }.
%
%   Type: Discrete, Bounded
%   Restrictions:
%     0<=x<=n (integer)
%     n>=1    (integer)
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 2
    error('renconcdf:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Expand size if necessary
    x=x+zeros(size(n));
    n=n+zeros(size(x));
catch
    error('renconcdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(n,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

y(x<0 | x>n | x~=round(x) | n<1 | n~=round(n))=NaN;

%     0<=x<=n (integer)
%     n>=1    (integer)

k=(x>=0 & x<=n & x==round(x) & n>=1 & n==round(n));
if any(k)
    xk=x(k); nk=n(k);
    
    %use iterative sums
    
end




end