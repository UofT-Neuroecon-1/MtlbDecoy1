function y=pascalpdf(x,n,p)
%PASCALPDF Pascal's Distribution probability density function
%   Y = PASCALPDF(X,N,P) returns the probability density function of the
%   Pascal Distribution of the number of trials with success probability
%   P before N successes occur, at the values of X.
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      x >= n
%      n >= 0 (integer)
%      0 <= p < =1
%
%   The size of Y is the size of the input variable X.
%

%   Mike Sheppard
%   Last Modified 18-Jun-2011


if nargin < 3
    error('pascalpdf:TooFewInputs',...
          'Requires three input argument.');
end


[errorcode, x, n,p] = distchck(3,x,n,p);

if errorcode > 0
    error('pascalpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN for out of range parameters.
x(x<n | x~=round(x))=NaN; n(x<n | n<0 | n~=round(n))=NaN; p(p<0 | p>1)=NaN;

%PASCALPDF(X,N,P)=NBINPDF(X-N,N,P)
%Let nbinpdf handle all error checking
y=nbinpdf(x-n,n,p);

%Round off
y(y<0)=0;

end