function y = renconpdf(x,n)
%RENCONPDF Rencontres probability density function
%   Y = RENCONPDF(x,n) returns the probability density function of the
%   Recontres Distribution of having X fixed points in a uniformly
%   distributed random permutation of { 1, ..., N }.
%
%   Type: Discrete, Bounded
%   Restrictions:
%        0 <= X <= N (integer)
%        N >= 1      (integer)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011


if nargin < 2
    error('renconpdf:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Expand size if necessary
    x=x+zeros(size(n));
    n=n+zeros(size(x));
catch
    error('renconpdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(n,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN for out of range parameters.
y(x<0 | x>n | x~=round(x) | n<1 | n~=round(n))=NaN;


k=(x>=0 & x<=n & x==round(x) & n>=1 & n==round(n));
if any(k)
    xk=x(k); nk=n(k);
    %Use gammainc. The imaginary part should be zero, but may not be so due
    %to round-off. Take the real part only.
    term1=real(gammainc(-1,nk-xk+1,'upper'));
    term2=1./exp(1+gammaln(xk+1));
    y(k)=term1.*term2;
end


end