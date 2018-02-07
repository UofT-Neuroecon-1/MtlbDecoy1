function y = solitpdf(x,n)
%SOLITPDF Soliton probability density function
%   Y = SOLITPDF(x,n) returns the probability density function of the
%   Soliton Distribution with range N, at the values in x.
%
%   Type: Discrete, bounded
%   Restrictions:
%        1 <= X <= N  (both integers)
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 2
    error('solitpdf:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Expand size if necessary
    x=x+zeros(size(n));
    n=n+zeros(size(x));
catch
    error('stzpdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Return NaN for out of range parameters.
x(x<1 | x>n | x~=round(x))=NaN;
n(n~=round(n))=NaN;

% Initialize Y to zero.
if isa(x,'single') || isa(n,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

k=(x>=1 & x==round(x) & n>=1 & n==round(n) & x<=n);
if any(k)
    y(k)=1./(x(k).*(x(k)-1));
end

%Correct for x=1
x1=(x==1 & x==round(x) & n>=1 & n==round(n) & x<=n);
if any(x1)
    y(x1)=1./n(x1);
end


end