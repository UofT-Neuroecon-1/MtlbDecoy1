function x = hmtinv(p,n)
%HMTINV Inverse of the Heads-Minus-Tails cumulative distribution function
%   X = HMTINV(P,N) returns the inverse cumulative distribution function
%   of the Heads-Minus-Tails Distribution, given a fair coin is tossed 2N
%   number of times, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Since the Heads-Minus-Tails distribution is discrete, HMTINV
%   returns the least integer X such that the Heads-Minus-Tails cdf
%   evaluated at X, equals or exceeds P.
%
%   Type: Discrete, bounded, {0,...,N}
%   Restrictions:
%        N >= 1    (integer)
%
%   See also HMTPDF, HMTCDF, HMTSTAT, HMTFIT, HMTLIKE, HMTRND, 
%            HMTSF, HMTHAZ
%

%   Mike Sheppard
%   Last Modified 24-Dec-2011




if nargin < 2
    error('hmtinv:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Expand size if necessary
    p=p+zeros(size(n));
    n=n+zeros(size(p));
catch
    error('hmtinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN for out of range parameters.
n(n~=round(n) | n<1)=NaN;  p(p<0 | p>1)=NaN;

k=1:numel(p);
k1=find(p<0 | p>1);
k2=find(p==1);

%Initialize X to zero.
if isa(p,'single') || isa(n,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end



cumdist=x;
x(k2)=n(k2);  %Equals to N
x(k1)=NaN;
k([k1(:);k2(:)])=[];
if isempty(k), return; end

cumdist(k)=hmtpdf(0,n(k));

count=0;

k=k(cumdist(k)<p(k));
while ~isempty(k)
    x(k)=x(k)+1;
    count=count+1;
    cumdist(k)=cumdist(k)+hmtpdf(count,n(k));
    k=k(cumdist(k)<p(k));
end


end