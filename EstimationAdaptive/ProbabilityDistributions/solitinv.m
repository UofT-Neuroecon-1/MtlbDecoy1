function k = solitinv(p,N)
%SOLITINV Soliton Distribution
%   K = SOLITINV(p,n) returns the probability density of the
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 2
    error('solitcdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

try
    %Expand size if necessary
    p=p+zeros(size(N));
    N=N+zeros(size(p));
catch
    error('solitcdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize P to NaN.
if isa(p,'single') || isa(N,'single')
    k=zeros(size(p),'single');
else
    k=zeros(size(p));
end

term=N./(1+N.*(1-p));

%round off those that need rounding
%Usually round up, unless within eps of answer
k1=abs(term-floor(term));
k2=(k1<10*eps);
if any(k2)
    term(k2)=floor(term(k2));
end
k=ceil(term);

end