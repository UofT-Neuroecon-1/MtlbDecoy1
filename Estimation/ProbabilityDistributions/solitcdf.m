function p = solitcdf(k,N)
%SOLITCDF Soliton Distribution
%   P = SOLITCDF(k,n) returns the probability density of the
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 2
    error('solitcdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

try
    %Expand size if necessary
    k=k+zeros(size(N));
    N=N+zeros(size(k));
catch
    error('solitcdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize P to NaN.
if isa(k,'single') || isa(N,'single')
    p=zeros(size(k),'single');
else
    p=zeros(size(k));
end

%Boundary case
k1=(k==1); if any(k1), p(k1)=1./N(k1); end
%Else
k2=(k>1&k<=N);
if any(k2)
    p(k2)=(1./N(k2))+(1-(1./k(k2)));
end

p(k<1)=0;
p(k>N)=1;




end