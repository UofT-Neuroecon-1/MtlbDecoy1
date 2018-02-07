function y = sgompdf(x,b,n)
%SGOMPDF Shifted Gompertz probability density function
%   Y = SGOMPDF(x,b,n) returns the probability density function of the
%   Shifted Gompertz Distribution with scale parameter B,
%   shape parameter N, at the values in X.
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%        X, B, N > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 10-Jun-2011


if nargin < 3
    error('sgompdf:TooFewInputs',...
        'Requires three input arguments.');
end


[errorcode,x,b,n] = distchck(3,x,b,n);

if errorcode > 0
    error('sgompdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN for out of range parameters.
y(x<0 | b<0 | n<0)=NaN;

%Out of bounds
b(b<0)=NaN; n(n<0)=NaN;

k=(x>0 & b>0 & n>0);
if any(k)
    term1=exp(-b(k).*x(k));
    term2=exp(-n(k).*term1);
    y(k)=b(k).*term1.*term2.*(1+n(k).*(1-term1));
end

%round off
y(y<0)=0;


end