function p = sgomcdf(x,b,n)
%SGOMCDF Shifted Gompertz cumulative distribution function
%   Y = SGOMPDF(x,b,n) returns the shifted gompertz probability density
%   function with scale parameter B and shape parameter N at the values in
%   X.
%
%   X,B,N are all real numbers great than 0.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 10-Jun-2011


if nargin < 3
    error('sgomcdf:TooFewInputs',...
          'Requires three input arguments.');
end


[errorcode,x,b,n] = distchck(3,x,b,n);

if errorcode > 0
    error('sgomcdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

%Out of bounds
b(b<0)=NaN; n(n<0)=NaN;

term1=exp(-b.*x);
term2=exp(-n.*term1);
p=(1-term1).*term2;


end