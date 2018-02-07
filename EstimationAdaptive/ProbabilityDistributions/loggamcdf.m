function p = loggamcdf(x,a,b,u)
%LOGGAMCDF Log-Gamma cumulative distribution function
%   P = LOGGAMCDF(X,A,B,U) returns the Log-Gamma cumulative distribution
%   function with shape parameters A and B and location parameter U. 
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     X>=U
%     A,B>0
%     U>=0
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 4
   error('loggamcdf:TooFewInputs','Requires four input arguments.');
end

[errorcode, x,a,b,u] = distchck(4,x,a,b,u);

if errorcode > 0
    error('loggamcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(u,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

k=(a>0 & b>0 & u>=0 & x>=u);
if any(k)
    term1=log(1+x(k)-u(k))./b(k);
    p(k)=gammainc(term1,a(k));
end

%Else not valid or zero
k1=(a>0 & b>0 & u>=0 & x<u);
p(k1)=0; p(~(k|k1))=NaN;


end