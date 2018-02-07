function [m,v] = loggamstat(a,b,u)
%LOGGAMSTAT Mean and variance for the Log-gamma distribution
%   [M,V] = LOGGAMSTAT(A,B,U) returns the mean and variance for the 
%   Log-Gamma distribution with shape parameters A and B and location
%   parameter U.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     A,B>0
%     U>=0
%
%   Mean is valid only when B<1
%   Variance is valid only when B<1/2
%
%
%   The size of the outputs is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011


if nargin < 3
   error('loggamstat:TooFewInputs','Requires three input arguments.');
end

[errorcode,a,b,u] = distchck(3,a,b,u);

if errorcode > 0
    error('loggamstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize m to zero.
if isa(a,'single') || isa(b,'single') || isa(u,'single')
   m = zeros(size(a),'single');
else
   m = zeros(size(a));
end
v=m;

k=(a>0 & b>0 & b<1);  %Mean valid for B<1
if any(k)
m(k)=-1+((1-b(k)).^(-a(k)))+u(k);
end
m(~k)=Inf;

kv=(a>0 & b>0 & b<.5); %Variance valid for B<1/2
if any(kv)
    term1=(1-2.*b(kv)).^(-a(kv));
    term2=(1-b(kv)).^(-2.*a(kv));
    v(kv)=term1-term2;
end
v(~kv)=Inf;


end