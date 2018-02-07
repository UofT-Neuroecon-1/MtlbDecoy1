function x = kumarinv(p,a,b)
%KUMARINV Inverse of the Kumaraswamy cumulative distribution function
%   X = KUMARINV(P,A,B) returns the inverse of the cumulative distribution
%   function of the Kumaraswamy distribution with shape parameters A and B.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Type: Continuous, bounded
%   Restrictions:
%      A , B > 0
%
%   Note: The notation for the Kumaraswamy distribution used here is 
%   F(X<x) = 1 - (1-x^A)^B  
%
%   See also KUMARPDF, KUMARCDF, KUMARSTAT, KUMARFIT, KUMARLIKE, KUMARRND
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin < 3
    error('kumarinv:TooFewInputs',...
          'Requires two input arguments.');
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('kumarinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k=(p>=0 & p<=1 & a>0 & b>0);

if any(k)
    x(k)=(1-((1-p(k)).^(1./b(k)))).^(1./a(k));
end

%Edge cases and round off
x(~k)=NaN;

end