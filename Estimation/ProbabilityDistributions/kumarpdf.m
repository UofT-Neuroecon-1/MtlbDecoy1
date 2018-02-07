function y = kumarpdf(x,a,b)
%KUMARPDF Kumaraswamy probability density function
%   Y = KUMARPDF(X,A,B) returns the probability density function of the
%   Kumaraswamy Distribution with shape parameters A and B.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, bounded, (0,1)
%   Restrictions:
%      A , B > 0
%
%   NOTE: The notation for the Kumaraswamy distribution used here is 
%   F(X<x) = 1 - (1-x^A)^B  
%
%   See also KUMARCDF, KUMARINV, KUMARSTAT, KUMARFIT,
%            KUMARLIKE, KUMARRND, KUMARSF, KUMARHAZ
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011


if nargin ~= 3
    error('kumarpdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('kumarpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (a>0 & a<Inf) & (b>0 & b<Inf);
okvar = (x>=0) & (x<=1);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    y(ok) = a.*b.*(x.^(a-1)).*(1-(x.^a)).^(b-1);
end


%Round off
y(y<0)=0;

end