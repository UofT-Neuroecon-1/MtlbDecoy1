function p = kumarcdf(x,a,b)
%KUMARCDF Kumaraswamy cumulative distribution function
%   P = KUMARCDF(X,A,B) returns the cumulative distribution function of
%   the Kumaraswamy distribution with shape parameters A and B.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, bounded, (0,1)
%   Restrictions:
%      A , B > 0
%
%   NOTE: The notation for the Kumaraswamy distribution used here is 
%   F(X<x) = 1 - (1-x^A)^B  
%
%   See also KUMARPDF, KUMARINV, KUMARSTAT, KUMARFIT,
%            KUMARLIKE, KUMARRND, KUMARSF, KUMARHAZ
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin ~= 3
    error('kumarcdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('kumarcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (a>0 & a<Inf) & (b>0 & b<Inf);
okvar = (x>0) & (x<1);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<=0)=0;
p(okparam & x>=1)=1;


if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    p(ok) = 1 - (1-(x.^a)).^b;
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end