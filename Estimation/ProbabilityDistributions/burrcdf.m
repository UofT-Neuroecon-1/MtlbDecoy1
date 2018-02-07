function p = burrcdf(x,c,k)
%BURRCDF Burr cumulative distribution function
%   P = BURRCDF(X,C,K) returns the cumulative distribution function of the 
%   Burr distribution with shape parameters C and K, at the values in X.  
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        C, K > 0
%
%   Note: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K   
%
%   See also BURRPDF, BURRINV, BURRSTAT, BURRFIT, 
%            BURRLIKE, BURRRND, BURRSF, BURRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011

if nargin ~= 3
    error('burrcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x c k] = distchck(3,x,c,k);

if errorcode > 0
    error('burrcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(c,'single') || isa(k,'single') 
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<k & k<Inf);
okvar = (0<x & x<Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<=0)=0;
p(okparam & x==Inf)=1;

if any(ok),
    x = x(ok); c = c(ok); k = k(ok);
    p(ok)=1-((1+(x.^c)).^(-k));
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end