function y = gbenfordpdf(x,b,n)
%GBENFORDPDF Generalized Benford probability density function
%   Y = GBENFORDPDF(X,B,N) returns the probability density function of the 
%   generalized Benford Distribution with base B, at the N'th digit, 
%   for the values of X.
%    
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   A set of numbers is said to satisfy Generalized Benford's Law if the 
%   N'th digit X={0,...,B-1} occurs with probability GBENFORDPDF(X,B,N)
%
%   Distribution: Discrete, bounded, {0,...,B-1}
%   Restrictions:
%        N >= 2         (N integer)  [Digit Position]
%        B >= 2         (B integer)  [Base]
%
%   See also GBENFORDCDF, GBENFORDINV, GBENFORDSTAT, GBENFORDFIT, 
%            GBENFORDLIKE, GBENFORDRND, GBENFORDSF, GBENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011


if nargin ~= 3
    error('gbenfordpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end


[errorcode x b n] = distchck(3,x,b,n);

if errorcode > 0
    error('gbenfordpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b)) & (n>=2 & n<Inf & n==round(n));
okvar = (0 <= x & x < b) & (x==round(x));
ok=(okparam & okvar);    
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x=x(ok); b=b(ok); n=n(ok);
    term1=gammaln((1+b.^n+x)./b);
    term2=gammaln(((b.^n)+(b.*x))./(b.^2));
    term3=gammaln((b+(b.^n)+(b.*x))./(b.^2));
    term4=gammaln((b.^n+x)./bk);
    y(ok)=((term1+term2)-(term3+term4))./log(b);
end

%Catch round off
y(y<0)=0;

end