function p = gbenfordcdf(x,b,n)
%GBENFORDCDF Generalized Benford cumulative distribution function
%   P = GBENFORDCDF(X,B,N) returns the cumulative distribution function 
%   of the Generalized Benford distribution with parameters of base B, 
%   at the N'th digit, at the values of X.
%    
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   A set of numbers is said to satisfy Generalized Benford's Law if the 
%   N'th (N>=2) digit X={0,...,B-1} occurs with probability 
%   GBENFORDPDF(X,B,N)
%
%   Distribution: Discrete, bounded, {0,...,B-1}
%   Restrictions:
%        N >= 2         (N integer)  [Digit Position]
%        B >= 2         (B integer)  [Base]
%
%   See also GBENFORDPDF, GBENFORDINV, GBENFORDSTAT, GBENFORDFIT, 
%            GBENFORDLIKE, GBENFORDRND, GBENFORDSF, GBENFORDHAZ
%

%   Mike Sheppard
%   Last Modified: 14-May-2012


if nargin ~= 3
    error('gbenfordcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end


[errorcode x b n] = distchck(3,x,b,n);

if errorcode > 0
    error('gbenfordcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b)) & (n>=2 & n<Inf & n==round(n));
x=floor(x); %round down for CDF for discrete
okvar = (0 <= x & x < b);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<0)=0;
p(okparam & x>=b)=1;


scalarbn=isscalar(b);
if any(ok)
    x=x(ok); b=b(ok); n=n(ok);
    val=max(x(:));
    I=(0:val)';
    if scalarbn
        tmp=cumsum(gbenfordpdf(I,b(1),n(1)));
        p(ok)=tmp(x+1);
    else
        n1=sum(ok); n2=numel(I);
        compare=repmat(I,1,n1);
        index=repmat(x,n2,1);
        bbig=repmat(b,n2,1);
        nbig=repmat(n,n2,1);
        y0=gbenfordpdf(compare,bbig,nbig);
        y0(compare>index)=0;
        p(ok)=sum(y0,1);
    end
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end