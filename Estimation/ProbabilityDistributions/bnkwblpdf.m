function y = bnkwblpdf(x,a,b)
%BNKWBLPDF Benktander-Weibull probability density function
%   Y = BNKWBLPDF(X,A,B) returns the probability density function of the 
%   Benktander-Weibull Distribution with parameters A and B at the 
%   values in X.
%
%   The size of Y is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input. 
%
%   Distribution: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A > 0
%        0 < B <= 1
%
%   Note: The Benktander-Weibull Distribution is also known as the
%   Benktander Distribution of Type II
%
%   See also BNKWBLCDF, BNKWBLINV, BNKWBLSTAT, BNKWBLFIT, 
%            BNKWBLLIKE, BNKWBLRND, BNKWBLSF, BNKWBLHAZ
%

%   Mike Sheppard
%   Last Modified: 13-May-2011

if nargin ~= 3
    error('bnkwblpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('bnkwblpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<=1);
okvar = (1 <= x & x < Inf);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x = x(ok); a = a(ok); b = b(ok);
    term1=exp(a.*(1-(x.^b))./b);
    term2=x.^(-2+b);
    term3=(1-b+a.*(x.^b));
    y(ok)=term1.*term2.*term3;
end

%Catch round off
y(y<0)=0;

end