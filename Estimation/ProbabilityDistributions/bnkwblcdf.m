function p = bnkwblcdf(x,a,b)
%BNKWBLCDF Benktander-Weibull cumulative distribution function
%   P = BNKWBLCDF(X,A,B) returns the cumulative distribution function of 
%   the Benktander-Weibull distribution with parameters A and B at the
%   values in X.
%
%   The size of P is the common sizes of the input arguments. A scalar input   
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
%   See also BNKWBLPDF, BNKWBLINV, BNKWBLSTAT, BNKWBLFIT, 
%            BNKWBLLIKE, BNKWBLRND, BNKWBLSF, BNKWBLHAZ
%

%   Mike Sheppard
%   Last Modified: 16-Dec-2011

if nargin ~= 3
    error('bnkwblcdf:TooFewInputs',...
          'Requires at least three input argument.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('bnkwblcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<=1);
okvar = (1 <= x & x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<1)=0;
p(okparam & x==Inf)=1;


if any(ok)
    x = x(ok); a = a(ok); b = b(ok);
    p(ok) = 1 - ((exp(a.*(1-(x.^b))./b)).*(x.^(-1+b)));
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end