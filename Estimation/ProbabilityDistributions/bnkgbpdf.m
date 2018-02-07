function y = bnkgbpdf(x,a,b)
%BNKGBPDF Benktander-Gibrat probability density function
%   Y = BNKGBPDF(X,A,B) returns the probability density function of the
%   Benktander-Gibrat Distribution with parameters A and B at the values
%   in X.
%
%   The size of Y is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input. 
%  
%   Distribution: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%        A(A+1) >= 2B
%        A , B > 0
%
%   Note: The Benktander-Gibrat Distribution is also known as the
%   Benktander Distribution of Type I
%
%   See also BNKGBCDF, BNKGBINV, BNKGBSTAT, BNKGBFIT, 
%            BNKGBLIKE, BNKGBRND, BNKGBSF, BNKGBHAZ
%

%   Mike Sheppard
%   Last Modified 15-Dec-2011

if nargin ~= 3
    error('bnkgbpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('bnkgbpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (a.*(a+1) >= 2*b);
okvar = (1 <= x & x < Inf);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;


if any(ok),
    x = x(ok); a = a(ok); b = b(ok);
    term1=exp(-b.*(log(x)).^2);
    term2=x.^(-2-a);
    term3=-2.*b./a;
    term4=1+a+2.*b.*log(x);
    term5=1+(2.*b.*log(x)./a);
    y(ok)=term1.*term2.*(term3+term4.*term5);
end

%Catch round off
y(y<0)=0;

end