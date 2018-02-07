function p = bnkgbcdf(x,a,b)
%BNKGBCDF Benktander-Gibrat cumulative distribution function
%   P = BNKGBCDF(X,A,B) returns the cumulative distribution function
%   of the Benktander-Gibrat distribution with parameters A and B,
%   at the values in X.
%
%   The size of P is the common sizes of the input arguments. A scalar input   
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
%   See also BNKGBPDF, BNKGBINV, BNKGBSTAT, BNKGBFIT, 
%            BNKGBLIKE, BNKGBRND, BNKGBSF, BNKGBHAZ
%  

%   Mike Sheppard
%   Last Modified 15-Dec-2011

if nargin ~= 3
    error('bnkgbcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('bnkgbcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (a.*(a+1) >= 2*b);
okvar = (1 <= x & x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<1)=0;
p(okparam & x==Inf)=1;


if any(ok)
    x = x(ok); a = a(ok); b = b(ok);
    term1=exp(-b.*(log(x)).^2);
    term2=x.^(-1-a);
    term3=1+(2.*b.*log(x)./a);
    p(ok)=1-(term1.*term2.*term3);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end