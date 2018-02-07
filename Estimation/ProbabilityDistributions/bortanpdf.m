function y = bortanpdf(x,alpha,n)
%BORTANPDF Borel-Tanner probability density function
%   Y = BORTANPDF(X,A,N) returns the probability density function of the
%   Borel-Tanner Distribution with shape parameters ALPHA and N, at the
%   values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, semi-bounded, {N,...,Inf}
%   Restrictions:
%        1 <= N         (N integer)
%        0 < ALPHA < 1
%
%   See also BORTANCDF, BORTANINV, BORTANSTAT, BORTANFIT,
%            BORTANLIKE, BORTANRND, BORTANSF, BORTANHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011



if nargin~=3
    error('bortanpdf:TooFewInputs','Requires three input argument.');
end

[errorcode x alpha n] = distchck(3,x,alpha,n);

if errorcode > 0
    error('bortanpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(alpha,'single') || isa(n,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<alpha & alpha<1) & (n>=1 & n==round(n));
okvar = (n <= x) & (x < Inf) & (x==round(x));
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    alpha=alpha(ok); n=n(ok); x=x(ok);
    numlog=(-n+x).*log(alpha)+(-alpha.*x)+log(n)+(-1-n+x).*log(x);
    denlog=gammaln(-n+x+1);
    y(ok)=exp(numlog-denlog);
end

%Catch round off
y(y<0)=0;


end