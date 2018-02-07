function p = betaprcdf(x,a,b)
%BETAPRCDF Beta Prime cumulative distribution function.
%   P = BETAPRCDF(X,A,B) returns the cumulative distribution function of
%   the Beta Prime Distribution with shape parameters A and B, at the
%   values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRINV, BETAPRSTAT, BETAPRFIT,
%            BETAPRLIKE, BETAPRRND, BETAPRSF, BETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 14-May-2012


if nargin~=3
    error('betaprcdf:TooFewInputs','Requires three input arguments.');
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('betaprcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    p=zeros(size(x),'single');
else
    p=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf);
okvar = (0 < x) & (x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<=0)=0;
p(okparam & x==Inf)=1;

if any(ok)
    x=x(ok); a=a(ok); b=b(ok);
    p(ok) = betainc(x./(x+1),a,b);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end