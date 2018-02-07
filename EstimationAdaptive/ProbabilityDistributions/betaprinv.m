function x = betaprinv(p,a,b)
%BETAPRINV Inverse of the Beta Prime cumulative distribution function
%   X = BETAPRINV(P,A,B) returns the inverse cumulative distribution
%   function of the Beta Prime distribution with shape parameters A and B,
%   at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRCDF, BETAPRSTAT, BETAPRFIT,
%            BETAPRLIKE, BETAPRRND, BETAPRSF, BETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin ~= 3
    error('betaprinv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode, p, a, b] = distchck(3,p,a,b);

if errorcode > 0
    error('betaprinv:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end

%Initialize X to 0.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf);
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=0;
x(okparam & p==1)=Inf;

if any(ok)
    p=p(ok); a=a(ok); b=b(ok);
    %Use transformation of betainv
    z=betainv(p,a,b);
    x(ok)=z./(1-z);  %tranformation to Beta Prime
end

end
