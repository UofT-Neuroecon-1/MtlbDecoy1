function y = gbetaprpdf(x,a,b,p,q)
%GBETAPDF Generalized Beta Prime probability density function.
%   Y = GBETAPRPDF(X,A,B,P,Q) returns the probability density function of
%   the Generalized Beta Prime Distribution with shape parameters
%   A, B, and P, and scale parameter Q at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%      A, B, P, Q > 0
%
%   Note: The Generalized Beta Prime Distribution is also known as the
%   Generalized Inverted Beta Distribution, or Generalized Beta
%   Distribution of the Second Kind
%
%   See also GBETAPRCDF, GBETAPRINV, GBETAPRSTAT, GBETAPRFIT,
%            GBETAPRLIKE, GBETAPRRND, GBETAPRSF, GBETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin ~= 5
    error('gbetaprpdf:TooFewInputs','Requires five input arguments.');
end

[errorcode x a b p q] = distchck(5,x,a,b,p,q);

if errorcode > 0
    error('gbetaprpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single') || isa(p,'single') || isa(q,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<p & p<Inf) & (0<q & q<Inf);
okvar = (0 < x) & (x < Inf);
ok=(okparam & okvar);    
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    
    a = a(ok); b = b(ok); x = x(ok); p = p(ok); q = q(ok);
    
    % Compute logs
    xq=x./q;
    xqp=xq.^p;
    smallxqp = (xqp)<0.1;
    
    loga = (a.*p-1).*log(xq);
    logb = zeros(size(x));
    logb(smallxqp) = (-a(smallxqp)-b(smallxqp)) .* log1p(xqp(smallxqp));
    logb(~smallxqp) = (-a(~smallxqp)-b(~smallxqp)) .* log(1+xqp(~smallxqp));
    
    y(ok) = (p./q).*exp(loga+logb - betaln(a,b));
end

%Catch round off
y(y<0)=0;

end