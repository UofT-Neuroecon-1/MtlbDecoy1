function y = betaprpdf(x,a,b)
%BETAPRPDF Beta Prime probability density function.
%   Y = BETAPRPDF(X,A,B) returns the probability density function of the
%   Beta Prime Distribution with shape parameters A and B, at the values
%   in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRCDF, BETAPRINV, BETAPRSTAT, BETAPRFIT,
%            BETAPRLIKE, BETAPRRND, BETAPRSF, BETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin ~= 3
    error('betaprpdf:TooFewInputs','Requires three input arguments.');
end

[errorcode, x, a, b] = distchck(3,x,a,b);

if errorcode > 0
    error('betaprpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf);
okvar = (0 < x) & (x < Inf);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;


if any(ok)
    a = a(ok); b = b(ok); x = x(ok);
    
    % Compute logs
    smallx = x<0.1;
    
    loga = (a-1).*log(x);
    logb = zeros(size(x));
    logb(smallx) = (-a(smallx)-b(smallx)) .* log1p(x(smallx));
    logb(~smallx) = (-a(~smallx)-b(~smallx)) .* log(1+x(~smallx));
    
    y(ok) = exp(loga+logb - betaln(a,b));
    
end

%Catch round off
y(y<0)=0;

end
