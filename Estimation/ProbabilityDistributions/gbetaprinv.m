function x = gbetaprinv(y,a,b,p,q)
%GBETAPRINV Inverse of the Generalized Beta Prime cdf
%   X = GBETAPRINV(Y,A,B,P,Q) returns the inverse cdf of the Generalized 
%   Beta Prime distribution with shape parameters A, B, and P, and scale
%   parameter Q at the values in Y.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%        A, B, P, Q > 0
%
%   Note: The Generalized Beta Prime Distribution is also known as the
%   Generalized Inverted Beta Distribution, or Generalized Beta
%   Distribution of the Second Kind
%
%   See also GBETAPRPDF, GBETAPRCDF, GBETAPRSTAT, GBETAPRFIT, 
%            GBETAPRLIKE, GBETAPRRND, GBETAPRSF, GBETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011


if nargin < 5
    error('gbetaprinv:TooFewInputs',...
          'Requires at least five input arguments.');
end

[errorcode, y, a, b, p, q] = distchck(5,y,a,b,p,q);

if errorcode > 0
    error('gbetaprinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (0<p & p<Inf) & (0<q & q<Inf);
okvar = (0 < y) & (y < Inf);
k=(okparam & okvar);
x(~k)=NaN;
x(okparam & y==0)=0;
x(okparam & y==1)=Inf;

%Use transformation of Beta Prime
tempx=betaprinv(y(k),a(k),b(k));
x(k)=q(k).*(tempx.^(1./p(k))); %Transformation to Generalized Beta Prime

end
