function y = garguspdf(x,p,chi,c)
%GARGUSPDF Generalized ARGUS probability density function
%   Y = GARGUSPDF(X,P,CHI,C) returns the probability density function of 
%   the Generalized ARGUS distribution with power P, curvature CHI,
%   and cut-off C.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        P > -1
%        CHI, C > 0
%
%   See also GARGUSCDF, GARGUSINV, GARGUSSTAT, GARGUSFIT, 
%            GARGUSLIKE, GARGUSRND, GARGUSSF, GARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 3
    error('garguspdf:TooFewInputs',...
        'Requires at least three input arguments.');
end

if nargin == 3
    c=1;
end


[errorcode x p chi c] = distchck(4,x,p,chi,c);

if errorcode > 0
    error('garguspdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(p,'single') || isa(c,'single') || isa(chi,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
okvar=(0 < x & x < c);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;


if any(ok),
    x=x(ok); p=p(ok); c=c(ok); chi=chi(ok);
    rad=1-((x.^2)./(c.^2));
    term1=exp((-1/2).*(chi.^2).*rad);
    term2=(x./(c.^2)).*(rad.^p);
    term3=(2.^(-p)).*(chi.^(2.*(p+1)));
    term4=(gamma(1+p))-(gammainc(chi.^2/2,1+p,'upper').*gamma(1+p));
    y(ok)=(term3./term4).*term2.*term1;
end

%Catch round off
y(y<0)=0;

end

