function p = arguscdf(x,chi,c)
%ARGUSCDF ARGUS cumulative distribution function
%   P = ARGUSCDF(X,CHI,C) returns the cumulative distribution function
%   of the ARGUS distribution with curvature parameter CHI, cut-off
%   parameter C, evaluated at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSINV, ARGUSSTAT, ARGUSFIT,
%            ARGUSLIKE, ARGUSRND, ARGUSSF, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 2
    error('arguscdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    c=1;
end

[errorcode x chi c] = distchck(3,x,chi,c);

if errorcode > 0
    error('arguscdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(c,'single') || isa(chi,'single')
    p = zeros(size(x),'single');
else
    p = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
okvar=(0 < x & x < c);
ok = (okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<0)=0;
p(okparam & x>c)=1;

if any(ok)
    x=x(ok); chi=chi(ok); c=c(ok);
    rad=1-((x.^2)./(c.^2));
    p(ok)=1-(psi(chi.*(rad.^(1/2)))./psi(chi));
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end


function y_psi=psi(x)
y_psi=normcdf(x)-x.*normpdf(x)-.5;
end
