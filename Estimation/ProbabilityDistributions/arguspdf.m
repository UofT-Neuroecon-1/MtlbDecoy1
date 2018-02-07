function y = arguspdf(x,chi,c)
%ARGUSPDF ARGUS probability density function
%   Y = ARGUSPDF(X,CHI,C) returns the probability density function of the
%   ARGUS Distribution with curvature CHI, cut-off parameter C, evaluated
%   at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSCDF, ARGUSINV, ARGUSSTAT, ARGUSFIT, 
%            ARGUSLIKE, ARGUSRND, ARGUSSF, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 2
    error('arguspdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

if nargin == 2
    c=1;
end


[errorcode x chi c] = distchck(3,x,chi,c);

if errorcode > 0
    error('arguspdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(c,'single') || isa(chi,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
okvar=(0 < x) & (x < c);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x=x(ok); c=c(ok); chi=chi(ok);
    rad=1-((x.^2)./(c.^2));
    y(ok)=(chi.^3./(sqrt(2*pi)*psi(chi))).*(x./(c.^2)).*...
        (rad.^(1/2)).*exp((-1/2).*(chi.^2).*rad);
end

%Catch round off
y(y<0)=0;


end



function y_psi=psi(x)
y_psi=normcdf(x)-x.*normpdf(x)-.5;
end
