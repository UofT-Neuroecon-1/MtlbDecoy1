function [m,v] = argusstat(chi,c)
%ARGUSSTAT Mean and variance for the ARGUS distribution
%   [M,V] = ARGUSSTAT(CHI,C) returns the mean and variance for the ARGUS
%   distribution with curvature parameter CHI, and cut-off parameter C.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default value for C is 1.
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        CHI, C > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSFIT, 
%            ARGUSLIKE, ARGUSRND, ARGUSSF, ARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 11-Dec-2011


if nargin < 1
    error('argusstat:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin == 1
    c=1;
end


try
    I = besseli(1,(chi.^2)/4); %Modified Bessel function the first kind, order 1
    m = c.*sqrt(pi/8).*chi.*exp(-(chi.^2)/4).*I./psi(chi);
    v = ((c.^2).*(1-(3./(chi.^2))+(chi.*normpdf(chi)./psi(chi))))-(m.^2);
catch err
    error('argusstat:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<c & c<Inf) & (0<chi & chi<Inf);
m(~okparam)=NaN;
v(~okparam)=NaN;

end

function y_psi=psi(x)
y_psi=normcdf(x)-x.*normpdf(x)-.5;
end

