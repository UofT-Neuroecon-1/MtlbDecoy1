function y = maxwcdf(x,a)
%MAXWCDF Maxwell–Boltzmann cumulative distribution function
%   Y = MAXWPDF(X,A) returns the Maxwell-Boltzmann cumulative distribution
%   function with scale parameter A, at the values in X.

%   Mike Sheppard
%   Last Modified 24-Mar-2011


if nargin < 2
    error('maxwcdf:TooFewInputs','Requires at least two input arguments');
end

% Return NaN for out of range parameters.
x(x<=0)=NaN; a(a<=0)=NaN;

try
    nx=x.^2./(2*a.^2);
    y=gammainc(nx,3/2);
catch
    error('maxwcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end