function p = uquadcdf(x,a,b)
%UQUADCDF U-quadratic cumulative distribution function
%   P = UQUADCDF(X,A,B) returns the U-quadratic cumulative distribution function
%   with lower limit A and upper limit B, evaluated at the values in X
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('uquadcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('uquadcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Return NaN for out of range parameters.
a(a>b) = NaN;
b(a>b) = NaN;

try
   alpha=12./((b-a).^3); %vertical scale
   beta=(b+a)./2; %gravitational balance center, offset
   p=(alpha./3).*(((x-beta).^3)+((beta-a).^3));
catch
    error('uquadcdf:InputSizeMismatch');
end



end