function y = uquadpdf(x,a,b)
%UQUADPDF U-quadratic probability density function
%   Y = UQUADPDF(X,A,B) returns the probability density function of the
%   U-quadratic Distribution with lower limit A and upper limit B, 
%   evaluated at the values in X
%
%   UQUADPDF(X) uses the default domain of a=0, b=1
%
%   Type: Continuous, bounded
%   Restrictions:
%        A <= X <= B
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if (nargin ~= 1) && (nargin~=3)
    error('uquadpdf:TooFewInputs',...
          'Requires either one or three input arguments.'); 
end

if nargin==1
    a=0; b=1;
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('uquadpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Return NaN for out of range parameters.
a(a>b) = NaN;
b(a>b) = NaN;

try
   alpha=12./((b-a).^3); %vertical scale
   beta=(b+a)./2; %gravitational balance center, offset
   y=alpha.*((x-beta).^2);
catch
    error('uquadpdf:InputSizeMismatch');
end

%Catch out of bounds
y(x<a)=0;
y(x>b)=0;

end