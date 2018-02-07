function [m,v]=uquadstat(a,b)
%UQUADSTAT Mean and variance of the U-quadratic distribution
%   [M,V]=UQUADSTAT(A,B) returns the mean and variance of the U-quadratic
%   distribution with lower limit A and upper limit B.
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 2
    error('uquadstat:TooFewInputs',...
          'Requires two input arguments.'); 
end


try
    %Expand size if necessary
    a=a+zeros(size(b));
    b=b+zeros(size(a));
catch
    error('uquadstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Return NaN for out of range parameters.
a(a>b) = NaN; b(a>b) = NaN;

m=(a+b)/2;
v=(3/20).*((b-a).^2);


end