function y = wgnsccdf(x,R)
%WGNSCCDF Wigner semicircle cumulative distribution function
%   y = WGNSCCDF(x,R) returns the CDF of the Wigner semicircle distribution
%   with radius R

%   Mike Sheppard
%   Last Modified 13-Dec-2011


if nargin < 2
    error('wgnsccdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

% Return NaN for out of range parameters.
x(x<-R | x>R)=NaN;

try
    y = (1/2) + (x.*sqrt(1-(x./R).^2)./(pi*R)) + (asin(x./R)./pi);
catch
    error('wgnsccdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


end
