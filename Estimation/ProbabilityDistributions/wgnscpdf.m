function y = wgnscpdf(x,R)
%WGNSCPDF Wigner semicircle probability density function
%   y = WGNSCPDF(x,R) returns the PDF of the Wigner semicircle distribution
%   with radius R

%   Mike Sheppard
%   Last Modified 22-Apr-2011


if nargin < 2
    error('wgnscpdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

% Return NaN for out of range parameters.
x(x<-R | x>R)=NaN;

try
    y = 2.*sqrt(1-(x./r).^2)./(pi*R);
catch
    error('wgnscpdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


end
