function [m,v]=szkstat(a,b)
%SZKSTAT Mean and variance for the Suzuki distribution
%   [M,V] = SZKSTAT(A,B) returns the mean and variance of the Suzuki
%   distribution with shape parameters A and B.
%
%   If X|S ~ Rayeigh(S) and S~Lognormal(A,B) then X~Suzuki(A,B)
%

%   Mike Sheppard
%   Last Modified 13-Dec-2011




if nargin < 2
    error('szkstat:TooFewInputs',...
          'Requires two input arguments.');
end


% Return NaN for out of range parameters.
b(b<0)=NaN;

try
    r=1; m1=(2.^(r/2)) .* exp(r.*a + (((r.*b).^2)/2)) .* gamma(1+(r/2));
    r=2; m2=(2.^(r/2)) .* exp(r.*a + (((r.*b).^2)/2)) .* gamma(1+(r/2));
    m=m1;
    v=m2-(m1).^2;    
catch
    error('szkstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');    
end


end