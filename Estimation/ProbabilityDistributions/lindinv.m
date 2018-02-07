function x = lindinv(p,s)
%LINDINV Inverse of the Lindley cumulative distribution function (cdf)
%   X = LINDINV(P,S) returns the inverse of the Lindley distribution with
%   shape parameter S at the values in P.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   USES LAMBERTW FUNCTION, NEWTON ELSE
%
%   Mike Sheppard
%   Last Modified: 20-May-2011


if nargin < 2
    error('lindinv:TooFewInputs',...
        'Requires at least two input arguments.');
end


% Return NaN for out of range parameters.
p(p<0 | p>1)=NaN; s(s<0)=NaN;

try
    x = (-1-s-lambertw(-1,exp(-1-s).*(-1+p).*(1+s)))./(s);
catch
    error('lindinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end



end