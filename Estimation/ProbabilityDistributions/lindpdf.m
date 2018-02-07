function y = lindpdf(x,s)
%LINDPDF Lindley probability density function.
%   Y = LINDPDF(X,S) returns the Lindley probability density
%   function with shape parameters S at the values in X.
%
%   Both X and S can be any positive real number
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 19-May-2011


if nargin < 2
   error('lindpdf:TooFewInputs','Requires two input arguments.');
end


% Return NaN for out of range parameters.
x(x<=0)=NaN; s(s<=0)=NaN;

try
    y=exp(-x.*s).*(1+x).*(s.^2)./(1+s);
catch
    error('lindpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end