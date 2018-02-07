function y = zetapdf(x,p)
%ZETAPDF Zeta probability density function
%   Y = ZETAPDF returns the Zipf density function with
%   with parameter P
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 25-May-2011


if nargin < 2
    error('zetapdf:TooFewInputs',...
          'Requires at least two input arguments.'); 
end


% Return NaN for out of range parameters.
p(p<=0)=NaN; x(x<1 | x~=round(x))=NaN;

try
   y=(x.^(-1-p))./zeta(1+p);
catch
    error('zetapdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');    
end


end