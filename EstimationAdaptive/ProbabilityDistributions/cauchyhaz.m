function h = cauchyhaz(x,a,b)
%CAUCHYHAZ Cauchy hazard function
%   H = CAUCHYHAZ(X,A,B) returns the hazard function of the 
%   Cauchy Distribution with location parameter A and scale
%   parameter B, at the values in X. 
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, unbounded, (-Inf,Inf)
%   Restrictions:
%        B > 0
%
%   Note: The Cauchy Distribution is also known as the Lorentz Distribution
%
%   See also CAUCHYPDF, CAUCHYCDF, CAUCHYINV, CAUCHYSTAT, 
%            CAUCHYFIT, CAUCHYLIKE, CAUCHYRND, CAUCHYSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin<1
    error('cauchyhaz:TooFewInputs','Input argument X is undefined.');
end
if nargin<2, a=0; end
if nargin<3, b=1; end


try
    h = cauchypdf(x,a,b) ./ cauchysf(x,a,b);
catch
    error('cauchyhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end