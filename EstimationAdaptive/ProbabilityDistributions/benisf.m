function y = benisf(x,a,b,s)
%BENISF Benini survival function
%   Y = BENISF(X,A,B,S) returns the survival function of the Benini 
%   distribution with shape parameters A and B and scale parameter S, 
%   at the values in X.
%
%   The size of Y is the common sizes of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [S,Inf)
%   Restrictions:
%        A, B, S > 0
%
%   Note: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENICDF, BENIINV, BENISTAT, BENIFIT, BENILIKE,
%            BENIRND, BENIHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 4
    error('benisf:TooFewInputs',...
          'Requires at least four input argument.'); 
end

try
    y=1-benicdf(x,a,b,s);
catch
    error('benisf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end