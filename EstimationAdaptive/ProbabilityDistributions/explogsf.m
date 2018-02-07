function s = explogsf(x,p,b)
%EXPLOGSF Exponential-logarithmic survival function
%   S = EXPLOGSF(X,P,B) returns the survival function of the 
%   Exponential-Logarithmic Distribution with parameters P and B, 
%   at the values in X.
%
%   The size of S is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        0 < P < 1
%        B > 0
%
%   See also EXPLOGPDF, EXPLOGCDF, EXPLOGINV, EXPLOGSTAT, EXPLOGFIT, 
%            EXPLOGLIKE, EXPLOGRND, EXPLOGHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('explogsf:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    s = 1-explogcdf(x,p,b);
catch
    error('explogsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end