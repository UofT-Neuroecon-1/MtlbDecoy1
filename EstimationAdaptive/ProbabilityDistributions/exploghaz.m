function h = exploghaz(x,p,b)
%EXPLOGHAZ Exponential-logarithmic hazard function
%   H = EXPLOGHAZ(X,P,B) returns the hazard function of the 
%   Exponential-Logarithmic Distribution with parameters P and B, 
%   at the values in X.
%
%   The size of H is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        0 < P < 1
%        B > 0
%
%   See also EXPLOGPDF, EXPLOGCDF, EXPLOGINV, EXPLOGSTAT, 
%            EXPLOGFIT, EXPLOGLIKE, EXPLOGRND, EXPLOGSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('exploghaz:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    h = explogpdf(x,p,b) ./ explogsf(x,p,b);
catch
    error('exploghaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end