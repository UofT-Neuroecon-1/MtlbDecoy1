function s = poisssf(x,lambda)
%POISSSF Poisson survival function
%   S = POISSSF(X,LAMBDA) computes the survival function of the Poisson
%   distribution with parameter LAMBDA at the values in X.
%
%   The size of S is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also POISSPDF, POISSCDF, POISSINV, POISSSTAT, POISSFIT, POISSLIKE, 
%            POISSRND, POISSHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011
 
if nargin < 2, 
    error(message('poisssf:TooFewInputs')); 
end

try
    s = 1 - poisscdf(x,lambda);
catch
    error('poisssf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
