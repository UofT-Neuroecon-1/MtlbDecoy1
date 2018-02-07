function s = binosf(x,n,p)
%BINOSF Binomial survival function.
%   S=BINOSF(X,N,P) returns the survival function of the binomial 
%   distribution with parameters N and P at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also BINOPDF, BINOCDF, BINOINV, BINOSTAT, 
%            BINOFIT, BINOLIKE, BINORND, BINOHAZ
%

%   Mike Sheppard
%   Last Modified 26-Dec-2011

if nargin < 3
    error('binosf:TooFewInputs',...
        'Requires three input argument.');
end


try
    s = 1 - binocdf(x,n,p);
catch
    error('binosf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end

