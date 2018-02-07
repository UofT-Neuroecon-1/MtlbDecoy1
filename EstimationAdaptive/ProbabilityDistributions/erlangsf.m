function s = erlangsf(x,k,lambda)
%ERLANGSF Erlang survival function
%   S = ERLANGSF(X,K,LAMBDA) returns the survival function of the 
%   Erlang Distribution with shape parameter K and rate LAMBDA
%
%   The size of S is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        K >= 1 (integer)
%        LAMBDA > 0
%
%   See also ERLANGPDF, ERLANGCDF, ERLANGINV, ERLANGSTAT, ERLANGFIT, 
%            ERLANGLIKE, ERLANGRND, ERLANGHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('erlangsf:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    s=1-erlangcdf(x,k,lambda);
catch
    error('erlangsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end