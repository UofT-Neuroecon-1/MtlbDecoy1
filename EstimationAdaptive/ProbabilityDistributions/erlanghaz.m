function h = erlanghaz(x,k,lambda)
%ERLANGHAZ Erlang hazard function
%   H = ERLANGHAZ(X,K,LAMBDA) returns the hazard function of the 
%   Erlang Distribution with shape parameter K and rate LAMBDA
%
%   The size of H is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        K >= 1 (integer)
%        LAMBDA > 0
%
%   See also ERLANGPDF, ERLANGCDF, ERLANGINV, ERLANGSTAT, 
%            ERLANGFIT, ERLANGLIKE, ERLANGRND, ERLANGSF
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011

if nargin < 3
    error('erlanghaz:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    h = erlangpdf(x,k,lambda) ./ erlangsf(x,k,lambda);
catch
    error('erlanghaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end