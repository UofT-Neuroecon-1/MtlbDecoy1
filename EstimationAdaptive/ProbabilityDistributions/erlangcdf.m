function p = erlangcdf(x,k,lambda)
%ERLANGCDF Erlang cumulative distribution function
%   P = ERLANGCDF(X,K,LAMBDA) returns the cumulative distribution function
%   of the Erlang Distribution with shape parameter K and rate LAMBDA
%
%   The size of P is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        K >= 1 (integer)
%        LAMBDA > 0
%
%   See also ERLANGPDF, ERLANGINV, ERLANGSTAT, ERLANGFIT, 
%            ERLANGLIKE, ERLANGRND, ERLANGSF, ERLANGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin < 3
    error('erlangcdf:TooFewInputs',...
          'Requires at least three input argument.'); 
end


try
    % Return NaN if any arguments are outside of their respective limits.
    k(k~=round(k))=NaN;
    
    %Let gamcdf catch any other errors
    p=gamcdf(x,k,1./lambda);
    
    %Catch round off
    p(p<0)=0; p(p>1)=1;

catch err
        error('erlangcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end