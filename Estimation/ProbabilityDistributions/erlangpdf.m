function y = erlangpdf(x,k,lambda)
%ERLANGPDF Erlang probability density function
%   Y = ERLANGPDF(X,K,LAMBDA) returns the probability density function of 
%   the Erlang Distribution with shape parameter K and rate LAMBDA, at the
%   values in X.
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        K >= 1 (integer)
%        LAMBDA > 0
%
%   See also ERLANGCDF, ERLANGINV, ERLANGSTAT, ERLANGFIT, 
%            ERLANGLIKE, ERLANGRND, ERLANGSF, ERLANGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 3
    error('erlangpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

try
    % Return NaN if any arguments are outside of their respective limits.
    k(k~=round(k))=NaN;
    
    %Let gampdf catch any other errors
    y=gampdf(x,k,1./lambda);
    
    %Catch round off
    y(y<0)=0;

catch err
        error('erlangpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end