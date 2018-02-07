function x = erlanginv(p,k,lambda)
%ERLANGINV Inverse of the Erlang cumulative distribution function
%   X = ERLANGINV(P,K,LAMBDA) returns the inverse cumulative distribution
%   function of the Erlang Distribution with shape parameter K and
%   rate LAMBDA
%
%   The size of X is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%      K >= 1 (integer)
%      LAMBDA > 0
%
%   See also ERLANGPDF, ERLANGCDF, ERLANGSTAT, ERLANGFIT, 
%            ERLANGLIKE, ERLANGRND, ERLANGSF, ERLANGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 3
    error('erlanginv:TooFewInputs',...
        'Requires three input arguments.');
end


try
    % Return NaN if any arguments are outside of their respective limits.
    k(k~=round(k))=NaN;
    
    %Let gaminv catch any other errors
    x=gaminv(p,k,1./lambda);
    
catch err
    error('erlanginv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


end