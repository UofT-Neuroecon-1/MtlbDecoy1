function [m,v] = erlangstat(k,a)
%ERLANGSTAT Mean and variance for the Erlang distribution
%   [M,V] = ERLANGSTAT(K,A) returns the mean and variance of the Erlang
%   distribution with shape parameter K and rate A.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.   
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        K >= 1 (integer)
%        LAMBDA > 0
%
%   See also ERLANGPDF, ERLANGCDF, ERLANGINV, ERLANGFIT, 
%            ERLANGLIKE, ERLANGRND, ERLANGSF, ERLANGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin < 2
    error('erlangstat:TooFewInputs',...
          'Requires at least two input argument.');
end

% Return NaN for out of range parameters.
k(k~=round(k)) = NaN;

try
    m = k ./ a;
    v = k ./ (a.^2);
catch
    error('erlangstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


end
