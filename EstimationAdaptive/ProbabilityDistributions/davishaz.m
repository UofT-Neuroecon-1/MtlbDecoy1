function h = davishaz(x,b,n,u)
%DAVISCDF Davis probability density function
%   Y = DAVISCDF(X,B,N,U) returns the Davis probability density function
%   with with scale parameter B, shape parameter N, and location 
%   parameter U, at the values in X.
%
%   Restrictions:
%   B: any positive real number
%   N: great than 1
%   U: any non-negative real number
%   Valid for X>U

%   The size of Y is the common size of the input arguments. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%

%   Mike Sheppard
%   Last Modified 9-May-2011

if nargin < 4
    error('davishaz:TooFewInputs',...
          'Requires at least four input argument.'); 
end

try
    yt = davispdf(x,b,n,u);
    st = davissf(x,b,n,u);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch
    error('davishaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end