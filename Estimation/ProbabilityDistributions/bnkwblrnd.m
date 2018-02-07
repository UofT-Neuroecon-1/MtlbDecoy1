function r = bnkwblrnd(a,b,varargin)
%BNKGBRND Random arrays from the Benktander-Weibull distribution
%   R = BNKWBLRND(A,B) returns an array of random numbers chosen from the
%   Benktander-Weibull distribution with parameters A and B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BNKWBLRND(A,B,M,N,...) or R = BNKWBLRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        A > 0
%        0 < B <= 1
%
%   Note: The Benktander-Weibull Distribution is also known as the
%   Benktander Distribution of Type II
%
%   See also BNKWBLPDF, BNKWBLCDF, BNKWBLINV, BNKWBLSTAT, 
%            BNKWBLFIT, BNKWBLLIKE
%

%   Mike Sheppard
%   Last Modified 22-Mar-2011


if nargin < 2
    error('bnkwblrnd:TooFewInputs',...
          'Requires at least two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use bnkwblinv
r = bnkwblinv(rand(varargin{:}),a,b);

end