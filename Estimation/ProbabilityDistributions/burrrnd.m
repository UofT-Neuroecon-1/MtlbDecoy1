function r = burrrnd(c,k,varargin)
%BURRRND Random arrays from the Burr distribution
%   R = BURRRND(C,K) returns an array of random numbers chosen from the
%   Burr distribution with shape parameters C and K.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BURRRND(C,K,M,N,...) or R = BURRRND(C,K,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        C, K > 0
%
%   Note: The Burr Distribution is also known as the Burr Type XII
%   Distribution, and is a special case of the Sing-Maddala Distribution.
%   The notation for the distribution used here is F(X<x) = 1 - (1+x^C)^-K
%
%   See also BURRPDF, BURRCDF, BURRINV, BURRSTAT, BURRFIT, BURRLIKE
%

%   Mike Sheppard
%   Last Modified 25-Mar-2011


if nargin < 2
    error('burrrnd:TooFewInputs',...
          'Requires at least two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use burrinv
r = burrinv(rand(varargin{:}),c,k);


end