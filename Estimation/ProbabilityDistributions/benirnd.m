function r = benirnd(a,b,s,varargin)
%BENIRND Random arrays from the Benini distribution
%   R = BENIRND(A,B,S) returns an array of random numbers chosen from the
%   Benini distribution with shape parameters A, B and scale parameter S.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BENIRND(A,B,S,M,N,...) or R = BENIRND(A,B,S,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        A, B, S > 0
%
%   Note: The Benini Distribution is also known as log-Rayleigh Distribution
%
%   See also BENIPDF, BENICDF, BENIINV, BENISTAT, BENIFIT, BENILIKE
%

%   Mike Sheppard
%   Last Modified 13-Mar-2011


if nargin < 3
    error('benirnd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use beniinv
r = beniinv(rand(varargin{:}),a,b,s);


end
