function r = betaprrnd(a,b,varargin)
%BETAPRRND Random arrays from the Beta Prime Distribution.
%   R = BETAPRRND(A,B) returns an array of random numbers chosen from the
%   Beta Prime Distribution with shape parameters A and B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BETAPRRND(A,B,M,N,...) or R = BETAPRRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        A, B > 0
%
%   Note: The Beta Prime Distribution is also known as the
%   Inverted Beta Distribution, or Beta Distribution of the Second Kind
%
%   See also BETAPRPDF, BETAPRCDF, BETAPRINV, BETAPRSTAT, 
%            BETAPRFIT, BETAPRLIKE
%

%   Mike Sheppard
%   Last Modified 15-Mar-2011


if nargin < 2
    error('betaprrnd:TooFewInputs',...
          'Requires at least two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%X~Gamma(a,1) and Y~Gamma(b,1) then X/Y~BetaPr(a,b)
gam1=gamrnd(a,1,varargin{:});
gam2=gamrnd(b,1,varargin{:});
r = gam1./gam2;

end