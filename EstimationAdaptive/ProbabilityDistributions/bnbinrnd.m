function r = bnbinrnd(n,a,b,varargin)
%BNBINRND Random arrays from the Beta Negative Binomial distribution
%   R = BNBINRND(N,A,B) returns an array of random numbers chosen from the
%   Beta Negative Binomial Distribution with parameters A and B, 
%   with N successes.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BNBINRND(N,A,B,MM,NN,...) or R = BNBINRND(N,A,B,[MM,NN,...])
%   returns an MM-by-NN-by-... array.
%
%   Type: Discrete, semi-bounded
%   Restrictions:
%         A, B > 0
%         N >= 1     (N integer)
%
%   See also BNBINPDF, BNBINCDF, BNBININV, BNBINSTAT, BNBINFIT, BNBINLIKE
%

%   Mike Sheppard
%   Last Modified 13-Mar-2011


if nargin < 3
    error('bnbinrnd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use binornd and betarnd
n(round(n)~=n | n<1)=NaN;
r = nbinrnd(n,betarnd(a,b,varargin{:}),varargin{:});


end