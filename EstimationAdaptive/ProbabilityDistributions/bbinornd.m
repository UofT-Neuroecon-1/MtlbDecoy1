function r = bbinornd(n,a,b,varargin)
%BBINORND Random arrays from the Beta Binomial distribution
%   R = BBINORND(N,A,B) returns an array of random numbers chosen from the
%   Beta Binomial distribution with parameters A and B, with N trials.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BBINORND(N,A,B,MM,NN,...) or R = BBINORND(N,A,B,[MM,NN,...])
%   returns an MM-by-NN-by-... array.
%
%   Type: Discrete, bounded
%   Restrictions:
%         A, B > 0
%         N >= 1       (N integer)
%
%   See also BBINOPDF, BBINOCDF, BBINOINV, BBINOSTAT, BBINOFIT, BBINOLIKE
%

%   Mike Sheppard
%   Last Modified 13-Mar-2011


if nargin < 3
    error('bbinornd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use binornd and betarnd
r = binornd(n,betarnd(a,b,varargin{:}),varargin{:});


end