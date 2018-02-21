function r = gbetaprrnd(a,b,p,q,varargin)
%GBETAPRRND Random arrays from the Generalized Beta Prime Distribution.
%   R = GBETAPRRND(A,B,P,Q) returns an array of random numbers chosen from
%   the Generalized Beta Prime Distribution with shape parameters 
%   A, B, and P, and scale parameter Q.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = GBETAPRRND(A,B,P,Q,M,N,...) or R = GBETAPRRND(A,B,P,Q,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        A, B, P, Q > 0
%
%   Note: The Generalized Beta Prime Distribution is also known as the
%   Generalized Inverted Beta Distribution, or Generalized Beta
%   Distribution of the Second Kind
%
%   See also GBETAPRPDF, GBETAPRCDF, GBETAPRINV, GBETAPRSTAT, 
%            GBETAPRFIT, GBETAPRLIKE, GBETAPRSF, GBETAPRHAZ
%

%   Mike Sheppard
%   Last Modified 23-Mar-2011
%   Modified by Remi daviet 16-Feb-2018


if nargin < 4
    error('gbetaprrndrnd:TooFewInputs',...
          'Requires at least four input argument.'); 
end

%Transformation of Beta Prime Distribution
r = betaprrnd(a,b,varargin{:});
r=q.*(r.^(1./p));


end