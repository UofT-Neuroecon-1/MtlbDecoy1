function r = dagrnd(p,a,b,varargin)
%DAGRND Random arrays from the Dagum Distribution
%   R = DAGRND(P,A,B) returns an array of random numbers chosen from the
%   Dagum Distribution with shape parameters P and A and scale parameter B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = DAGRND(P,A,B,M,N,...) or R = DAGRND(P,A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        P, A , B > 0
%
%   Note: The Dagum Distribution is also known as the Burr III
%   Distribution, and is also a special case of the Generalized 
%   Beta Prime Distribution.
%
%   See also DAGPDF, DAGCDF, DAGINV, DAGSTAT, DAGFIT, DAGLIKE
%

%   Mike Sheppard
%   Last Modified 8-May-2011

if nargin < 3
    error('dagrnd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Special case of the Generalized Beta Prime Distribution
r = gbetaprrnd(p,1,a,b,varargin{:});

end
