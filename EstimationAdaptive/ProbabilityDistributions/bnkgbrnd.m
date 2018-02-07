function r = bnkgbrnd(a,b,varargin)
%BNKGBRND Random arrays from the Benktander-Gibrat distribution
%   R = BNKGBRND(A,B) returns an array of random numbers chosen from the
%   Benktander-Gibrat distribution with parameters A and B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BNKGBRND(A,B,M,N,...) or R = BNKGBRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        A(A+1) >= 2B
%        A , B > 0
%
%   Note: The Benktander-Gibrat Distribution is also known as the
%   Benktander Distribution of Type I.
%
%   See also BNKGBPDF, BNKGBCDF, BNKGBINV, BNKGBSTAT, BNKGBFIT, BNKGBLIKE
%

%   Mike Sheppard
%   Last Modified: 22-Mar-2011


if nargin < 2
    error('bnkgbrnd:TooFewInputs',...
          'Requires at least two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Use bnkgbinv
r = bnkgbinv(rand(varargin{:}),a,b);


end
