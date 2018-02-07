function r = invgamrnd(a,b,varargin)
%INVGAMRND Random arrays from the Inverse-Gamma Distribution
%   R = INVGAMRND(A,B) returns an array of random numbers chosen from the
%   Inverse-Gamma Distribution with shape parameter A and scale
%   parameter B.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = INVGAMRND(A,B,M,N,...) or R = INVGAMRND(A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded, (0,Inf)
%   Restrictions:
%     A, B > 0
%
%   See also INVGAMPDF, INVGAMCDF, INVGAMINV, INVGAMSTAT, 
%            INVGAMFIT, INVGAMLIKE, INVGAMSF, INVGAMHAZ
%

%   Mike Sheppard
%   Last Modified 25-Jun-2011



if nargin < 2
    error('invgamrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


%Use invgaminv
r = invgaminv(rand(varargin{:}),a,b);

end