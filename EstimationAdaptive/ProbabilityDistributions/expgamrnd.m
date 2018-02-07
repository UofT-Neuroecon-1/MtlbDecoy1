function r = expgamrnd(k,t,mu,varargin)
%EXPGAMRND Random arrays from the Exp-Gamma Distribution.
%   R = EXPGAMRND(K,T,MU) returns an array of random numbers chosen from
%   the Exp-Gamma distribution with shape parameter K, scale paramter T,
%   and location parameter MU.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = EXPGAMRND(K,T,MU,M,N,...) or R = EXPGAMRND(K,T,MU,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, unbounded
%   Restrictions:
%        K, T > 0
%
%   Note: The Exp-Gamma Distribution is also known as the Generalized
%   Extreme Value Distribution, and also should not be confused with the
%   Log-Gamma Distribution.
%
%   See also EXPGAMPDF, EXPGAMCDF, EXPGAMINV, EXPGAMSTAT, 
%            EXPGAMFIT, EXPGAMLIKE
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 3
    error('expgamrnd:TooFewInputs',...
          'Requires at least three input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar

%Using expgaminv
r=expgaminv(rand(varargin{:}),k,t,mu);


end