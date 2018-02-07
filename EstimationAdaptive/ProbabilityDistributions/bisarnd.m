function r = bisarnd(gamma,beta,mu,varargin)
%BISARND Random arrays from the Birnbaum-Saunders Distribution
%   R = BISARND(GAMMA,BETA,MU) returns an array of random numbers chosen 
%   from the Birnbaum-Saunders Distribution with shape parameter GAMMA,
%   scale parameter BETA, and location parameter MU.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BISARND(GAMMA,BETA,MU,M,N,...) or R = BISARND(GAMMA,BETA,MU,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        GAMMA, BETA > 0
%
%   Note: The Birnbaum-Saunders Distribution is also known as the Fatigue
%   Life Distribution. Alternative definitions have BETA defined as
%   1/LAMBDA with LAMBDA being the scale parameter. The definition used
%   here is the notation (X-MU)/BETA
%
%   See also BISAPDF, BISACDF, BISAINV, BISASTAT, BISAFIT, BISALIKE
%

%   Mike Sheppard
%   Last Modified: 5-Dec-2011


if nargin < 3
    error('bisarnd:TooFewInputs','Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


% Return NaN for out of range parameters.
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;


%Modified from addbisa.m in statistics toolbox
%-----
plusminus = 2.*(rand(varargin{:})>.5) - 1; % plus or minus one, w.p. 1/2
gamz = gamma.*randn(varargin{:});
r = 0.5.*beta .* (2 + gamz.^2 + plusminus.*gamz.*sqrt(4 + gamz.^2));
%-----

%Translate r, but keep same variable name
r=(r+mu);


end
