function r = nakarnd(mu,omega,varargin)
%NAKARND Random arrays from the Nakagami Distribution
%   R = NAKARND(MU,OMEGA) returns an array of random numbers chosen from
%   the Nakagami Distribution with shape parameter MU and scale parameter
%   OMEGA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = NAKARND(MU,OMEGA,M,N,...) or R = NAKARND(MU,OMEGA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     MU, OMEGA > 0
%
%   See also NAKAPDF, NAKACDF, NAKAINV, NAKASTAT, NAKAFIT, NAKALIKE
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011

if nargin < 2
    error('nakarnd:TooFewInputs',...
          'Requires at two two input argument.'); 
end

if isempty(varargin), varargin={1}; end %Scalar


% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

%Modified from addnaka.m in statistics toolbox
%-----
r = sqrt(gamrnd(mu,omega./mu,varargin{:}));
%-----

end