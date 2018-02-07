function r = logirnd(mu,sigma,varargin)
%LOGIRND Random arrays from the Logistic Distribution
%   R = LOGIRND(MU,SIGMA) returns an array of random numbers chosen from
%   the Logistic Distribution with mean MU and scale parameter SIGMA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = LOGIRND(MU,SIGMA,M,N,...) or R = LOGIRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   See also LOGIPDF, LOGICDF, LOGIINV, LOGISTAT, LOGIFIT, LOGILIKE
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 2
    error('logirnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;


%Derived from addlogi.m in statistics toolbox
%-----
p = rand(varargin{:});
r = log(p./(1-p)).*sigma + mu;
%-----

end