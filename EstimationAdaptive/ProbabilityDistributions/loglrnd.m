function r = loglrnd(mu,sigma,varargin)
%LOGLRND Random arrays from the Log-Logistic Distribution
%   R = LOGLRND(MU,SIGMA) returns an array of random numbers chosen from
%   the Log-Logistic Distribution with log-mean MU and log-scale SIGMA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = LOGLRND(MU,SIGMA,M,N,...) or R = LOGLRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     SIGMA>0
%
%   See also LOGLPDF, LOGLCDF, LOGLINV, LOGLSTAT, LOGLFIT, LOGLLIKE
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 2
    error('loglrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;


%Derived from addlogi.m in statistics toolbox
%-----
p = rand(varargin{:});
r = exp(log(p./(1-p)).*sigma + mu);
%-----

end