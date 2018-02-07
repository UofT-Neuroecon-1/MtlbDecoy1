function r = ricernd(s,sigma,varargin)
%RICERND Random arrays from the Rician Distribution
%   R = RICERND(S,SIGMA) returns an array of random numbers chosen 
%   from the Rician Distribution with noncentrality parameter S and scale
%   parameter SIGMA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = RICERND(S,SIGMA,M,N,...) or R = RICERND(S,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%     S >= 0
%     SIGMA > 0
%
%   See also RICEPDF, RICECDF, RICEINV, RICESTAT, RICEFIT, RICELIKE
%

%   Mike Sheppard
%   Last Modified 14-Dec-2011


if nargin < 2
    error('ricernd:TooFewInputs','Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


% Return NaN for out of range parameters.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;


%Modified from addrice.m in statistics toolbox
%-----
r = sigma .* sqrt(ncx2rnd(2, (s./sigma).^2, varargin{:}));
%-----


end