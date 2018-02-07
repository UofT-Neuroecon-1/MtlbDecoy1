function r = tlsrnd(mu,sigma,nu,varargin)
%TLSRND Random arrays from the T Location-Scale Distribution
%   R = TLSRND(MU,SIGMA,NU) returns an array of random numbers chosen 
%   from the T Location-Scale Distribution with location parameter MU, 
%   scale parameter SIGMA, and NU degrees of freedom.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = TLSRND(MU,SIGMA,NU,M,N,...) or R = TLSRND(MU,SIGMA,NU,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA, NU > 0
%
%   See also TLSPDF, TLSCDF, TLSINV, TLSSTAT, TLSFIT, TLSLIKE
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 3
    error('tlsrnd:TooFewInputs','Requires at least three input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar


% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;


%From addtls.m in statistics toolbox
%-----
r = mu + sigma.*trnd(nu,varargin{:});
%-----


end