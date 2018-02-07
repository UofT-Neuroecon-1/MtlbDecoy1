function r = laplrnd(mu,sigma,varargin)
%LAPLRND Random arrays from the Laplace Distribution
%   R = LAPLRND(MU,SIGMA) returns an array of random numbers chosen from
%   the Laplace Distribution with location MU and scale SIGMA.
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = LAPLRND(MU,SIGMA,M,N,...) or R = LAPLRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, unbounded
%   Restrictions:
%     SIGMA>0
%
%   See also LAPLPDF, LAPLCDF, LAPLINV, LAPLSTAT, LAPLFIT, LAPLLIKE
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin < 2
    error('laplrnd:TooFewInputs',...
          'Requires at least two input arguments.');
end

if isempty(varargin), varargin={1}; end %Scalar

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

%Use LAPLINV
r = laplinv(rand(varargin{:}),mu,sigma);

end