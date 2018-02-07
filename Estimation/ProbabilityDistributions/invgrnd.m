function r = invgrnd(mu,lambda,varargin)
%INVGRND Random arrays from the Inverse Gaussian Distribution
%   R = INVGRND(MU,LAMBDA) returns an array of random numbers chosen from
%   the Inverse Gaussian Distribution with mean MU and scale parameter
%   LAMBDA.
%
%   The size of R is the common size of the parameters if all are arrays.
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = INVGRND(MU,LAMBDA,M,N,...) or R = INVGRND(MU,LAMBDA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%        MU, LAMBDA > 0
%
%   Note: The Inverse Gaussian Distribution is also known as the
%   Wald Distribution.
%
%   See also INVGPDF, INVGCDF, INVGINV, INVGSTAT, 
%            INVGFIT, INVGLIKE, INVGSF, INVGHAZ
%

%   Mike Sheppard
%   Last Modified: 12-Dec-2011

if nargin < 2
    error('invgrnd:TooFewInputs','Requires at least two input arguments.');
end

try
    %Expand size if necessary
    mu=mu+zeros(size(lambda));
    lambda=lambda+zeros(size(mu));
catch
    error('invgrnd:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

if isempty(varargin), varargin={1}; end %Scalar

sizeOut=cell2mat(varargin);

if isscalar(mu), mu = repmat(mu,sizeOut); end

c = mu.*chi2rnd(1,sizeOut);
r = (mu./(2.*lambda)) .* (2.*lambda + c - sqrt(4.*lambda.*c + c.^2));
invert = (rand(sizeOut).*(mu+r) > mu);
r(invert) = mu(invert).^2 ./ r(invert);


end