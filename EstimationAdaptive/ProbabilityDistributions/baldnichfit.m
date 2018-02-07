function phat = baldnichfit(x)
%BALDNICHFIT Parameter estimates for Balding-Nichols distribution
%   PHAT = BALDNICHFIT(X) Returns the maximum likelihood estimates of the
%   parameters of the Balding-Nichols distribution given the data in the
%   vector X.
%
%   PHAT are the estimates for F and p, respectively.
%
%   See also BALDNICHPDF, BALDNICHCDF, BALDNICHINV, BALDNICHSTAT, 
%            BALDNICHLIKE, BALDNICHRND
%

%   Mike Sheppard
%   Last Modified 3-Dec-2011


if nargin < 1
    error('baldnichfit:TooFewInputs',...
          'Requires at least one input argument.'); 
end

if nargin < 2 
    alpha = 0.05;
end

%Reparametrization of Beta Distribution
%USE DELTA METHOD FOR NONLINEAR REPARAMETERIZATION
phat = betafit(x(:),alpha);
phat=[1./(1+sum(phat)) phat(1)./sum(phat)];


end