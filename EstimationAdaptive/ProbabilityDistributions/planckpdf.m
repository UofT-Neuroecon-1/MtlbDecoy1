function y = planckpdf(x)
%PLANCKPDF Planck's Radiation probability density function
%   Y = PLANCKPDF(X) returns the probability density function of
%   Planck's Radiation Function at the values in X.
%
%   The Planck's Radiation Function is given by
%   f(x) = (15/pi^4) * 1/((x^5)(exp(1/x)-1))
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%      x >= 0
%
%   The size of Y is the size of the input variable X.
%

%   Mike Sheppard
%   Last Modified 16-Jun-2011


if nargin < 1
    error('planckpdf:TooFewInputs',...
          'Requires one input argument.');
end


% Return NaN for out of range parameters.
x(x<0)=NaN;

%Planck's Function
y = (15./pi.^4) .* (1./((x.^5).*(exp(1./x)-1)));

%Round off
y(y<0)=0;

end