function y = maxwjuttpdf(gamma,theta)
%MAXWJUTTPDF Maxwell–Juttner probability density function
%   Y = MAXWJUTTPDF(GAMMA,THETA) returns the probability density function
%   of the Maxwell-Juttner Distributionn with parameter THETA.
%
%   The Maxwell-Juttner Distribution represents the distribution for
%   speeds in a relativistic Maxwellian gas, at the values of
%   Lorentz Factor GAMMA
%
%   GAMMA is the Lorentz factor defined 1/sqrt(1-beta^2) with beta=v/c
%   THETA is defined as kT/mc2
%
%   Type: Continuous, semi-bounded
%   Restrictions:
%      GAMMA, THETA > 0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011



if nargin < 2
    error('maxwjuttpdf:TooFewInputs',...
        'Requires two input arguments.');
end

% Return NaN for out of range parameters.
gamma(gamma<=0)=NaN;
theta(theta<=0)=NaN;

try
    %Solve for beta
    bk=sqrt(1-(1./gamma).^2);
    K2 = besselk(2,1./theta);
    y=(gamma.^2.*bk./(theta.*K2)).*exp(-gamma./theta);
catch
    error('maxwjuttpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end