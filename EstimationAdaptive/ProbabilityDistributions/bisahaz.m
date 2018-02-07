function h = bisahaz(x,gamma,beta,mu)
%BISAHAZ Birnbaum-Saunders hazard function
%   H=BISAHAZ(X,GAMMA,BETA,MU) returns the hazard function of the
%   Birnbaum-Saunders Distribution with shape parameter GAMMA,
%   scale parameter BETA, and location parameter MU, evaluated at the
%   values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for BETA and MU are 1 and 0, respectively.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%        GAMMA, BETA > 0
%
%   Note: The Birnbaum-Saunders Distribution is also known as the Fatigue
%   Life Distribution. Alternative definitions have BETA defined as
%   1/LAMBDA with LAMBDA being the scale parameter. The definition used
%   here is the notation (X-MU)/BETA
%
%   See also BISAPDF, BISACDF, BISAINV, BISASTAT, BISAFIT, BISALIKE, 
%            BISARND, BISASF
%

%   Mike Sheppard
%   Last Modified: 20-Dec-2011


if nargin < 2
    error('bisahaz:TooFewInputs',...
        'Requires at least two input arguments.');
end
if nargin < 3, beta = 1; end
if nargin < 4, mu = 0; end


try
    h = bisapdf(x,gamma,beta,mu) ./ bisasf(x,gamma,beta,mu);
catch
    error('bisahaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end
