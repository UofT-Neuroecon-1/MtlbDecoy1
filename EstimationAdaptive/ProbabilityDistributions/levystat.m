function [m,v] = levystat(mu,sigma)
%LEVYSTAT Mean and variance of the Lévy Distribution
%   [M,V] = LEVYSTAT(MU,SIGMA) returns the mean and variance of the Lévy
%   Distribution with location parameter MU and dispersion parameter SIGMA. 
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, (MU, Inf)
%   Restrictions:
%     SIGMA > 0
%
%   NOTE: The mean and variance is INF for all inputs
%
%   See also LEVYPDF, LEVYCDF, LEVYINV, LEVYFIT, 
%            LEVYLIKE, LEVYRND, LEVYSF, LEVYHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin ~= 2
    error('levystat:TooFewInputs',...
          'Requires two input arguments.'); 
end


try
    bothinf=inf(size(mu))+inf(size(sigma)); %check dimensionality
    m = bothinf; v=bothinf;
catch
    error('levystat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end

