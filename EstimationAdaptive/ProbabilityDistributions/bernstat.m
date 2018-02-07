function [m,v] = bernstat(p)
% BERNSTAT Mean and variance of the Bernoulli distribution.
%   [M,V] = BERNSTAT(P) returns the mean and variance of the
%   Bernoulli distribution with parameter P.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNFIT, 
%            BERNLIKE, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011

if nargin < 1
    error('bernstat:TooFewInputs',...
          'Requires one input argument.'); 
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<=p & p<=1);
p(~okparam)=NaN;

%Mean and Variance 
m=p;
v=p.*(1-p);

end