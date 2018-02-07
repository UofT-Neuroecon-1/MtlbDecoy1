function [phat, pci] = bernfit(x,alpha)
%BERNFIT Parameter estimates and confidence intervals for Bernoulli data
%   PHAT = BERNFIT(X) returns estimates of the probability of success for
%   the Bernoulli distribution. X is a vector of indicator variables of 0
%   or 1 indicating failure or success.
%
%   [PHAT, PCI] = BERNFIT(X,ALPHA) gives MLEs and 100(1-ALPHA) percent
%   confidence intervals given the data. By default, the optional parameter
%   ALPHA=0.05 corresponding to 95% confidence intervals.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNSTAT,
%            BERNLIKE, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 14-May-2012

if nargin < 2
    alpha = 0.05;
end

%As Benoulli is built in to MLE, let it do all the work
try
    [phat, pci] = mle(x,'distribution','Bernoulli','alpha',alpha);
catch err
    error('bernmle:InvalidBernoulliData',...
        'Bernoulli data must either be 1 or 0');
end

end