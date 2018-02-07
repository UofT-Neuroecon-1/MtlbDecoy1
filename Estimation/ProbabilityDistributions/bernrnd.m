function r = bernrnd(p,varargin)
% BERNRND Random arrays from the Bernoulli distribution.
%   R = BERNRND(P) returns an array of random numbers chosen from a
%   Bernoulli distribution with parameter P. 
%
%   The size of R is the common size of the parameters if all are arrays.  
%   If any parameter is a scalar, the size of R is the size of the other
%   parameters.
%
%   R = BERNRND(P,M,N,...) or R = BERNRND(P,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNSTAT,
%            BERNFIT, BERNLIKE, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 4-Dec-2011

if isempty(varargin), varargin={1}; end %Scalar

%Return NaN for out of range parameters.
p(p<0 | p>1)=NaN;

try
    r = (rand(varargin{:})<=p);
    if isnan(p), r=NaN(size(r)); end
catch err
    error('bernrnd:InputSizeMismatch',...
         'Size information is inconsistent.');
end

end
