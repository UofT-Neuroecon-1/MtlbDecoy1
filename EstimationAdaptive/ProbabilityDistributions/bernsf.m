function s = bernsf(x,p)
%BERNSF Bernoulli survival function.
%   S=BERNSF(X,P) returns the survival function of the Bernoulli 
%   distribution with parameter P at the values in X.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNSTAT, 
%            BERNFIT, BERNLIKE, BERNRND, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 2
    error('bernsf:TooFewInputs',...
          'Requires two input arguments.'); 
end


try
    s = 1 - berncdf(x,p);
catch err
    error('bernsf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end