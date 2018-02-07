function h = bernhaz(x,p)
%BERNHAZ Bernoulli hazard function.
%   H=BERNHAZ(X,P) returns the Bernoulli hazard function with 
%   parameter P at the values in X.
%
%   The size of H is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNINV, BERNSTAT, 
%            BERNFIT, BERNLIKE, BERNRND, BERNSF
%

%   Mike Sheppard
%   Last Modified 21-Dec-2011


if nargin < 2
    error('bernhaz:TooFewInputs',...
          'Requires two input arguments.'); 
end



try
    yt = bernpdf(x,p);
    st = bernsf(x,p);
    h = yt ./ (yt+st);  % +yt term in denominator for discrete r.v.
catch err
    error('bernhaz:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


end