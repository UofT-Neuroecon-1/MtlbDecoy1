function y = bernpdf(x,p)
% BERNPDF Bernoulli probability density function
%   Y = BERNPDF(X,P) returns the probability density function of the
%   Bernoulli Distribution with parameter P at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNCDF, BERNINV, BERNSTAT, BERNFIT,
%            BERNLIKE, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 13-May-2012


if nargin ~=2
    error('bernpdf:TooFewInputs',...
        'Requires two input arguments.');
end

try
    %Match dimensions
    x=x+zeros(size(p));
    p=p+zeros(size(x));
catch err
    error('bernpdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<=p & p<=1);
okvar = (x==0 | x==1);
ok = (okparam & okvar);
y(okparam & ~okvar)=0;
y(~okparam)=NaN;

if any(ok)
    x=x(ok); p=p(ok);
    y(ok)=x.*p+(1-x).*(1-p); %Catch both cases at once
end

%Catch round off
y(y<0)=0;


end