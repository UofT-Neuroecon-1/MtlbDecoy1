function y = berncdf(x,p)
%BERNCDF Bernoulli cumulative distribution function.
%   Y=BERNCDF(X,P) returns the cumulative distribution function of the
%   Bernoulli distribution with parameter P at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%        0 <= P <= 1
%
%   See also BERNPDF, BERNINV, BERNSTAT, BERNFIT,
%            BERNLIKE, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified 16-Dec-2011


if nargin < 2
    error('berncdf:TooFewInputs',...
        'Requires two input arguments.');
end



try
    %Match dimensions
    x=x+zeros(size(p));
    p=p+zeros(size(x));
catch err
    error('berncdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<=p & p<=1);
okvar = (x==0 | x==1);
ok=(okparam & okvar);

y=zeros(size(x));
y(ok)=x(ok).*1+(1-x(ok)).*(1-p(ok)); %Catch both cases at once

y(~okparam)=NaN;
y(okparam & x<0)=0;
y(okparam & x>1)=1;

%Catch round off
y(y<0)=0; y(y>1)=1;


end