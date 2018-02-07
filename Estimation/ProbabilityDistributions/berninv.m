function x = berninv(y,p)
%BERNINV Inverse of the Bernoulli cumulative distribution function
%   X = BERNINV(Y,P) returns the inverse cumulative distribution function
%   of the Bernoulli distribution with parameter P.
%
%   Since the Bernoulli distribution is discrete, BERNINV returns the
%   least integer X such that the Bernoulli cdf evaluated at X, equals
%   or exceeds Y.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Discrete, bounded, {0,1}
%   Restrictions:
%      0 <= P <= 1
%
%   See also BERNPDF, BERNCDF, BERNSTAT, BERNFIT,
%            BERNLIKE, BERNRND, BERNSF, BERNHAZ
%

%   Mike Sheppard
%   Last Modified: 14-May-2012

if nargin ~= 2
    error('berninv:TooFewInputs',...
        'Requires two input arguments.');
end


try
    %Match dimensions
    y=y+zeros(size(p));
    p=p+zeros(size(y));
catch err
    error('berninv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(y,'single') || isa(p,'single')
    x = zeros(size(y),'single');
else
    x = zeros(size(y));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<=p & p<=1);
okvar = (0<=y & y<=1);
ok=(okparam & okvar);
x(~ok)=NaN;

if any(ok)
    y=y(ok); p=p(ok);
    x(ok)=(0)*(y<=(1-p))+(1)*(y(k)>(1-p)); %Catch both cases at once
end



end