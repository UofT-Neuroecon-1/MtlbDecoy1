function x = chiinv(p,v)
%CHIINV Inverse of the Chi cumulative distribution function (cdf).
%   X = CHIINV(P,V)  returns the inverse cumulative distribution function
%   of the Chi distribution with V degrees of freedom, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%      V > 0
%
%   See also CHIPDF, CHICDF, CHISTAT, CHIFIT,
%            CHILIKE, CHIRND, CHISF, CHIHAZ
%

%   Mike Sheppard
%   Last Modified 16-May-2012

if nargin ~= 2,
    error('chiinv:TooFewInputs','Requires two input arguments.');
end


try
    p=p+zeros(size(v)); %Match dimensionality
    v=v+zeros(size(p));
catch err
    error('chiinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(v,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<v & v<Inf);
okvar = (0 < p) & (p < 1);
ok=(okparam & okvar);
x(~ok)=NaN; 
x(okparam & p==0)=0;
x(okparam & p==1)=Inf;

if any(ok)
    p=p(ok); v=v(ok);
    %Call the Inverse incomplete gamma function.
    x(ok) = sqrt(2*gammaincinv(p,v/2));
end


end