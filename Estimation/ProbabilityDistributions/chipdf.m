function y = chipdf(x,v)
%CHIPDF Chi probability density function
%   Y = CHIPDF(X,V) returns the probability density function of the
%   Chi Distribution with V degrees of freedom, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        V > 0
%
%   See also CHICDF, CHIINV, CHISTAT, CHIFIT,
%            CHILIKE, CHIRND, CHISF, CHIHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011

if nargin ~= 2,
    error('chipdf:TooFewInputs','Requires two input arguments.');
end

try
    %Match dimensions
    x=x+zeros(size(v));
    v=v+zeros(size(x));
catch err
    error('chipdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<v & v<Inf);
okvar = (0 <= x) & (x < Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x=x(ok); v=v(ok);
    y(ok)=pow2(1-(v/2)).*(x.^(v-1)).*exp(-(x.^2)/2)./gamma(v/2);
end


%Catch round off
y(y<0)=0;

end

