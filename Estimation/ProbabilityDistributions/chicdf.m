function p = chicdf(x,v)
%CHICDF Chi cumulative distribution function.
%   P = CHICDF(X,V) returns the cumulative distribution function of the
%   Chi Distribution with V degrees of freedom, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        V > 0
%
%   See also CHIPDF, CHIINV, CHISTAT, CHIFIT,
%            CHILIKE, CHIRND, CHISF, CHIHAZ
%

%   Mike Sheppard
%   Last Modified 17-Dec-2011

if   nargin ~= 2,
    error('chicdf:TooFewInputs','Requires two input arguments.');
end


try
    x=x+zeros(size(v));
    v=v+zeros(size(x));
catch err
    error('chicdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<v & v<Inf);
okvar = (0 <= x) & (x < Inf);
ok=(okparam & okvar);
p(~okparam)=NaN;
p(okparam & x<0)=0;
p(okparam & x==Inf)=1;

if any(ok)
    x = x(ok); v = v(ok);
    p(ok) = gammainc(x.^2/2,v/2);
end

%Catch round off
p(p<0)=0; p(p>1)=1;

end
