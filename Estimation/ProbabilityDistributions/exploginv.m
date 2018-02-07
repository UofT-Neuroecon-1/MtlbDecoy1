function x = exploginv(y,p,b)
%EXPLOGINV Inverse of the Exponential-logarithmic cdf
%   X = EXPLOGINV(Y,P,B) returns the inverse cdf of the 
%   Exponential-logarithmic distribution with parameters P and B, 
%   at the values in Y.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%      0 < P < 1
%      B > 0
%
%   See also EXPLOGPDF, EXPLOGCDF, EXPLOGSTAT, EXPLOGFIT,
%            EXPLOGLIKE, EXPLOGRND, EXPLOGSF, EXPLOGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 3
    error('exploginv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode y p b] = distchck(3,y,p,b);

if errorcode > 0
    error('exploginv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(y,'single') || isa(p,'single') || isa(b,'single')
    x=zeros(size(y),'single');
else
    x=zeros(size(y));
end


% Return NaN if any arguments are outside of their respective limits.
okparam = (0<p & p<1) & (0<b & b<Inf);
okvar = (0<y & y<1);
ok=(okparam & okvar);
x(~okparam)=NaN;
x(okparam & y==0)=0;
x(okparam & y==1)=Inf;

if any(ok)
    y=y(ok); p=p(ok); b=b(ok);
    x(ok) = (log(-((p-1).*p.^y)./(p.^y-p)))./b;
end


end