function y = davispdf(x,b,n,mu)
%DAVISPDF Davis probability density function
%   Y = DAVISPDF(X,B,N,MU) returns the probability density function
%   of the Davis distribution with scale parameter B, shape parameter N, 
%   and location parameter MU, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Distribution: Continuous, semi-bounded, [MU,Inf)
%   Restrictions:
%        B > 0
%        N > 1
%        MU >=0
%
%   See also DAVISCDF, DAVISINV, DAVISSTAT, DAVISFIT, 
%            DAVISLIKE, DAVISRND, DAVISSF, DAVISHAZ
%

%
%  Note: Requires "Special Functions math library"
%  http://www.mathworks.com/matlabcentral/fileexchange/978
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 4
    error('davispdf:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode x b n mu] = distchck(4,x,b,n,mu);

if errorcode > 0
    error('davispdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(b,'single') || isa(n,'single') || isa(mu,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>0 & b<Inf) & (n>1 & n<Inf) & (mu>=0 & mu<Inf);
okvar = (mu <= x & x < Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & ~okvar)=0;

if any(ok)
    x = x(ok); b = b(ok); n = n(ok); mu = mu(ok);
    num=(b.^n).*(x-mu).^(-1-n);
    denexp=-1+exp(b./(x-mu));
    den=denexp.*gamma(n).*zeta(n);
    y(ok) = num./den;
end

%Catch round off
y(y<0)=0;

end