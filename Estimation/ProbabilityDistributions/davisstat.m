function [m,v] = davisstat(b,n,mu)
%DAVISSTAT Mean and variance for the Davis distribution
%   [M,V] = DAVISSTAT(B,N,MU) returns the mean and variance of the Davis
%   distribution with scale parameter B, shape parameter N, and location
%   parameter MU.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Continuous, semi-bounded, [MU,Inf)
%   Restrictions:
%        B > 0
%        N > 1
%        MU >=0
%
%   See also DAVISPDF, DAVISCDF, DAVISINV, DAVISFIT, 
%            DAVISLIKE, DAVISRND, DAVISSF, DAVISHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011


if nargin < 3
    error('davisstat:TooFewInputs',...
          'Requires at least three input argument.'); 
end

[errorcode b n mu] = distchck(3,b,n,mu);

if errorcode > 0
    error('davisstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize M and V to zeros
if isa(b,'single') || isa(n,'single') || isa(mu,'single')
   m = zeros(size(b),'single');
else
   m = zeros(size(b));
end
v=m;

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>0 & b<Inf) & (n>1 & n<Inf) & (mu>=0 & mu<Inf);
m(~okparam)=NaN; v(~okparam)=NaN;

if any(okparam)
    b=b(okparam); n=n(okparam); mu=mu(okparam);
    %Mean
    n(n<=2)=NaN;
    m(okparam)=mu+((b.*zeta(-1+n)) ./ ((-1+n).*zeta(n)));
    %Variance
    n(n<=3)=NaN;
    num=b.^2.*(-(-2+n).*(zeta(-1+n)).^2 + (-1+n).*zeta(-2+n).*zeta(n));
    den=(-2+n).*(-1+n).^2.*(zeta(n)).^2;
    v(okparam)=num./den;
end



end
