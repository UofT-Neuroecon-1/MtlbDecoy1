function [m,v]= dncfstat(nu1,nu2,delta1,delta2)
%NCFSTAT Mean and variance for the doubly noncentral F distribution.
%   [M,V] = NCFSTAT(NU1,NU2,DELTA1,DELTA2) returns the mean and variance
%   of the doubly noncentral F pdf with NU1 and NU2 degrees of freedom and
%   noncentrality parameters, DELTA1 and DELTA2.
%
%   The Variance is not given, and will return a NaN
%

%   Mike Sheppard
%   Last Modified 11-May-2011


if nargin < 4, 
    error('dncfstat:TooFewInputs','Requires four input arguments.'); 
end

[errorcode, nu1, nu2, delta1, delta2] = distchck(4,nu1,nu2,delta1,delta2);

if errorcode > 0
    error('dncfstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize the mean and variance to NaN.
if isa(nu1,'single') || isa(nu2,'single') || isa(delta1,'single') || isa(delta2,'single')
   m = zeros(size(delta1),'single');
else
   m = zeros(size(delta1));
end
v = m;

% Return NaN for mean and variance if NU2 is less than 2.
k = find(nu2 <= 2);
m(k) = NaN;
v(k) = NaN;

% Return NaN for variance if NU2 is less than or equal to 4.
v(nu2 <= 4) = NaN;

k = find(nu2 > 2);
if any(k)
   term1=-(2.^(-2+(nu2(k)/2))).*exp(-delta2(k)/2).*nu2(k);
   term2=(1+(delta1(k)./nu1(k))).*((-delta2(k)).^(-nu2(k)/2)).*delta2(k);
   term3=gamma(-1+(nu2(k)/2))-((1-gammainc(-delta2(k)/2,-1+(nu2(k)/2))).*gamma(-1+(nu2(k)/2)));
   m(k) = term1.*term2.*term3;
   
end

%Mean should be real
k2=(imag(m)<100*eps);
if any(k2)
    m(k2)=real(m(k2));
end



end
