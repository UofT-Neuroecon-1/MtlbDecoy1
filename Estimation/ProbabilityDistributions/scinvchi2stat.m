function [m,v] = scinvchi2stat(nu,s)
%SCINVCHI2STAT mean and variance for the Scaled Inverse Chi Square distribution
%   [M,V] = SCINVCHI2STAT(nu,s) returns the mean and variance of the Scaled
%   Inverse Chi Square distribution with degrees of freedom NU, and scale
%   parameter S^2.
%
%   Type:
%      Continuous, Semi-bounded
%   Restrictions:
%      nu,s>0
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%

%   Mike Sheppard
%   Last Modified 12-Jun-2011

if nargin < 2
    error('scinvchi2stat:TooFewInputs',...
        'Requires two input arguments.');
end

try
    nu=nu+zeros(size(s));
    s=s+zeros(size(nu));
catch
    error('scinvchi2stat:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end

% Initialize Y to zero.
if isa(nu,'single') || isa(s,'single')
    m=zeros(size(nu),'single');
else
    m=zeros(size(nu));
end
v=m;

k2=(nu>2);
k4=(nu>4);
if any(k4)  %Variance valid for mu>4
    nuk=nu(k4); sk=s(k4);
    v(k4)=(2.*(nuk.^2).*(sk.^4))./(((nuk-2).^2).*(nuk-4));
end
if any(k2)  %Mean valid for mu>2
    m(k2)=(nu(k2).*(s(k2)).^2)./(nu(k2)-2);
end


end