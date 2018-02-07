function [mn,v]= stzstat(nu);
%TSTAT  Mean and variance for the Student's Z distribution.
%   [MN,V] = TSTAT(NU) returns the mean and variance of
%   Student's Z distribution with NU degrees of freedom.
%

%   Mike Sheppard
%   Last Modified 6-Jun-2011


if nargin < 1,   
    error('stzstat:TooFewInputs');       
end

% Initialize the mean and variance to zero.
if isa(nu,'single')
   mn = zeros(size(nu),'single');
else
   mn = zeros(size(nu));
end

k = find(nu <= 0);
if any(k)
    mn(k) = NaN;
end
v = mn;

% The mean of the z distribution is zero unless there is only one 
% degree of freedom. In that case the mean does not exist.
mn(nu <= 1) = NaN;

% The variance of the z distribution is undefined for one and two
% or three degrees of freedom.
v(nu <= 3) = NaN;

k = (nu > 3);
v(k)=1./(nu(k)-3);
end


