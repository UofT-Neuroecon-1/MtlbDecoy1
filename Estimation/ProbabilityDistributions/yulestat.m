function [m,v] = yulestat(a)
%YULESTAT Mean and variance of the Yule–Simon distribution
%   [M,V] = YULESTAT(A) returns the mean and variance of the Yule-Simon
%   distribution with scale parameter A.
%
%   Note: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1
%
%   The sizes of M and V are the common size of the input arguments.  A
%   scalar input functions as a constant matrix of the same size as the
%   other inputs.
%

%   Mike Sheppard
%   Last Modified 6-Apr-2011


if nargin < 1
    error('yulestat:TooFewInputs',...
          'Requires at least one input argument.');
end

%Initialize
m=zeros(size(a));
v=m;

k1=a>1; k2=a>2;
if any(k1 | k2) 
    try
m(k1)=a(k1)./(a(k1)-1);
v(k2)=a(k2).^2 ./ ((a(k2)-1).^2.*(a(k2)-2));
    catch
     error('yulestat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
    end
end

m=m-1;  %Support on x>=1
end