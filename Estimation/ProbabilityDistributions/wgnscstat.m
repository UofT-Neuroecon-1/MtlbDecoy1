function [m,v] = wgnscstat(R)
%WGNSCSTAT Mean and variance of the Wigner semicircle distribution
%   [M,V] = wgnscstat(R) returns the mean and variance of the Wigner
%   semicircle distribution with radius R.
%

if nargin < 1
    error('wgnscstat:TooFewInputs',...
          'Requires at least one input argument.');  
end


if isa(R,'single')
    m=zeros(size(R),'single');
else
    m=zeros(size(R));
end

v=(R.^2)/4;


end
