function [m,v]=lindstat(s)
%LINDSTAT Mean and variance for the Lindley Distribution
%   [M,V] = LINDSTAT(s) returns the mean and variance of the Lindley
%   distribution with shape parameter S.
%

%   Mike Sheppard
%   Last Modified 20-May-2011


if nargin < 1
    error('lindstat:TooFewInputs',...
          'Requires at least one input arguments.'); 
end


m=(2+s)./(s.*(1+s));
v=(2./(s.^2))-(1./((1+s).^2));

m(s<0)=NaN;
v(s<0)=NaN;


end