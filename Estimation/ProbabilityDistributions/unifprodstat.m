function [m,v]=unifprodstat(n)
%UNIFPRODSTAT Mean and variance of the Uniform Product Distribution
%   [M,V]=UNIFPRODSTAT(N) returns the mean and variance of the Uniform
%   Product Distribution of N uniform distributions
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.

%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('uquadstat:TooFewInputs',...
          'Requires one input arguments.'); 
end


%Initialize M to NaN.
if isa(n,'single')
    m=zeros(size(n),'single');
else
    m=zeros(size(n));
end
v=m;


try
    m=2.^(-n);
    v=(3.^(-n))-(m.^2);
catch
    error('unifprodstat:InputSizeMismatch');
end

% Return NaN for out of range parameters.
k1=(n~=round(n));
m(k1)=NaN; v(k1)=NaN;


end