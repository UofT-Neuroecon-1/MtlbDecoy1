function y = unifdiffpdf(x)
%UNIFdiffPDF Uniform Difference probability density function
%   Y = UNIFDIFFPDF(X) returns the probability density function of the
%   Uniform Difference Distribution, evaluated at the values in X
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 1
    error('unifdiffpdf:TooFewInputs',...
          'Requires one input argument.'); 
end

%Initialize y to zeros.
if isa(x,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

k1=(x>=-1)&(x<=1);
y(k1)=1-abs(x(k1));


end