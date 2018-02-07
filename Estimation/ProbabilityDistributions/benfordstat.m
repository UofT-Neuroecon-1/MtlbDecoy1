function [m,v] = benfordstat(b)
%BENFORDSTAT Mean and variance for the Benford distribution
%   [M,V] = BENFORDSTAT(b) returns the mean and variance of the Benford
%   distribution of base B. Default value for B is 10. (Base 10)
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Discrete, bounded, {1,...,B-1}
%   Restrictions:
%        B >= 2    (B integer)  [Base]
%
%   See also BENFORDPDF, BENFORDCDF, BENFORDINV, BENFORDFIT,
%            BENFORDLIKE, BENFORDRND, BENFORDSF, BENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011


if nargin < 1
    b=10;
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (b>=2 & b<Inf & b==round(b));
b(~okparam)=NaN;

%Initialize M and V to zero.
if isa(b,'single')
    m=zeros(size(b),'single');
else
    m=zeros(size(b));
end
v=m;


%Although there is exact formulas for the mean and variance, the variance
%computations takes just as long as computing both by the definition

%Vectorized expression for mean:
%m=b-(gammaln(b+1)./log(b)); %vectorized, and using gammaln

for k=1:numel(b)
    x=1:(b(k)-1);
    y=benfordpdf(x,b(k));
    m(k)=sum(x.*y);
    v(k)=sum(x.^2.*y)-(m(k)).^2;
end


end