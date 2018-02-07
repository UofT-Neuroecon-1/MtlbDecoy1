function [m,v] = gbenfordstat(b,n)
%GBENFORDSTAT Mean and variance for the Generalized Benford distribution
%   [M,V] = GBENFORDSTAT(b,n) returns the mean and variance of the
%   Generalized Benford distribution of base B at the N'th digit.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Type: Discrete, bounded, {0,...,B-1}
%   Restrictions:
%        B >= 2     (B integer) [Base]
%        N >= 2     (N integer) [Digit Position]
%
%   See also GBENFORDPDF, GBENFORDCDF, GBENFORDINV, GBENFORDFIT,
%            GBENFORDLIKE, GBENFORDRND, GBENFORDSF, GBENFORDHAZ
%

%   Mike Sheppard
%   Last Modified 12-Dec-2011


if nargin < 2
    error('gbenfordstat:TooFewInputs',...
        'Requires at least two input arguments.');
end


try
    %Expand size if necessary
    b=b+zeros(size(n));
    n=n+zeros(size(b));
catch
    error('gbenfordstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


sz=size(b);
nm=numel(b);

%use gbenfordpdf to create entire pdf and then compute m and v
m=zeros(nm,1); v=zeros(nm,1); %pre-allocate
for i=1:nm
    x=0:b(i)-1;
    p=gbenfordpdf(x,b(i),n(i));
    m(i)=dot(x,p);
    v(i)=dot((x-m(i)).^2,p);
end

m=reshape(m,sz);
v=reshape(v,sz);

end