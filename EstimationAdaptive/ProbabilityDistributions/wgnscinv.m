function x = wgnscinv(p,R)
%WGNSCINV Inverse of the Wigner semicircle cumulative distribution function
%   X = WGNSCINV(P,R) returns the inverse cumulative distribution function 
%   of the Wigner semicircle distribution with radius R, evaluated at the values in P.

%   Mike Sheppard
%   Last Modified 22-Apr-2011


if nargin < 2
    error('wgnscinv:TooFewInputs',...
          'Requires at least two input arguments.'); 
end


if isa(p,'single') ||  isa(R,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


try
  x = R.*( -1 + 2.*betainv(p,3/2,3/2) );
catch
  error('wgnscinv:InputSizeMismatch',...
         'Non-scalar arguments must match in size.');
end

k = (p<0 | p>1);
x(k)=NaN;

end
