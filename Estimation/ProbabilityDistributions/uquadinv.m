function x = uquadinv(p,a,b)
%UQUADINV Inverse of the U-quadratic cumulative distribution function
%   X = UQUADINV(P,A,B) returns the inverse of the U-quadratice cumulative
%   density function with lower limit A and upper limit B, evaluated at the
%   values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 4-Jun-2011


if nargin < 3
    error('uquadinv:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('uquadinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


%Initialize X to 0.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
else
    x=zeros(size(p));
end


% Return NaN for out of range parameters.
a(a>b) = NaN;
b(a>b) = NaN;
p(p<0 | p>1)=NaN;


try
   alpha=12./((b-a).^3); %vertical scale
   beta=(b+a)./2; %gravitational balance center, offset
   temp=((3.*p./alpha)-((beta-a).^3));
   %Split between below or above beta
   kn=(temp<0); kp=(temp>0);
   if any(kn)
       tn=temp(kn); bn=beta(kn);
       x(kn)=bn-(-tn).^(1/3);
   end
   if any(kp)
       tp=temp(kp); bp=beta(kp);
       x(kp)=bp+(tp).^(1/3);
   end
catch
    error('uquadinv:InputSizeMismatch');
end



end