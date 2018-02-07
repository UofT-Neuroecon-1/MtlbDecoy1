function [m,v]=vanhstat(a,pa,b,pb)
%VANHSTAT Mean and variance of the Van Houtom distribution
%   [M,V]=VANHSTAT(A,PA,B,PB) returns the mean and variance of the Van
%   Houtom distribution with parameters P(X=A)=PA, P(X=B)=PB and A<=B
%
%   The size of the output is the common size of the input arguments. 
%   A scalar input functions as a constant matrix of the same size as 
%   the other inputs.
%

%   Mike Sheppard
%   Last Modified 3-Jun-2011




if nargin < 4
    error('vanhstat:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode a pa b pb] = distchck(4,a,pa,b,pb);

if errorcode > 0
    error('vanhstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize m to NaN.
if isa(a,'single') || isa(pa,'single') || isa(b,'single') || isa(pb,'single')
    m=zeros(size(a),'single');
else
    m=zeros(size(a));
end
v=m;


%Mean
m=a.*pa+b.*pb+(1-pa-pb).*(a+b)./2;
%Variance
term1=(a.^2).*pa+(b.^2).*pb;
term2=(((a+b).*(1-pa-pb))+(2.*a.*pa)+(2.*b.*pb))/4;
term3=((b.*((2.*b)-1).*(b-1))-(a.*((2.*a)-1).*(a+1)))/6;
v=term1-term2+term3;


%Error checking
k1=(a>b)|(pa<0)|(pa>1)|(pb<0)|(pb>1);
if any(k1)
    m(k1)=NaN;
    v(k1)=NaN;
end


end