function x = moyalinv(p,u,s)
%MOYALINV Inverse of the Moyal cumulative distribution function
%   X = MOYALCDF(P,U,S) returns the inverse Moyal cumulative distribution
%   function with location parameter U and scale parameter S.
%
%   Type: Continuous, unbounded
%   Restrictions:
%      U any real number
%      S>0
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 20-Jun-2011



if nargin < 3
    error('moyalinv:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode p u s] = distchck(3,p,u,s);

if errorcode > 0
    error('moyalinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize X to zero.
if isa(p,'single') || isa(u,'single') || isa(s,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k=(s>0 & p>0 & p<1);
if any(k)
    pk=p(k); uk=u(k); sk=s(k); 
    sq2=1.414213562373095048801688724209698078569671875376948073176;
    term=erfcinv(pk);
    z=-2.*log(term.*sq2);
    x(k)=uk+(sk.*z);
end

%End cases
x(s>0 & p==0)=-Inf;
x(s>0 & p==1)=Inf;


end