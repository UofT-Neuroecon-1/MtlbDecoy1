function p = moyalcdf(x,u,s)
%MOYALCDF Moyal cumulative distribution function
%   P = MOYALCDF(X,U,S) returns the Moyal cumulative distribution function
%   with location parameter U and scale parameter S.
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
    error('moyalcdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x u s] = distchck(3,x,u,s);

if errorcode > 0
    error('moyalcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(u,'single') || isa(s,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

k=(s>0);
if any(k)
    xk=x(k); uk=u(k); sk=s(k); 
    sq2=1.414213562373095048801688724209698078569671875376948073176;
    z=(xk-uk)./sk;
    term=exp(-z/2)./sq2;
    p(k)=erfc(term);
end


end