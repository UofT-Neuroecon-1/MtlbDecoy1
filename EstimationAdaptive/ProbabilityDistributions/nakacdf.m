function p = nakacdf(x,mu,omega)
%NAKACDF Nakagami cumulative distribution function
%   P = NAKACDF(X,MU,OMEGA) returns the cumulative distribution function 
%   of the Nakagami Distribution with shape parameter MU and scale 
%   parameter OMEGA, at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Distribution: Continuous, semi-bounded, (0, Inf)
%   Restrictions:
%      MU, OMEGA > 0
%
%   See also NAKAPDF, NAKAINV, NAKASTAT, NAKAFIT,
%            NAKALIKE, NAKARND, NAKASF, NAKAHAZ
%

%   Mike Sheppard
%   Last Modified 10-Dec-2011



if nargin ~= 3
    error('nakacdf:TooFewInputs',...
          'Requires three input arguments.');
end

[errorcode x mu omega] = distchck(3,x,mu,omega);

if errorcode > 0
    error('nakacdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(mu,'single') || isa(omega,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end


% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

%From addnaka.m in statistics toolbox
%-----
x(x<0) = 0;
p = gamcdf(x.^2, mu, omega./mu);
%-----


%Catch round off
p(p<0)=0; p(p>1)=1;

end