function y = explogcdf(x,p,b)
%EXPLOGCDF Exponential-logarithmic cumulative distribution function
%   Y = EXPLOGCDF(X,P,B) returns the cumulative distribution function
%   of the Exponential-Logarithmic Distribution with parameters P and B, 
%   at the values in X.
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Distribution: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        0 < P < 1
%        B > 0
%
%   See also EXPLOGPDF, EXPLOGINV, EXPLOGSTAT, EXPLOGFIT, 
%            EXPLOGLIKE, EXPLOGRND, EXPLOGSF, EXPLOGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin ~= 3
    error('explogcdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x p b] = distchck(3,x,p,b);

if errorcode > 0
    error('explogcdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zeros.
if isa(x,'single') || isa(p,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<p & p<1) & (0<b & b<Inf);
okvar = (0<=x & x<Inf);
ok=(okparam & okvar);
y(~okparam)=NaN;
y(okparam & x<0)=0;
y(okparam & x==Inf)=1;
  
if any(ok)
    x = x(ok); p = p(ok); b = b(ok);
    num=log(1-(1-p).*exp(-b.*x));
    den=log(p);
    y(ok)=1-(num./den);
end

%Catch round off
y(y<0)=0; y(y>1)=1;

end