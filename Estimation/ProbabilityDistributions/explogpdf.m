function y = explogpdf(x,p,b)
%EXPLOGPDF Exponential-logarithmic probability density function
%   Y = EXPLOGPDF(X,P,B) returns the probability density function of the
%   Exponential-Logarithmic Distribution with parameters P and B, 
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
%   See also EXPLOGCDF, EXPLOGINV, EXPLOGSTAT, EXPLOGFIT, 
%            EXPLOGLIKE, EXPLOGRND, EXPLOGSF, EXPLOGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin < 3
    error('explogpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode x p b] = distchck(3,x,p,b);

if errorcode > 0
    error('explogpdf:InputSizeMismatch',...
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
y(okparam & ~okvar)=0;
    
if any(ok)
    x = x(ok); p = p(ok); b = b(ok);
    expterm=(1-p).*exp(-b.*x);
    num=b.*expterm;
    den=1-expterm;
    y(ok)=(num./den).*(1./(-log(p)));
end

%Catch round off
y(y<0)=0;

end