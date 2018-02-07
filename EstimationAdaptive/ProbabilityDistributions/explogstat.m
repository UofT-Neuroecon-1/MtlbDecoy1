function [m,v] = explogstat(p,b)
%ERLANGSTAT Mean and variance for the Exponential-logarithmic distribution
%   [M,V] = EXPLOGSTAT(P,B) returns the mean and variance of the 
%   Exponential-Logarithmic Distribution with parameters P and B.
%
%   The sizes of M and V are the common size of the input arguments.  
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.   
%
%   Type: Continuous, semi-bounded, [0,Inf)
%   Restrictions:
%        0 < P <1
%        B > 0
%
%   See also EXPLOGPDF, EXPLOGCDF, EXPLOGINV, EXPLOGFIT, 
%            EXPLOGLIKE, EXPLOGRND, EXPLOGSF, EXPLOGHAZ
%

%   Mike Sheppard
%   Last Modified 18-Dec-2011

if nargin < 2
    error('explogstat:TooFewInputs',...
          'Requires at least two input argument.');
end


try
    % Return NaN if any arguments are outside of their respective limits.
    p(~(0<p & p<1))=NaN;
    b(~(0<b & b<Inf))=NaN;
    
    pol2=polylog(2,1-p);
    pol3=polylog(3,1-p);
    m =-pol2./(b.*log(p));
    v = (-2*pol3./(b.^2.*log(p)))-(m.^2);
catch
    error('explogstat:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end


end

function a=polylog(n,z)
%Good for when 0<z<1, for which this is the case
%Scalar n, and vector z
a=zeros(size(z)); k=1; term=z;
while max(term(:))>(eps)
    term=(z.^k)./(k.^n);
    a=a+term;
    k=k+1;
end

end
