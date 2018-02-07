function [m,v] = zetastat(p)
%ZETASTAT Mean and variance for the Zeta distribution
%   [M,V] = ZETASTAT(P,N) returns the mean and variance of the Zeta
%   distribution with parameters P.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 25-May-2011



if nargin < 1
    error('zetastat:TooFewInputs',...
          'Requires at least one input argument.'); 
end



% Initialize M to NaN
if isa(p,'single')
    m = zeros(size(p),'single');
else
    m = zeros(size(p));
end
v=m; 


%Valid cases for infinite case
k2=(p>1);
if any(k2)
    %For each k1, create harmonic number vector
    %(Should be able to do it without a for-loop)
    kin=find(k2);
    %Use zeta function
    zp=zeta(p(k2));
    zpp1=zeta(1+p(k2));
    zpm1=zeta(-1+p(k2));
    m(k2)=zp./zpp1;
    v(k2)=(zpm1./zpp1)-(m(k2)).^2;
    
    k3=(p<2)&(p>=1);
    v(k3)=Inf; %Correct for p<=2 var=NaN
end




end