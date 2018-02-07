function [m,v] = zipfstat(p,N)
%ZIPFSTAT Mean and variance for the Zipf distribution
%   [M,V] = ZIPFSTAT(P,N) returns the mean and variance of the Zipf
%   distribution with parameters P and range N.
%
%   If N is not given default is Infinity
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.


%   Mike Sheppard
%   Last Modified 24-May-2011



% 
% if nargin < 1
%     error('zipfstat:TooFewInputs',...
%           'Requires at least one input argument.'); 
% end

if nargin==1, N=Inf; end


try
    %Expand size if necessary
    p=p+zeros(size(N));
    N=N+zeros(size(p));
catch
        error('zipfstat:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize M to NaN
if isa(p,'single') || isa(N,'single')
    m = zeros(size(p),'single');
else
    m = zeros(size(p));
end
v=m; 


%Valid cases for finite case
k1=isfinite(N)&(p>0);
if any(k1)
    %For each k1, create harmonic number vector
    %(Should be able to do it without a for-loop)
    kin=find(k1);
    for i=1:length(kin)
        hmnp(i)=sum(1./((1:N(kin(i))).^p(kin(i))));
        hmn1p(i)=sum(1./((1:N(kin(i))).^(1+p(kin(i)))));
        hmnpm1(i)=sum(1./((1:N(kin(i))).^(-1+p(kin(i)))));
    end
    keyboard
    m(k1)=hmnp./hmn1p;
    v(k1)=((hmnpm1.*hmn1p)-((hmnp).^2))./(hmn1p.^2);
end

%Valid cases for infinite case
k2=(~isfinite(N))&(p>1);
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
    
    k3=(~isfinite(N))&(p<2)&(p>=1);
    v(k3)=Inf; %Correct for p<=2 var=NaN
end




end