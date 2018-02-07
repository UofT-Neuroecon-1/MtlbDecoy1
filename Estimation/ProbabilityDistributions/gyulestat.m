function [m,v] = gyulestat(rho,alpha)
%GYULESTAT Mean and variance of the Generalized Yule–Simon distribution
%   [M,V] = GYULESTAT(RHO,ALPHA) returns the mean and variance of the Yule-Simon
%   distribution with scale parameter A.
%
%   Note: Some define the support of Yule-Simon to be x>=0
%   This function uses x>=1
%
%   The sizes of M and V are the common size of the input arguments.  A
%   scalar input functions as a constant matrix of the same size as the
%   other inputs.
%

% Superstars without talent?
% The Yule distribution controversy
% Laura Spierdijka and Mark Voorneveldb, c,1
% a Dept. of Econometrics, University of Groningen, The Netherlands
% b Dept. of Econometrics and Operations Research, Tilburg University, The Netherlands
% c Dept. of Economics, Stockholm School of Economics, Sweden

%   Mike Sheppard
%   Last Modified 26-May-2011


if nargin < 2
    error('gyulestat:TooFewInputs',...
        'Requires at least one input argument.');
end

try
    %Expand size if necessary
    rho=rho+zeros(size(alpha));
    alpha=alpha+zeros(size(rho));
catch
    error('gyulestat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


%Initialize
if isa(rho,'single') || isa(alpha,'single')
    m=zeros(size(rho),'single');
else
    m=zeros(size(rho));
end
v=m;
sm=m;


if any(rho==1)
    k1=(rho==1);
    alp=alpha(k1);
    m(k1)=-log(alp)./(1-alp);
    sm(k1)=(2*(1-alp)+alp.*log(alp))./(alp.*(1-alp));
elseif any(rho==2)
    k2=(rho==2);
    alp=alpha(k2);
    m(k2)=-log(alp)./(1-alp);
    sm(k2)=2*(alp-1-2*log(alp))./(1-(alp.^2));
else
    kn=(rho~=1)&(rho~=2);
    alp=alpha(kn);
    rhon=rho(kn);
    m(kn)=rhon.*(1-(alp.^(rhon-1)))./((rhon-1).*(1-(alp.^rhon)));
    term1=rhon./(1-(alp.^rhon));
    num=(-2*(alp.^(rhon-2)).*(rhon-1))+((alp.^(rhon-1)).*(rhon-2))+rhon;
    den=(rhon-1).*(rhon-2);
    sm(kn)=term1.*(num./den);
end

v=sm-m.^2;

end
