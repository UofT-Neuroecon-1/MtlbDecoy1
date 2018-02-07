function logL = arguslike(params,data)
%ARGUSLIKE ARGUS distribution log-likelihood function.
%   LOGL = ARGUSLIKE(PARAMS,DATA) returns the negative of the ARGUS
%   distribution log-likelihood function for the parameters
%   PARAMS(1) = CHI, and PARAMS(2) = C.
%
%   Type: Continuous, bounded
%   Restrictions:
%        0 <= DATA <= C
%        CHI > 0
%
%   See also ARGUSPDF, ARGUSCDF, ARGUSINV, ARGUSSTAT, ARGUSFIT, ARGUSRND
%

%   Mike Sheppard
%   Last Modified 26-Apr-2011

if nargin < 2
    error('arguslike:TooFewInputs','Requires at least two input arguments.');
end

if isscalar(params)
    params=[params 1];
end

chi=params(1); c=params(2); x=data(:);

%Return NaN for out of range parameters or data.
chi(chi<0)=NaN; c(c<0)=NaN; x(x<0|x>c)=NaN;

repterm=(1-(x./c).^2);
logterm1=3.*log(chi)-(1/2).*log(2*pi)-logpsi(chi);
logterm2=log(x)-2*log(c)+(1/2).*log(repterm);
logterm3=(-1/2).*chi.^2.*repterm;
logL=-sum(logterm1+logterm2+logterm3);

end


function logy_psi=logpsi(x)
logy_psi=log(normcdf(x)-x.*normpdf(x)-.5);
end
