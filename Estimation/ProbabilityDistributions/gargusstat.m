function [m,v] = gargusstat(p,chi,c)
%ARGUSSTAT Mean and variance for the Generalized ARGUS distribution
%   [M,V] = GARGUSSTAT(P,CHI,C) returns the mean and variance for the
%   Generalized ARGUS distribution with power P, curvature CHI,
%   and cut-off C.
%
%   The sizes of M and V are the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as
%   the other inputs.
%
%   Default value for C is 1.
%
%   Type: Continuous, bounded, (0,C)
%   Restrictions:
%        P > -1
%        CHI, C > 0
%
%   See also GARGUSPDF, GARGUSCDF, GARGUSINV, GARGUSFIT,
%            GARGUSLIKE, GARGUSRND, GARGUSSF, GARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 3
    error('gargusstat:TooFewInputs',...
        'Requires at least one input argument.');
end

if nargin == 2
    c=1;
end



[errorcode p chi c] = distchck(3,p,chi,c);

if errorcode > 0
    error('gargusstat:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize M to zero.
if isa(p,'single') || isa(c,'single') || isa(chi,'single')
    m = zeros(size(p),'single');
else
    m = zeros(size(p));
end
v=m;

% Return NaN if any arguments are outside of their respective limits.
okparam = (-1<p & p<Inf) & (0<c & c<Inf) & (0<chi & chi<Inf);
m(~okparam)=NaN; v(~okparam)=NaN;

if any(okparam)
    p=p(okparam); chi=chi(okparam); c=c(okparam);
    
    %First moment
    %PFQ is not vectorized in parameter space, that is it does not know if
    %the parameter P is a vector parameter, or a vector of singular
    %parameters
    H1F1Reg=zeros(size(p));
    %Create unique matrix, to speed of operations
    %If multiple [p,chi] pairs, do one calculation each
    [pchiM,ii,jj]=unique([p(:) chi(:)],'rows');
    for kk=1:numel(ii)
        indx=(jj==kk); p_t=pchiM(kk,1); chi_t=pchiM(kk,2);
        H1F1Reg(indx)=pfq(1+p_t,(5/2)+p_t,-chi_t.^2/2)./gamma((5/2)+p_t);
    end
    
    term1=(2.^(-2-p)).*c.*((chi.^2).^(1+p)).*sqrt(pi).*gamma(1+p);
    num1m=term1.*H1F1Reg;
    den1m=(gamma(1+p))-(gammainc(chi.^2/2,1+p,'upper').*gamma(1+p));
    M1=num1m./den1m;
    
    %Second moment
    term2=(gamma(2+p))-(gammainc(chi.^2/2,2+p,'upper').*gamma(2+p));
    term3=2.*((chi.^2)-2.*(1+p));
    term4=(2.^(-p)).*(chi.^4).*((chi.^2).^p).*exp(-chi.^2/2);
    term5=(c.^2).*gamma(1+p);
    term6=2.*chi.^2.*gamma(2+p);
    num2m=term5.*(term4+term3.*term2);
    den2m=term6.*den1m;
    M2=num2m./den2m;
    
    %Substitute
    m(okparam)=M1;
    v(okparam)=M2-(M1.^2);
end


end

