function y = garguscdf(x,p,chi,c)
%GARGUSCDF Generalized ARGUS cumulative distribution function
%   Y = GARGUSCDF(X,P,CHI,C) returns the cumulative distribution function
%   of the Generalized ARGUS distribution with power P, curvature CHI,
%   and cut-off C, evaluated at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default value for C is 1
%
%   Distribution: Continuous, bounded, (0,C)
%   Restrictions:
%        P > -1
%        CHI, C > 0
%
%   Note: The Generalized ARGUS Distribution reduces to the ARGUS
%   distribution with p=1/2.
%
%   See also GARGUSPDF, GARGUSINV, GARGUSSTAT, GARGUSFIT,
%            GARGUSLIKE, GARGUSRND, GARGUSSF, GARGUSHAZ
%

%   Mike Sheppard
%   Last Modified 20-Dec-2011


if nargin < 3
    error('garguscdf:TooFewInputs',...
        'Requires at least three input arguments.');
end

if nargin==3
    c=1;
end

[errorcode x p chi c] = distchck(4,x,p,chi,c);

if errorcode > 0
    error('garguscdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


% Initialize Y to zero.
if isa(x,'single') || isa(p,'single') || isa(c,'single') || isa(chi,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (-1<p & p<Inf) & (0<c & c<Inf) & (0<chi & chi<Inf);
ok = (okparam & (0 < x & x < c));
y(~okparam)=NaN;
y(okparam & x<=0)=0;
y(okparam & x>=c)=1;

if any(ok)
    x=x(ok); p=p(ok); c=c(ok); chi=chi(ok);
    x=x(:); p=p(:); c=c(:); chi=chi(:); %ALL COLUMN VECTORS
    x=x./c; %CDF was defined normalized to [0,1]
    x_orig=x; p_orig=p; c_orig=c; k_orig=chi;
    k=chi; %rename for shorthand
    
    
    %STEP 1:
    %The series expansion of CDF(x;p,chi) can be simplified by first
    %multiplying by a proportionality constant that is independent of x,
    %that is, f(p,chi);
    
    term1=(gamma(1+p))-(gammainc(k.^2/2,1+p,'upper').*gamma(1+p));
    term2=((k.^2).^(1+p)).*exp(-k.^2/2);
    term3=2.^p;
    coef_factor=(term1./term2).*term3;
    
    %STEP 2:
    %Define G(x;p,chi) = CDF(x;p,chi) * coef_factor
    %The series expansion of G(x) is now semi-easily retractable into even
    %powers of x.
    gy=zeros(size(x)); goon=ones(size(x)); n=0; indx=(1:numel(gy));
    
    while any(goon)
        
        %Keep only those that need
        x=x_orig(indx); p=p_orig(indx); c=c_orig(indx); k=k_orig(indx);
        
        n=n+1; %x^(2n) term
        
        %In theory, the terms should be
        %[Coef] * x^(2n)
        %However, due to P and K there is another term
        %[Coef] * [POL(n;k,p)] * x^(2n)
        %Where POL(n;k,p) is a polynomial expansion of K and P that depends on the
        %variable n
        %All three terms are expressible:
        
        %First, the coefficient
        coef_n=1./(2^n*factorial(n));  % = 1 / (2n)!!
        
        %Second, polynomial terms of k and p, POL(n;k,p)
        %Each term of this polynomial expansion is of the form of
        %(secondary coef) * k^m * FF(p,r)
        %Where FF stands for falling factorial of p of order r
        %FF(x,n)=x*(x-1)*(x-2)*...*(x-n+1)
        %This computes all terms, as individual vectors, and then sums them
        m=n-1;
        abscoef2=round(exp(gammaln(m+1)-gammaln((0:m)+1)-gammaln(m-(0:m)+1)).*(2.^(0:m)));
        coef2=(1-2*mod(0:m,2)).*(abscoef2);
        expk=(2*m:-2:0); %exponents of k
        pterms=FF(p,[0:m])'; %Falling factorials of p
        %polyterm=sum(coef2.*(k.^expk).*pterms); %not fully vectorized
        polyterm=sum(bsxfun(@times,coef2,pterms).*bsxfun(@power,k,expk),2);
        %Example:
        %The POL term for x^14:
        %POL(7;k,p)=k^12-12*k^10*p+60*k^8*FF(p,2)-160*k^6*FF(p,3)+240*k^4*FF(p,4)-192*k^2*FF(p,5)+64*FF(p,6);
        %The coefficients are from a triangle whose (i, j)-th entry is binomial(i, j)*2^j.
        
        %Third, the x^(2n) term
        xp2n=x.^(2*n);
        
        %Multiply everything together to get the new term in the sequence
        term_n=coef_n.*polyterm.*xp2n;
        
        %Check to see which ones were good
        goon=(abs(term_n)>eps) & isfinite(term_n);
        indx(~goon)=[];
        
        %Add them to previous terms
        if ~isempty(indx)
            gy(indx)=gy(indx)+term_n(goon);
        end
        
    end
    
    
    %STEP 3:
    %Once the series expression for G(x;p,chi) has been computed, divide by
    %coef_factor to return back to CDF via
    %G(x;p,chi) = CDF(x;p,chi) * coef_factor
    
    gy=gy./coef_factor;
    y(ok)=gy;
    
end

%Catch round off
y(y<0)=0; y(y>1)=1;


    function y=FF(x,n)
        %falling factorial
        %N is the number of terms
        %x*(x-1)*(x-2)*...*(x-n+1)
        for kk=1:numel(n)
            y(kk,:)=gamma(x+1)./gamma(x-n(kk)+1);
        end
    end
end