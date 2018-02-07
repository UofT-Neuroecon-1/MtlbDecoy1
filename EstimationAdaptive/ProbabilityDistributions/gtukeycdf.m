function p = gtukeycdf(x,lambda1,lambda2,mu,sigma1,sigma2)
%GTUKEYCDF Cumulative function of the Generalized Tukey Lambda distribution
%   P = GTUKEYCDF(X,LAMBDA1, LAMBDA2, MU, SIGMA1, SIGMA2) returns the CDF 
%   of the Generalized Tukey Lambda distribution with location
%   parameter MU, scale parameters SIGMA1 and SIGMA2, and shape parameters
%   LAMBDA1, LAMBDA2
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   GTUKEYCDF uses Newton's method to converge to the solution, as there is
%   no closed form expression other than for GTUKEYINV.


%   Mike Sheppard
%   Last Modified 30-May-2011




if nargin < 6
    error('gtukeycdf:TooFewInputs',...
          'Requires six input arguments.');
end


[errorcode, x, lambda1, lambda2, mu, sigma1, sigma2] = distchck(6,x,lambda1,lambda2,mu,sigma1,sigma2);

if errorcode > 0
    error('gtukeyinv:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L1=lambda1; L2=lambda2; u=mu; s1=sigma1; s2=sigma2;

% Initialize x to NaN
if isa(x,'single') || isa(L1,'single') || isa(L2,'single') || isa(u,'single') || isa(s1,'single') || isa(s2,'single')
   p = zeros(size(x),'single');
   seps=sqrt(eps('single'));
else
   p = zeros(size(x));
   seps=sqrt(eps);
end


%For the boundary cases
k0=(x==-Inf);
k1=(x==u-(s1./L1) & L1>0);
k2=(x==Inf);
k3=(x==u+(s2./L2) & L2>0);
if any(k0|k1)
    %k0 Unbounded at negative infinity
    %k1 Bounded at minimum
    p(k0|k1)=0;
end
if any(k2|k3)
    %k2 Unbounded at positive infinity
    %k3 Bounded at maximum
    p(k2|k3)=1;
end




%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;

%All other cases, that are not boundary cases
k=find(~(k0|k1|k2|k3)); %not any of the cases above
if isempty(k)
    return;
end
xk=x(k);


% Use 0.5 a starting guess
pk=0.5.*ones(size(xk));
if isa(x,'single')
    pk=single(pk);
end



h=ones(size(xk));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to p).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(pk)) && max(abs(h))>crit && count<count_limit),
    count=count+1;
    %The derivative of the inverse is explicitly known, subfunction below
    h=(gtukeyinv(pk,lambda1(k),lambda2(k),mu(k),sigma1(k),sigma2(k))-xk) ./ ...
        gtukeyinvderv(pk,lambda1(k),lambda2(k),mu(k),sigma1(k),sigma2(k));
    pnew=pk-h;
    
    %Make sure that the values stay inside the bounds (0<=p<=1)
    %Initially, Newton's Method may take big steps
    psmall = find(pnew<=0);
    plarge = find(pnew>=1);
    if any(psmall) || any(plarge)
        pnew(psmall) = abs(pk(psmall)) / 10;
        pnew(plarge) = 1 - abs(1-pk(plarge))/10;
    end
    pk=pnew;
end

%Return the converged value(s).
p(k)=pk;

if count==count_limit
    fprintf('\nWarning: TUKEYCDF did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


%Reshape to original input
p=reshape(p,size(x));


end


function dinv=gtukeyinvderv(p,L1,L2,u,s1,s2)
%The parameter u is never used in the derivative
%Parameters have passed error checking


%Seven cases
k1=(p>1)|(p<0);
k2=(p==0&L1>0)|(p==1&L2>0);
k3=(p>0&p<1)&(L1~=0&L2~=0);
k4=(L1<=0&p==0)|(L2<=0&p==1);
k5=(p>0&p<1)&(L2==0&L1~=0);
k6=(p>0&p<1)&(L1==0)&(L2==0);
k7=~(k1|k2|k3|k4|k5|k6);


if any(k1)
    dinv(k1)=NaN;
elseif any(k2)
    dinv(k2)=0;
elseif any(k3)
    pk=p(k3); L1k=L1(k3); L2k=L2(k3); s1k=s1(k3); s2k=s2(k3);
    dinv(k3)=(s1k.*(pk.^(L1k-1)))+(s2k.*((1-pk).^(L2k-1)));
elseif any(k4)
    dinv(k4)=0;
elseif any(k5)
    pk=p(k5); L1k=L1(k5); s1k=s1(k5); s2k=s2(k5);
    dinv(k5)=(s2k./(1-pk))+(s1k.*(pk.^(L1k-1)));
elseif any(k6)
    pk=p(k6); s1k=s1(k6); s2k=s2(k6);
    dinv(k6)=(s1k./pk)+(s2k./(1-pk));
else
    pk=p(k7); L2k=L2(k7); s1k=s1(k7); s2k=s2(k7);
    dinv(k7)=(s1k./pk)+(s2k.*((1-pk).^(L2k-1)));
end


dinv=dinv(:);
dinv=reshape(dinv,size(p));
end