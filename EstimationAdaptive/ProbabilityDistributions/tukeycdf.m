function p = tukeycdf(x,lambda,mu,sigma)
%TUKEYCDF Cumulative function of the Tukey Lambda distribution
%   P = TUKEYCDF(X,LAMBDA,MU,SIGMA) returns the CDF of the Tukey Lambda
%   distribution with shape parameter LAMBDA, location parameter MU, and
%   scale parameter SIGMA.
%
%   TUKEYCDF(X, LAMBDA) uses the default values for MU=0, SIGMA=1
%   TUKEYCDF(X, LAMBDA, MU) uses the default value SIGMA=1;
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   TUKEYCDF uses Newton's method to converge to the solution, as there is
%   no closed form expression other than for TUKEYINV.


%   Mike Sheppard
%   Last Modified 30-May-2011




if nargin < 2
    error('tukeycdf:TooFewInputs',...
          'Requires at least two input arguments.');
end

if nargin==3
    sigma=1;
end

if nargin==2
    mu=0;
    sigma=1;
end

[errorcode, x,lambda,mu, sigma] = distchck(4,x,lambda,mu,sigma);

if errorcode > 0
    error('tukeycdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L=lambda; u=mu; s=sigma;

% Initialize x to NaN
if isa(x,'single') || isa(L,'single') || isa(u,'single') || isa(s,'single')
   p = zeros(size(x),'single');
   seps=sqrt(eps('single'));
else
   p = zeros(size(x));
   seps=sqrt(eps);
end


%For the boundary cases
k0=(x==-Inf);
k1=(x==u-(s./L) & L>0);
k2=(x==Inf);
k3=(x==u+(s./L) & L>0);
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
    h=(tukeyinv(pk,lambda(k),mu(k),sigma(k))-xk) ./ tukeyinvderv(pk,lambda(k),mu(k),sigma(k));
    pnew=pk-h;
    
    %Make sure that the values stay inside the bounds (0<=p<=1)
    %Initially, Newton's Method may take big steps
    psmall = find(pnew<=0);
    plarge = find(pnew>=1);
    if any(psmall) || any(plarge)
        pnew(psmall) = pk(psmall) / 10;
        pnew(plarge) = 1 - (1-pk(plarge))/10;
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


function dinv=tukeyinvderv(p,L,u,s)
%The parameter u is never used in the derivative
%Parameters have passed error checking
%Five cases
k1=(p>1)|(p<0);
k2=(p==0&L>0)|(p==1&L>0);
k3=(p>0&p<1&L>0)|(p>0&p<1&L<0);
k4=(L<=0&p==0)|(L<=0&p==1);
k5=~(k1|k2|k3|k4);
if any(k1)
    dinv(k1)=NaN;
elseif any(k2)
    dinv(k2)=0;
elseif any(k3)
    Lk=L(k3); pk=p(k3); sk=s(k3);
    dinv(k3)=sk.*((pk.^(Lk-1))+((1-pk).^(Lk-1)));
elseif any(k4)
    dinv(k4)=0;
else
    %k5 (L==0)
    sk=s(k5); pk=p(k5); %Note: L and u are not used
    dinv(k5)=(sk./pk)+(sk./(1-pk));
end

dinv=dinv(:);
dinv=reshape(dinv,size(p));
end