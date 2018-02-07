function y = gtukeypdf(x,lambda1,lambda2,mu,sigma1,sigma2)
%GTUKEYPDF Probability density function for the Generalized Tukey Lambda Distribution
%   Y = GTUKEYPDF(X,LAMBDA1, LAMBDA2, MU, SIGMA1, SIGMA2) returns the 
%   probability density function of the Generalized Tukey Lambda
%   Distribution with location parameter MU, scale parameters SIGMA1
%   and SIGMA2, and shape parameters LAMBDA1, LAMBDA2, evaluated at the
%   values in X.
%
%   Type: Continuous, unbounded / bounded
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%


%   Mike Sheppard
%   Last Modified 5-Jun-2011



if nargin < 6
    error('gtukeypdf:TooFewInputs',...
          'Requires six input arguments.');
end


[errorcode, x, lambda1, lambda2, mu, sigma1, sigma2] = distchck(6,x,lambda1,lambda2,mu,sigma1,sigma2);

if errorcode > 0
    error('gtukeypdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end

%for ease of programming
L1=lambda1; L2=lambda2; u=mu; s1=sigma1; s2=sigma2;
sz=size(L1);

%Only the inverse is given
%CDF uses newton's method to find answer
%PDF needs 1/derivative for answer

%First, use CDF to find P
p = gtukeycdf(x,lambda1,lambda2,mu,sigma1,sigma2);

%Now with p. calculate derivative of the inverse at those points
dinv=gtukeyinvderv(p,L1,L2,u,s1,s2);

%PDF is 1/dinv;
y=1./dinv;

%Reshape
y=reshape(y,sz);

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