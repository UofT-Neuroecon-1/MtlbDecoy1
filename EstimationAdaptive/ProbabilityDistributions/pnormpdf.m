function y = pnormpdf(x,lambda,mu,sigma)
%PNORMPDF Power Normal probability density function (pdf).
%   Y = NORMPDF(X,LAMBDA, MU,SIGMA) returns the probability density
%   function of the Power Normal Distribution with Box-Cox variable LAMBDA,
%   mean MU, and standard deviation SIGMA, evaluated at the values in X.
%
%  X>0 

% http://www.udc.edu/docs/dc_water_resources/technical_reports/report_n_190.pdf


if nargin<4
    error('pnormpdf:TooFewInputs','Requires four input argument.');
end



%Function P_Lambda(X)
pl=((x.^lambda)-1)/lambda;
%pl=ln(x) if lambda==0
    
[errorcode x lambda mu sigma] = distchck(4,x,lambda,mu,sigma);

if errorcode > 0
    error('pnormpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end
   
% Initialize Y to zero.
if isa(x,'single') || isa(lambda,'single') || isa(mu,'single') || isa(sigma,'single') 
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


%First compute Box-Cox transformation
plx=zeros(size(y)); %same size
k0=(lambda==0);
plx(k0)=log(x(k0)); %If lambda==0
plx(~k0)=((x(~k0).^lambda(~k0))-1)/lambda(~k0); %Otherwise

%coefficient of variation
k=sigma./mu;
T=(1./(lambda.*sigma))+(1./k);
Kf=ones(size(y));
Kf(lambda>0)=normcdf(T(lambda>0));
Kf(lambda==0)=0;
Kf(lambda<0)=normcdf(-T(lambda<0));

%Exponential term
expterm=exp((-1/2)*((plx-mu)/sigma).^2);
%Power term
powx=x.^(lambda-1);
%Coef term
coef=1./(sqrt(2*pi).*sigma);

y=(1./Kf).*coef.*powx.*expterm;

end