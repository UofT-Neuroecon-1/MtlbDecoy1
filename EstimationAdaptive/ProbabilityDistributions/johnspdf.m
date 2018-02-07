function y = johnspdf(x,type,gamma,delta,mu,sigma)
%JOHNSPDF Johnson probability density function
%   Y = JOHNSPDF(X,TYPE,GAMMA,DELTA,MU,SIGMA) returns the Johnson
%   Distribution of type TYPE, with shape parameters GAMMA, DELTA; location
%   parameter MU, and scale parameter SIGMA, at the values in X.
%
%   Each distribution represents the distribution of the form 
%      Y ~ mu + sigma * G((X-gamma)/delta)
%   Where X is from the Normal Distribution
%
%   TYPE can be one of {'SB','SL','SU','SN'}
%      SB: Bounded, corresponding to G(x)=1/(1+exp(-x))
%      SL: Semi-bounded, corresponding to G(x)=exp(x)
%      SU: Unbounded, corresponding to G(x)=sinh(x)
%      SN: Normal, corresponding to G(x)=x
%
%   Restrictions:
%      SIGMA, DELTA > 0
%      Case 'SB' : MU < X < MU+SIGMA
%      Case 'SL' : X > MU
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin < 6
    error('johnspdf:TooFewInputs',...
          'Requires six input arguments.');
end

if  ~ismember(type,{'SB','SL','SU','SN'})
    error('johnspdf:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

[errorcode,x,gamma,delta,mu,sigma] = distchck(5,x,gamma,delta,mu,sigma);

if errorcode > 0
    error('johnspdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(gamma,'single') || isa(delta,'single') || isa(mu,'single') || isa(sigma,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

k_valid_all=(sigma>0 & delta>0);
sqr2pi=2.506628274631000502415765284811045253006986740609938316629;

switch type
    case 'SB' 
        k=(k_valid_all & x>mu & x<(mu+sigma));
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            xu=xk-mu; xus=-xk+mk+sk;
            term1=(gk+dk.*log(xu./xus)).^2;
            num=exp((-1/2).*term1).*dk.*sk;
            den=sqr2pi.*xu.*xus;
            y(k)=num./den;
        end
    case 'SL'
        k=(k_valid_all & x>mu);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            num=dk.*exp((-1/2).*((gk+sk.*log((xk-mk)./sk)).^2));
            den=sqr2pi.*(xk-uk);
            y(k)=num./den;
        end
    case 'SU'
        k=(k_valid_all);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            num=dk.*exp((-1/2).*((gk+sk.*asinh((xk-mk)./sk)).^2));
            den=sqr2pi.*sqrt((xk-uk).^2+sk.^2);
            y(k)=num./den;            
        end
    case 'SN'
        k=(k_valid_all);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            %Transformation of the Normal Distribution
            mu_norm=gk-(mk.*dk./sk);
            sigma_norm=dk./sk;
            y(k)=normpdf(x,mu_norm,sigma_norm);        
        end
    otherwise
    error('johnspdf:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

%Round off and non valid inputs
y(~k_valid_all)=NaN;
y(y<0)=0;


end