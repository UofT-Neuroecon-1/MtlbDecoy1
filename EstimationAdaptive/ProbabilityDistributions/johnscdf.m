function p = johnscdf(x,type,gamma,delta,mu,sigma)
%JOHNSCDF Johnson cumulative distribution function
%   P = JOHNSCDF(X,TYPE,GAMMA,DELTA,MU,SIGMA) returns the Johnson
%   cumulative distribution function of type TYPE, with shape 
%   parameters GAMMA, DELTA; location parameter MU, and scale 
%   parameter SIGMA, at the values in X.
%
%   Each distribution represents the distribution of the from 
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
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin < 6
    error('johnscdf:TooFewInputs',...
          'Requires six input arguments.');
end

if  ~ismember(type,{'SB','SL','SU','SN'})
    error('johnscdf:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

[errorcode,x,gamma,delta,mu,sigma] = distchck(5,x,gamma,delta,mu,sigma);

if errorcode > 0
    error('johnscdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(gamma,'single') || isa(delta,'single') || isa(mu,'single') || isa(sigma,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end

k_valid_all=(sigma>0 & delta>0);
sqr2pi=2.506628274631000502415765284811045253006986740609938316629;
sqr2=1.414213562373095;

switch type
    case 'SB' 
        k=(k_valid_all & x>mu & x<(mu+sigma));
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            xu=xk-mu; xus=-xk+mk+sk;
            term1=(gk+dk.*log(xu./xus))./sqr2;
            %Split into two cases
            k2=(xk<(uk+(sk/2)));
            pk2(k2)=(1/2).*erfc(-term1(k2));
            pk2(~k2)=(1/2).*(1+erf(term1(~k2)));
            p(k)=pk2;
        end
        p(k_valid_all & x<mu)=0;
        p(k_valid_all & x>(mu+sigma))=1;
    case 'SL'
        k=(k_valid_all & x>mu);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            term1=(gk+sk.*log((xk-mk)./sk))./sqr2;
            %Split into two cases
            k2=(xk<=(mk+sk));
            pk2(k2)=(1/2).*erfc(-term1(k2));
            pk2(~k2)=(1/2).*(1+erf(term1(~k2)));
            p(k)=pk2;
        end
        p(k_valid_all & x<mu)=0;
    case 'SU'
        k=(k_valid_all);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            term1=(gk+sk.*asinh((xk-mk)./sk))./sqr2;
            p(k)=(1/2)+(1+erf(term1));         
        end
    case 'SN'
        k=(k_valid_all);
        if any(k)
            xk=x(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            term1=(gk+(dk.*(xk-mk)./sk))./sqr2;
            p(k)=(1/2).*erfc(-term1);    
        end
    otherwise
    error('johnscdf:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

%Round off and non valid inputs
p(~k_valid_all)=NaN;



end