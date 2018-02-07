function x = johnsinv(p,type,gamma,delta,mu,sigma)
%JOHNSINV Inverse of the Johnson cumulative distribution function
%   X = JOHNSINV(P,TYPE,GAMMA,DELTA,MU,SIGMA) returns the inverse of the 
%   Johnson cumulative distribution function of type TYPE, with shape 
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
%      0<=P<=1
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%

%   Mike Sheppard
%   Last Modified 26-Jun-2011



if nargin < 6
    error('johnsinv:TooFewInputs',...
          'Requires six input arguments.');
end

if  ~ismember(type,{'SB','SL','SU','SN'})
    error('johnsinv:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

[errorcode,p,gamma,delta,mu,sigma] = distchck(5,p,gamma,delta,mu,sigma);

if errorcode > 0
    error('johnsinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize x to zero.
if isa(p,'single') || isa(gamma,'single') || isa(delta,'single') || isa(mu,'single') || isa(sigma,'single')
   x = zeros(size(p),'single');
else
   x = zeros(size(p));
end

k_valid_all=(sigma>0 & delta>0 & p>=0 & p<=1);
sqr2pi=2.506628274631000502415765284811045253006986740609938316629;
sqr2=1.414213562373095;

switch type
    case 'SB' 
        k=(k_valid_all);
        if any(k)
            pk=p(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            x(k)=mk+(sk./(1+exp(-(-gk+sqr2.*erfinv(-1+2*pk))./dk)));
        end
    case 'SL'
        k=(k_valid_all);
        if any(k)
            pk=p(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            x(k)=mk+sk.*exp((-gk+sqr2.*erfinv(-1+2*pk))./dk);
        end
    case 'SU'
        k=(k_valid_all);
        if any(k)
            pk=p(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            x(k)=mk+sk.*sinh((-gk+sqr2.*erfinv(-1+2*pk))./dk);      
        end
    case 'SN'
        k=(k_valid_all);
        if any(k)
            pk=p(k); gk=gamma(k); dk=delta(k); mk=mu(k); sk=sigma(k);
            x(k)=mk+sk.*((-gk+sqr2.*erfinv(-1+2*pk))./dk);     
        end
    otherwise
    error('johnscdf:Type',...
          'TYPE must be one of the following: ''SB'' , ''SL'' , ''SU'' , or ''SN''.');
end

%Round off and non valid inputs
x(~k_valid_all)=NaN;


end