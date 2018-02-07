function y = poissbinocdf(x,p)
%POISSBINOCDF Poisson binomial cumulative distribution function
%   y = POISSBINOCDF(X,P) returns the Poisson Binomial cumulative distribution 
%   function of the sum of N independent Bernoulli trials given by the
%   probability vector P of length N; evaluated at values in X.
%
%   If P(i)=p for all i=1,...,N then the Poisson Binomial Distribution
%   simplifies to Binomial Distribution
%   POISSBINOCDF(X,P)=BINOCDF(X,N,p)
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      x>=0 (integer)
%      0<=P(i)<=1, for all i
%
%   The size of Y is the size of the input variable X.
%

%   Mike Sheppard
%   Last Modified 15-Jun-2011


if nargin < 2
    error('poissbinocdf:TooFewInputs',...
          'Requires two input arguments.');
end


%Convert X to column vector, P to row vector, for calculation
%purposes
sz=size(x);
x=x(:); p=(p(:))';


%Return error if any P is out of bounds
if any(p<0) || any(p>1)
    error('poissbinocdf:Input',...
          'Input vector P must be between 0 and 1.');
end


% Initialize p to zero.
if isa(x,'single') || isa(p,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


%Eliminate any zero elements of vector P
p(p==0)=[]; 
N=size(p,2);

%Return values for out of bounds
y(x<0)=0;
y(x>N)=1;

%Algorithm based on 
% Fernandez, M.; Williams, S.; , "Closed-Form Expression for the 
% Poisson-Binomial Probability Density Function," Aerospace and
% Electronic Systems, IEEE Transactions on , vol.46, no.2, 
% pp.803-817, April 2010

k=(x>=0 & x==round(x) & x<=N);
if any(k),
    alpha=prod(p);
    s=-(1-p)./p;
    S=poly(s);
    temp_P=alpha.*S;
    P=temp_P(N+1:-1:1); %Lists them all from 0 to N
    Q=cumsum(P);
    y(k)=Q(x(k)+1);     %Q(1)=CDF(X=0)
end

%Reshape to original input
y=reshape(y,sz);


end