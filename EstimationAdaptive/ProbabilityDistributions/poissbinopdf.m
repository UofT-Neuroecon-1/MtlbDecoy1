function y = poissbinopdf(x,p)
%POISSBINOPDF Poisson Binomial probability density function
%   Y = POISSBINOPDF(X,P) returns the probability density function of the
%   Poisson Binomial Distribution of the sum of N independent Bernoulli
%   trials with probability vector P of length N; evaluated at the
%   values in X.
%
%   Type: Discrete, Semi-bounded
%   Restrictions:
%      X >= 0      (integer)
%      0 <= P <= 1 (each element)
%
%   Note: For constant vector P the Poisson Binomial Distribution reduces
%   to the Binomial Distribution
%   POISSBINOPDF(X,P)=BINOPDF(X,N,p)
%
%   The size of Y is the size of the input variable X.
%

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 2
    error('poissbinopdf:TooFewInputs',...
          'Requires two input arguments.');
end


%Convert X to column vector, P to row vector, for calculation
%purposes
sz=size(x); x=x(:); p=(p(:))';

%Return error if any P is out of bounds
if any(p<0) || any(p>1)
    error('poissbinopdf:Input',...
          'All elements in input vector P must be between 0 and 1.');
end

% Initialize y to zero.
if isa(x,'single') || isa(p,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end

%Eliminate any zero elements of vector P
p(p==0)=[]; 
N=size(p,2);

%Return 0 for out of bounds
y(x<0 | x>N)=0;

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
    y(k)=P(x(k)+1);     %P(1)=P(X=0) and extract only inputed X values
end

%Reshape to original input
y=reshape(y,sz);


end