function y = polywblpdf(x,a,b)
%POLYWBLPDF Poly-Weibull probability density function
%   Y = POLYWBLPDF(X,A,B) returns the probability density function of the
%   Poly-Weibull Distribution of the minimum of N statistically independent
%   random variables having non-identical Weibull distributions, with scale
%   parameters A and shape parameters B, evaluated at values in X.
%
%   A and B are vectors of length N, representing the N non-identical
%   Weibull distributions.
%
%   Type: Continuous, Semi-bounded
%   Restrictions:
%      X > 0
%      A, B > 0  (every element)
%
%   The size of Y is the size of the input variable X.
%

%   Reference:
%   http://www.stat.purdue.edu/research/technical_reports/pdfs/1992/tr92-05c.pdf

%   Mike Sheppard
%   Last Modified 4-Jul-2011


if nargin < 3
    error('polywblpdf:TooFewInputs',...
        'Requires three input arguments.');
end


%Convert X to column vector, A and B to row vectors, for calculation
%purposes
sz=size(x);
x=x(:); a=(a(:))'; b=(b(:))';

%If A or B is given as a scalar, turn into vector with common size
%Also returns error if dimension mis-match
try
    a=a+zeros(size(b));
    b=b+zeros(size(a));
catch
    error('polywblcdf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end


%Return error if any A or B are negative
if any(a<0) || any(b<0)
    error('polywblcdf:Input',...
        'Inputs A and B are required to be positive.');
end


% Initialize y to zero.
if isa(x,'single') || isa(a,'single') || isa(b,'single')
    y=zeros(size(x),'single');
else
    y=zeros(size(x));
end


%If any x are negative, those y's return as NaN
y(x<0)=NaN;

k=(x>0);
if any(k),
    sk=sum(k);
    %Repmat each A and B to size of k
    aM=repmat(a,sk,1);
    bM=repmat(b,sk,1);
    %Repmat x to size of N
    xM=repmat(x(k),1,size(a,2));
    %Survival for each (a,b) pair
    z=(xM./aM).^bM;
    Sj=exp(-z);
    %Overall reliability
    S=prod(Sj,2);
    %Sum-matrix
    SM=(bM.*(xM.^(bM-1)))./(aM.^bM);
    SM=sum(SM,2);
    %Probability density
    y(k)=SM.*S;
end

%Reshape to original input
y=reshape(y,sz);



end