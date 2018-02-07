function s_hat=lindfit(x)
%LINDFIT Parameter estimate for Lindley distribution
%   S_HAT = LINDFIT(X) returns the maximum likelihood estimate of the
%   parameter S of the Lindley distribution given the data in X
%

%   Mike Sheppard
%   Last Modified: 20-May-2011


if nargin < 1
    error('lindfit:TooFewInputs',...
        'Requires at least one input arguments.');
end

x=x(:);
x(x<0)=[];

sx=sum(x);
n=length(x);

quad=sqrt(n^2+(6*n*sx)+(sx^2));
s_hat=(n-sx+quad)/(2*sx);



end