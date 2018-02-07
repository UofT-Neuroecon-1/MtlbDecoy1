function bhat = benfordfit(x)
%BENFORDFIT Parameter estimate for the Benford distribution
%   BENFORDFIT(X) Returns the maximum likelihood esimate of the parameter
%   of the Benford distribution given the data in the vector, X.
%
%   As B is an integer, the closest integer that matches the mean will
%   be selected. In this regard it is a Method of Moments estimator.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also BENFORDPDF, BENFORDCDF, BENFORDINV, BENFORDSTAT, 
%            BENFORDLIKE, BENFORDRND

%   Mike Sheppard
%   Last modified 7-May-2011


if nargin < 1
    error('benfordfit:TooFewInputs',...
          'Requires at least one input argument.'); 
end



%USE NEWTON'S METHOD TO CONVERGE TO NON-INTEGER SOLUTION
%---------
%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;
valk=mean(x(:)); %Value to be reached

% Use base 10 as starting point
bk=10;
h=1;
crit=sqrt(eps);

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(bk)) && count<count_limit),
    count=count+1;
    num= bk-(gammaln(bk+1)./log(bk));
    den= 1+ (gammaln(bk+1) / (bk*log(bk)^2)) - (psi(bk+1)/log(bk));
    h=(num-valk)./den;
    bnew=bk-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = find(bnew<2);
    if any(ksmall)
        bnew(ksmall)=2+eps;
    end
    bk=bnew;
end

%Return the converged value(s).

if count==count_limit
    fprintf('\nWarning: BENFORDFIT did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end

%--------


%With non-integer, find best integer
%Non integer, check above and below for closest match to estimated variance
bk = floor(bk)-2:floor(bk)+2; bk(bk<1)=[];
m2 = bk-(gammaln(bk+1)./log(bk));
[val,indx] = min((valk-m2).^2);
bhat=bk(indx);


end

