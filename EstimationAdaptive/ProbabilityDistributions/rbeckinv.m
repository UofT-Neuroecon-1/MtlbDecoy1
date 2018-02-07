function x = rbeckinv(p,s1,s2)
%RBECKINV Inverse of the Reduced Beckmann cumulative distribution function
%   X = RBECKINV(P,S1,S2) returns the inverse cumulative distribution 
%   function of the Reduced Beckmann distribution with standard deviations 
%   S1 and S2, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   The full Beckmann Distribution
%   If {x,y} follows the Bivariate Normal Distribution with means {u1,u2},
%   standard deviations {s1,s2} and correlation rho, then sqrt[x^2+y^2]
%   follows the Beckmann Distribution [u1,u2,s1,s2,rho]
%
%   The reduced Beckmann distribution assumes u1=u2=rho=0 so the bivariate
%   normal distributions are centered at (0,0) with covariance matrix
%   [s1^2 0; 0 s2^2]
%
%   Distribution: Continuous, semi-bounded, [0,Inf) 
%   Restrictions:
%      S1, S2 > 0
%
%   NOTE: RBECKINV uses Newton's method to converge to the solution.
%
%   See also RBECKPDF, RBECKCDF, RBECKSTAT, RBECKFIT, 
%            RBECKLIKE, RBECKRND, RBECKSF, RBECKHAZ
%

%   Mike Sheppard
%   Last Modified 7-May-2011

if nargin ~= 3
    error('rbeckinv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode p s1 s2] = distchck(3,p,s1,s2);

if errorcode > 0
    error('rbeckinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

%Initialize x to zero
if isa(p,'single') || isa(s1,'single') || isa(s2,'single')
    x=zeros(size(p),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(p));
    seps=sqrt(eps);
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<s1 & s1<Inf) & (0<s2 & s2<Inf);
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=0;
x(okparam & p==1)=Inf;



%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;

if all(~ok),  return; end
p=p(ok); s1=s1(ok); s2=s2(ok);


% Use the mean as a starting guess
xok = rbeckstat(s1,s2);
if isa(p,'single')
    xok=single(xok);
end

%Move starting values away from the boundaries
xok(xok==0)=seps;
h=ones(size(p));
crit=seps;

%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xok)) && count<count_limit),
    count=count+1;
    h=(rbeckcdf(xok,s1,s2)-p)./rbeckpdf(xok,s1,s2);
    xnew=xok-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = find(xnew<=0);
    if any(ksmall)
        xnew(ksmall)=xok(ksmall)/10;
    end
    xok=xnew;
end

%Return the converged value(s).
x(ok)=xok;

if count==count_limit
    fprintf('\nWarning: RBECKINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end


end