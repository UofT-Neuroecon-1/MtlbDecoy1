function x = bnkgbinv(p,a,b)
%BNKGBINV Inverse of the Benktander-Gibrat cumulative distribution function
%   X = BNKGBINV(P,A,B) returns the inverse cumulative distribution
%   function of the Benktander-Gibrat distribution with parameters A and B
%   at the values in P.
%
%   The size of X is the common sizes of the inputs. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   Distribution: Continuous, semi-bounded, [1,Inf)
%   Restrictions:
%      A(A+1) >= 2B
%      A , B > 0
%
%   Note: The Benktander-Gibrat Distribution is also known as the
%   Benktander Distribution of Type I. BNKGBINV uses Newton's method
%   to converge to the solution
%
%   See also BNKGBPDF, BNKGBCDF, BNKGBSTAT, BNKGBFIT, 
%            BNKGBLIKE, BNKGBRND, BNKGBSF, BNKGBHAZ
%

%   Mike Sheppard
%   Last Modified: 16-May-2012



if nargin ~= 3
    error('bnkgbinv:TooFewInputs',...
        'Requires three input arguments.');
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('bnkgbinv:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(p,'single') || isa(a,'single') || isa(b,'single')
    x=zeros(size(p),'single');
    seps=sqrt(eps('single'));
else
    x=zeros(size(p));
    seps=sqrt(eps);
end

% Return NaN if any arguments are outside of their respective limits.
okparam = (0<a & a<Inf) & (0<b & b<Inf) & (a.*(a+1) >= 2*b);
okvar = (0 < p & p < 1);
ok=(okparam & okvar);
x(~ok)=NaN;
x(okparam & p==0)=1;
x(okparam & p==1)=Inf;


%Newton's Method
%Permit no more than count_limit iterations
count_limit=100;
count=0;


p=p(ok); a=a(ok); b=b(ok);


% Use the mean as a starting guess
xok = 1+(1./a);
if isa(p,'single')
    xok=single(xok);
end


%Move starting values away from the boundaries
xok(xok==1)=1+seps;
h=ones(size(p));
crit=seps;


%Break out of the iteration loop for the following:
% 1) The last update is very small (compared to x).
% 2) The last update is very small (compared to 100*eps)
% 3) There are more than 100 iterations. This should NEVER happen.

while(any(abs(h)>crit*abs(xok)) && count<count_limit),
    count=count+1;
    h=(bnkgbcdf(xok,a,b)-p)./bnkgbpdf(xok,a,b);
    xnew=xok-h;
    
    %Make sure that the values stay inside the bounds.
    %Initially, Newton's Method may take big steps
    ksmall = find(xnew<=1);
    if any(ksmall)
        xnew(ksmall)=1+(1-xnew(ksmall))/10;
    end
    
    xok=xnew;
end

%Return the converged value(s).
x(ok)=xok;

if count==count_limit
    fprintf('\nWarning: BNKGBINV did not converge.\n');
    str='The last step was: ';
    outstr=sprintf([str,'%13.8f\n'],max(h(:)));
    fprintf(outstr);
end



end