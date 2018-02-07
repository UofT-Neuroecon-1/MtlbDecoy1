function s = batessf(x,n,a,b)
%BATESSF Bates survival function
%   S = BATESSF(x,n,a,b) returns the survival function of the
%   Bates Distribution of the mean of n independent random
%   variables uniformly distributed from A to B.
%
%   The size of S is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for A and B are 0 and 1, respectively.
%
%   Type: Continuous, bounded, [A,B]
%   Restrictions:
%        A < B
%        N > 1    (integer)
%
%   Note: Loss of accuracy occurs when N is greater than 25.
%
%   See also BATESPDF, BATESCDF, BATESINV, BATESSTAT, BATESFIT, BATESLIKE,
%            BATESRND, BATESHAZ
%

%   Mike Sheppard
%   Last Modified: 20-Dec-2011


if (nargin~=2)&&(nargin~=4)
    error('batessf:TooFewInputs',...
        'Requires either two or four input arguments.');
end

if nargin==2
    a=0;
    b=1;
end

try
    s=1-batescdf(x,n,a,b);
catch
    error('batessf:InputSizeMismatch',...
        'Requires non-scalar arguments to match in size.');
end

end
