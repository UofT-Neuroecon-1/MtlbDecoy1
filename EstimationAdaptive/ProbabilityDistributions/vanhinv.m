function x = vanhinv(y,a,pa,b,pb)
%VANHINV Inverse of the Van Houtom cumulative distribution function
%   X = VANHINV(Y,A,PA,B,PB) returns the inverse of the Van Houtom
%   cumulative distribution with parameters P(X=A)=PA, P(X=B)=PB and A<=B

%   Mike Sheppard
%   Last Modified 30-May-2011




if nargin < 5
    error('vanhinv:TooFewInputs',...
          'Requires four input arguments.'); 
end

[errorcode y a pa b pb] = distchck(5,y,a,pa,b,pb);

if errorcode > 0
    error('vanhinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end


% Initialize X to zero.
if isa(y,'single') || isa(a,'single') || isa(pa,'single') || ...
        isa(b,'single') || isa(pb,'single')
    x = zeros(size(y),'single');
else
    x = zeros(size(y));
end

% Return NaN if the arguments are outside their respective limits.
x((a>b)|(pa<0)|(pa>1)|(pb<0)|(pb>1))=NaN;

k = find(~((a>b)|(pa<0)|(pa>1)|(pb<0)|(pb>1)));
if any(k)
    term1=(-1+b(k)).*(pa(k)-y(k));
    term2=a.*(-1+pb(k)+y(k));
    term3=-1+pa(k)+pb(k);
    x(k) = ceil((term1+term2)./term3);
end

end