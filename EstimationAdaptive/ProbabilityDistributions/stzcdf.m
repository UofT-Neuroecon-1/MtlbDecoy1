function p = stzcdf(z,n)
%STZPDF Student's Z cumulative distribution function
%   P = STZPDF(z,n) returns the cumulative distribution function of the
%   Student's Z distribution
%
%   Requires 'hypergeometric2F1ODE' from the Matlab File Exchange

%   Mike Sheppard
%   Last Modified 1-Jul-2011


if nargin < 2
    error('stzcdf:TooFewInputs',...
        'Requires at least two input arguments.');
end

try
    %Expand size if necessary
    z=z+zeros(size(n));
    n=n+zeros(size(z));
catch
    error('stzcdf:InputSizeMismatch',...
        'Non-scalar arguments must match in size.');
end


% Initialize Y to zero.
if isa(z,'single') || isa(n,'single')
    p=zeros(size(z),'single');
else
    p=zeros(size(z));
end
h2f1=p;

sqrtpi=1.772453850905516;

%Function does range if given two input points; force only given points
%Also the (a,b,c) inputs must be scalars
[uniq_n,ignore,indx_n]=unique(n(:));
for i=1:length(uniq_n)
    ni=uniq_n(i); z4ni=z(indx_n==i); z4F=-(z4ni.^2);
    t_h2f1=evalh2f1((ni-1)/2,ni/2,(ni+1)/2,z4F);  %F(a,b,c;z)
    h2f1(indx_n)=t_h2f1;
end

h2f1=reshape(h2f1,size(z));  %truncate from 3rd position, 1st for data

num=(abs(z).^(1-n)).*gamma(n./2).*h2f1;
den=2.*sqrtpi.*gamma((n+1)/2);
p=num./den;
p(z>=0)=1-p(z>=0);

end


function out=evalh2f1(a,b,c,z)
%Function can not handle repeated values, or zeros, and the input has to be
%sorted and start with a zero.
orig_z=z; indx_z0=(z==0); indx_nz0=(z~=0);
%Delete any zeros
z(z==0)=[];
%Use only unique values
[z_uniq,ignore,z_indx]=unique(z(:));
%Values are sorted in ascending order, now sort in descending
z_sort=z_uniq(end:-1:1); z_indx=z_indx(end:-1:1);
%Create input index
z_min=max(z_sort(:));  %Closest to zero
z_input=[0; z_min/2; z_sort(:)];
%Solve for 2F1
[temp_h2f1 y]=hypergeometric2F1ODE(a,b,c,z_input);
temp_h2f1=temp_h2f1(3:end); %Ignore first two outputs
%Now place in correct spots
out(indx_nz0)=temp_h2f1(z_indx);
out(indx_z0)=0;
out=out(:);
end
