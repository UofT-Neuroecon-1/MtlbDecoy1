function y = voigtpdf(z,delta,sigma)
%VOIGTPDF Voigt Profile density function
%   Y = VOIGTPDF(Z,DELTA,SIGMA) returns the Voigt Profile density function
%   with scale parameters DELTA and SIGMA.
%
%   If SIGMA=0 the function reduces to the Cauchy Distribution
%
%   The size of Y is the common sizes of the inputs. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%

%   Mike Sheppard
%   Last Modified 3-Jun-2011


if nargin < 3
    error('voigtpdf:TooFewInputs',...
          'Requires three input arguments.'); 
end

[errorcode z delta sigma] = distchck(3,z,delta,sigma);

if errorcode > 0
    error('voigtpdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(z,'single') || isa(delta,'single') || isa(sigma,'single')
    y=zeros(size(z),'single');
else
    y=zeros(size(z));
end

%If sigma=0 then reduces to Cauchy distribution
k0=(isreal(z) & sigma==0); k1=(isreal(z) & sigma>0);

if any(k0)
    %Reduces to Cauchy Distribution
    y(k0) = cauchypdf(z(k0),0,delta(k0));
end

if any(k1)
i=1i; %pre-define imaginary unit
zk=z(k1); dk=delta(k1); sk=sigma(k1);
var=(zk+i.*dk)./(sk.*sqrt(2));
term1=1-i.*erfi(conj(var));
term2=exp(2*i*zk.*dk./(sk.^2));
term3=1+term1.*term2+(i*erfi(var));
fullterm=exp(-(var.^2)).*term3./(2*sqrt(2*pi).*sk);
%Should only be real, may not be so due to roundoff
y(k1)=real(fullterm);
end


% Return NaN if not real or not defined
y( ~isreal(z) | ~isreal(delta)| ~isreal(sigma) | delta<0 | sigma<0) = NaN;



end


function out=erfi(z)
%Matlab can not do the imaginary error function erfi=erf(iz)/i where z is
%complex. An alternative is to use the Faddeeva function, address given
%below, and transform it.
if ~exist('faddeeva.m','file')
    error('Requires faddeeva.m file, available on FEX #22207. See code for URL');
    %URL: http://www.mathworks.com/matlabcentral/fileexchange/22207-faddeeva-function-fft-based
end
i=1i; %pre-define imaginary unit
out=i*(1-(faddeeva(z)./exp(-(z.^2))));
end

