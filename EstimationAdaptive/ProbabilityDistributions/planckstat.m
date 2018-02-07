function [m,v]=planckstat()
%PLANCKSTAT Mean and variance for the Planck's Radiation Function
%   [M,V] = PLANCKSTAT() returns the mean and variance of the Planck's
%   Radiation Function
%
%   The Planck's Radiation Function is given by
%   f(x) = (15/pi^4) * 1/((x^5)(exp(1/x)-1))
%
%   Type: Continous, Semi-Bounded
%

%   Mike Sheppard
%   Last Modified 16-Jun-2011

apery=1.202056903159594285399738161511449990764986292; %Apery's Constant

m=30*apery/(pi^4);
m2=5/(2*pi^2); %Second moment
v=m2-m^2;


end