% FILE : total_resistance.m
% author : Brandon Dutcher
% signed off by: Nolan McCarthy
% date verified : 1 / 16 / 20


function f = total_resistance(thickness,env_temp)
%  calulate heat loss out of oxygen tank
%  (W/mK)
%
% Args:
%      thickness (float): insulation thickness (inches)
%      env temp (float): environmental temperature (Kelvin)
% 
% Returns:
%      Rreal (float): The thermal resistance of the tank insulation, SI units
% in SI units


load Air.dat;
OR = 8.63/2 %inches;
OR = (OR+thickness)*2.54 ./100;
OT = 200;
%env_temp = 320;
TR = 8.625*2.54/100/2; % in to m
L = 20*2.54/100/2; % in to m
k = 0.4;
Rcond = log(OR/TR)./(2*pi*L*k);
hconv = CylNuT(2,OR*2,Air,OT);
Rconv = 1./hconv./(L*OR*2*pi);
hrad = 0.6*5.67*10^-8*(OT^2+env_temp^2)*(OT+env_temp);
Rrad = 1./hrad./(L*OR*2*pi);
Rreal = Rcond + 1./(1./Rconv + 1./Rrad);
f = Rreal;
end