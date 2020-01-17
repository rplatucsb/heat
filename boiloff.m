% FILE : boiloff.m
% author : Sebastian Vargas, Nolan McCarthy
% signed off by: Nolan McCarthy
% date verified : 1 / 16 / 20

function f = boiloff(chosen_insulation_thickness,t1 )
%Calculate the rate of boil off of the tanks on the
%test stand.
%Args :
%    chosen_insulation_thickness ( float ) : thickness of insultion , inches
%    t1 ( float ) : temperature of the air, Ferenheit
%Returns :
%    boil_off_litres_hour ( float ) : the rate of boil off in the tanks in Litres/ boil_off_litres_hour

ambient_temperature = (t1-32)*5/9 +273.15; %conversion to kelvin
cryogenic_temperature = 152; %conversion to kelvin
Ts = cryogenic_temperature; %temperature of wall surface;
temperature_difference = (ambient_temperature - cryogenic_temperature);
latent_heat_vaporization_lox = 3.4099*1000; %J/mol
grams_liter_lox = 1140; %g/L
molar_mass_lox = 32; %g/mol
tank_length = 14 * .0254; %inches to meters
tank_internal_radius = (8.63/2 - .50)*.0254; %inches to meters
surface_area_inner = tank_length*2*pi*tank_internal_radius;
styrofoam_internal_radius = 8.63/2 * .0254; %inches to meters
thickness = chosen_insulation_thickness * .0254; %inches to meters
surface_area_outer = tank_length*2*pi*(styrofoam_internal_radius + thickness);
thermal_conductivity_styrofoam = .42; %W/mk
thermal_conductivity_tank = 16.2; %W/mk
thermal_emissivity_mylar=.044; %unitless
convection_coefficient_lox = 1000; %will prove to be very insignificant
convection_coefficient_air = 35; %W/m2C
sigma = 5.6704*10^(-8); %Stefan_Boltzmann constant W/m^2/K^4

Q = temperature_difference ./total_resistance(chosen_insulation_thickness,ambient_temperature);
boil_off_rate_mol = (Q/latent_heat_vaporization_lox); %mol/s
boil_off_rate = boil_off_rate_mol*molar_mass_lox/(grams_liter_lox);
boil_off_litres_hour = boil_off_rate*60*60;
f = boil_off_litres_hour;
end