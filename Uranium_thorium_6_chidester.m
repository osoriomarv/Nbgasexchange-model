function [helium_temporal, Dm_U, Dm_Th, U_left, Th_left, He3_ratio, helium_temp_0_Gya, helium_temp_4_Gya] = Uranium_thorium_6_chidester(He_Conc_metal, Tadi, Padi, frac, Upart_fit, Thpart_fit)
%output a function of planetary radius
Temp_adiabat = Tadi;
Pressure_adiabat = Padi;
He_conc_metal = He_Conc_metal;
frac = frac;
Dm_U = 10^Upart_fit;
Dm_Th = 10^Thpart_fit;

% Temp_adiabat = 2888; %K
% Pressure_adiabat = 10; %GPa
% He_conc_metal = 1e-4; 
% X_SiO2=.4;
% frac = 1;
% IW = -2;
% mol_O_metal= .1;
% Dfe = .3/.2;


% 3He/4He ratio during solar nebula - from RiMG
He_ratio=1.66e-4;
%current atmospheric ratio
atm_3_4_he=1.34e-6;

He_4=He_conc_metal/4; 
%initial he_3 solar nebula concentration
He_3=He_ratio.*He_4;

% Atomic mass
U_amu = 238.03806; 
Th_amu = 232.03806; 
% Uranium and thorium concentration
U_conc_py = 20.3; % for pyrolite PPB
Th_conc_py  = 79.5; % PPB 

% Conversion to mol g
U_py_mol_g = U_conc_py./1e6./U_amu;
Th_py_mol_g = Th_conc_py./1e6/Th_amu; 
t=4.64e9;
% Constants
half_life_U_238 = 4.468e9;      % Half-life of Uranium-238 in years
half_life_U_235 = 7e8;  %Half-life of Uranium-235 in years
half_life_Th_232 = 14.05e9;     % Half-life of Thorium in years
lambda_U_235 = log(2)/half_life_U_235;      % Decay constant for Uranium-235
lambda_U_238 = log(2)/half_life_U_238;      % Decay constant for Uranium-238
lambda_Th_232 = log(2)/half_life_Th_232;    % Decay constant for Thorium-232
% Given values
U_current = U_py_mol_g;   % Current amount of Uranium in mol/g
Th_current_232  = Th_py_mol_g;  % Current amount of Thorium in mol/g
U_current_235 = U_current * 0.00720; %mol_g
U_current_238 = U_current * 0.99275; %mol_g

% Exponential decay equation
Th_232 = Th_current_232 * exp(lambda_Th_232*t);
U_235 = U_current_235 * exp(lambda_U_235*t);
U_238 = U_current_238 * exp(lambda_U_238*t);

U_fraction_235 = U_235/U_238 + U_238;
U_fraction_238 = (1- U_fraction_235);

%Bulk Conc
%Potentially change the .66 to the fraction that is BSE
%Pass Frac in - that will determine the bulk
Conc_bulk_Th = Th_232 * (1 - (1/3 * frac)); 
Conc_bulk_U = (U_235+U_238) * (1 - (1/3 * frac));

for i=1:numel(Temp_adiabat)

    %1/6
    Conc_metal_U(i) = Conc_bulk_U * Dm_U(i) ./ ((1-(1/3*frac))* (1-Dm_U(i)) + Dm_U(i));
    Conc_metal_Th(i) = Conc_bulk_Th * Dm_Th(i) ./ ((1-(1/3 * frac)) * (1-Dm_Th(i)) + Dm_Th(i));

    Th_mol_metal(i) = Conc_metal_Th(i);
    U235_mol_metal(i) = Conc_metal_U(i) * U_fraction_235;
    U238_mol_metal(i) = Conc_metal_U(i) * U_fraction_238;

    time_vector=(0:1e6:4.5e9); 
     for tt=1:numel(time_vector)
        
        %half_life of U_238 is 4.47 GA 
        U238_left=U238_mol_metal(i).*exp(-1.54e-10.*time_vector); %mol/g
        % half_life of U_235 is 700 ma
        U235_left=U235_mol_metal(i).*exp(-9.72e-10.*time_vector); %mol/g
        %half life of Th is 14 Ga
        Th232_left=Th_mol_metal(i).*exp(-4.940e-11.*time_vector); %mol/g
        
        %Uranium 238 will produce 8 moles
        %uranium 235 will produce 7 moles
        %Th 232 will produce 6
        
        alpha_Th232_mol=(Th_mol_metal(i)-Th232_left).*(6); %mol
        alpha_U235_mol=(U235_mol_metal(i)-U235_left).*(7);%mol 
        alpha_U238_mol=(U238_mol_metal(i)-U238_left).*(8);%mol  
        
        %at the one mole rate - need to include the rest of the moles 
        %mols 
        He_4_temporal=He_4+alpha_Th232_mol+alpha_U235_mol+alpha_U238_mol;
        %mol in core from calculation
        %mol/mol 
        he3_ratio_temporal=He_3./He_4_temporal;
        
        helium_temporal = he3_ratio_temporal/atm_3_4_he;

        helium_temp_0_Gya = helium_temporal(1);
        helium_temp_4_Gya = helium_temporal(4501);
        
        He3_ratio = he3_ratio_temporal;

        U_left = U238_left + U235_left;
        Th_left = Th232_left;
        % % 
        % U_left = min(U238_left + U235_left, U235_mol_metal + U238_mol_metal);
        % Th_left = min(Th232_left, Th_mol_metal);
     end
end
% 
% 
figure(2); hold on
plot(time_vector, helium_temporal, 'b' , LineWidth=2)
xlabel('Time (Yrs)')
ylabel('3/4 Helium R/R_{A}')
title('Helium Isotope Ratio Over Time')
%legend('Helium 3/4 Ratio')
legend(sprintf('Helium 3/4 Ratio\nMainMO(.2, 5000, 1, 7,-2, 1,.25, .4)'))
ca = gca;
ca.FontSize = 16;
