%Main Code - Combines tandpcalc, mo_core_equi_T, IPNB and Uranium_thorium.

% Main_MO(5.94612878531627,4.36124254268669,0,0.450308544743621,0.212086431458656,0,5.19294475781180, 0.298159588641994,0.0350792060089926,0,4500,5000,10,1e6, .01, .75)

% Explore the logspace of Taccr and D_he
% 1e6 - 1e8 - do logspace
% .001 -.5 log space
% 
% 
function [SiO2molfrac, AlO1molfrac, FeOmolfrac, MgOmolfrac, CaOmolfrac, T_v_K, Tadi, Pressure_bottom_pa, Padi, U_left_logger, Th_left_logger, He3_ratio_logger, Dm_U_logger, Dm_Th_logger, DIW_LOG, core_grams, Bse_grams, He_atm_grams, helium_temp_0_Gya_logger, helium_temp_4_Gya_logger,He_Conc_metal_logger, He_Core_grams_logger,He_Bse_grams_logger]...
    = Main_MO(psi_si, radius,radius2,step,Taccretion, D_he, frac, adibat, moequi)

% function [SiO2molfrac, AlO1molfrac, FeOmolfrac, MgOmolfrac, CaOmolfrac, T_v_K, Tadi, Pressure_bottom_pa, Padi, U_left_logger, Th_left_logger, He3_ratio_logger, core_grams, Bse_grams, He_atm_grams, helium_temp_0_Gya_logger, helium_temp_4_Gya_logger,He_Conc_metal_logger, He_Core_grams_logger,He_Bse_grams_logger]...
%     = Main_MO(psi_si, radius,radius2,step,Taccretion, D_he, frac, adibat, moequi)


Taccretion = 10^(Taccretion);
D_he = 10^(D_he);
Radius = linspace(radius, radius2, step);

[n_SiO2, n_MgO, n_FeO, n_Al2O3, n_CaO, n_NiO, n_Fe, n_Ni, n_O_metal, n_Si_metal] = enstatitecomp_2(psi_si);

% Call tandpcalc and store the output vectors in variables
[T_v_K, Pressure_adiabat, core_grams, Bse_grams, He_atm_grams, px_he_pa, p_v_mpa, Pressure_bottom_pa] = tandpcalc_3(radius,radius2,step, Taccretion, frac);

%Preallocation of Values
dIW=zeros(size(step));
Padi = zeros(size(step));
Tb = zeros(size(step));
Tadi = zeros(size(step));
FeOmolfrac=zeros(size(step));
NiOmolfrac=zeros(size(step));
SiO2molfrac=zeros(size(step));
Femolfrac=zeros(size(step));
Nimolfrac=zeros(size(step));
Omolfracpick=zeros(size(step));
Simolfrac=zeros(size(step));
MgOmolfrac=zeros(size(step));
AlO1molfrac=zeros(size(step));
CaOmolfrac=zeros(size(step));
core_grams_logger=zeros(size(step));
Bse_grams_logger=zeros(size(step));
He_atm_grams_logger=zeros(size(step));
He_Core_grams_logger=zeros(size(step));
He_Bse_grams_logger=zeros(size(step));
He_atm_percent_logger=zeros(size(step));
He_BSE_percent_logger=zeros(size(step));
He_Core_percent_logger=zeros(size(step));
HeEarthsum_logger=zeros(size(step));
He_Conc_metal_logger=zeros(size(step));
U_left_logger=zeros(size(step));
Th_left_logger=zeros(size(step));
He3_ratio_logger=zeros(size(step));
PPM_TH_core=zeros(size(step));
PPM_U_core=zeros(size(step));
helium_temp_0_Gya_logger=zeros(size(step));
helium_temp_4_Gya_logger=zeros(size(step));


% For loop to iterate through values of tandpcalc using the Radius values
% which set the size and step.
count =0;
for i=1: numel(Radius)
    r_v_logger(i) = Radius(i);
    count = count;
    
    pbpi = Pressure_bottom_pa(1);
    % Calculation of Adiabatic Temperature and pressure
    Padi(i) = Pressure_adiabat(i);
    Tb(i) = T_v_K(i);

    Tbat = (Radius.*moequi) * adibat;
    Tadi(i) = Tb(i) + Tbat(i);
    %Export pressure
    px_he_pa(i) =px_he_pa(i);
    p_v_mpa(i) = p_v_mpa(i);
    Pressure_bottom_pa(i) = Pressure_bottom_pa(i);
    
    % Call upon mo_core_equi_T 
    [dIW(i),A(:,i)] = mo_core_equi_T_4(n_SiO2, n_MgO, n_FeO, n_Al2O3, n_CaO, n_NiO, n_Fe, n_Ni, n_O_metal, n_Si_metal, Tadi(i), Padi(i),frac);

    %Separate Values out
    dIW_logger(i) = dIW(i);
    FeOmolfrac(i)=A(1,i);
    NiOmolfrac(i)=A(2,i);
    SiO2molfrac(i)=A(3,i);
    Femolfrac(i)=A(4,i);
    Nimolfrac(i)=A(5,i);
    Omolfracpick(i)=A(6,i);
    Simolfrac(i)=A(7,i);
    MgOmolfrac(i)=A(8,i);
    AlO1molfrac(i)=A(9,i);
    CaOmolfrac(i)=A(10,i);
    

    %Set 0 for values
    TiO2 = 1e-6;
    Fe2O3 = 1e-6;
    Na2O = 1e-6;
    K2O = 1e-6;

    %Combine to a new Matrix
    B = [SiO2molfrac; AlO1molfrac; FeOmolfrac; MgOmolfrac; CaOmolfrac];

    % core components
    core_grams_logger(i) = core_grams(i);
    Bse_grams_logger(i) = Bse_grams(i);
    He_atm_grams_logger(i) = He_atm_grams(i);

    % Pass Values to IPNB
    [He_Core_grams, He_Bse_grams, He_atm_percent, He_BSE_percent, HE_Core_percent, HeEarthsum_grams, He_Conc_metal,  negativelnxi, ccSTPgHe, IP, fugacitycoeff, He_sol_g_g] = IPNB_MRK_5(SiO2molfrac(i), TiO2, AlO1molfrac(i), FeOmolfrac(i), Fe2O3, MgOmolfrac(i), CaOmolfrac(i), Na2O, K2O, Tb(i), pbpi,px_he_pa(i),Pressure_bottom_pa(i), core_grams(i), Bse_grams(i), D_he, He_atm_grams(i));
     
    
    %negativelnxi, ccSTPgHe_bar, IP, fugacitycoeff
     negativelnxi_logger(i) = negativelnxi;
     He_sol_g_g_logger(i) = He_sol_g_g;
     ccSTPgHe_logger(i) = ccSTPgHe;
     IP_logger(i) = IP;
     fugacitycoeff_logger(i) = fugacitycoeff;


    % Output separate values out
    He_Core_grams_logger(i) = He_Core_grams;
    He_Bse_grams_logger(i) = He_Bse_grams;
    He_atm_percent_logger(i) =He_atm_percent; 
    He_BSE_percent_logger(i) = He_BSE_percent;
    He_Core_percent_logger(i) = HE_Core_percent;
    HeEarthsum_logger(i) = HeEarthsum_grams;
    He_Conc_metal_logger(i) = He_Conc_metal;
    % 
    % FeOmolfrac(i)=A(1,i);
    % NiOmolfrac(i)=A(2,i);
    % SiO2molfrac(i)=A(3,i);
    % Femolfrac(i)=A(4,i);
    % Nimolfrac(i)=A(5,i);
    % Omolfracpick(i)=A(6,i);
    % Simolfrac(i)=A(7,i);
    % MgOmolfrac(i)=A(8,i);
    % AlO1molfrac(i)=A(9,i);
    % CaOmolfrac(i)=A(10,i);
    % 
    % [Ca, Al, Mg, Fe, Ni, O, Si] = silicatecomponents(CaOmolfrac(i), AlO1molfrac(i), MgOmolfrac(i), FeOmolfrac(i), NiOmolfrac(i), SiO2molfrac(i));
    % 
    % Camol = Ca(i);
    % Almol = Al(i);
    % Femol = Fe(i);
    % Mgmol = Mg(i);
    % Nimol = Ni(i);
    % Omol = O(i);
    % Simol = Si(i);
    % 
    % py_module = py.importlib.import_module('parfit_function');
    % py_parfit_function = py.getattr(py_module, 'parfit_function');
    % [logK_matrix] = py_parfit_function(Padi(i), Tadi(i), Omolfracpick, Femolfrac, Simolfrac, Omol, Femol, Almol, Camol, Mgmol, Simol,frac, He_Conc_metal);
    % 
    % logkmatrix = logK_matrix;
    % logkmatrix = double(logkmatrix)
    % 
    % 
    % Opart_fit = logkmatrix(1,2)
    % Sipart_fit = logkmatrix(2,2)
    % Upart_fit = logkmatrix(3,2)
    % Thpart_fit = logkmatrix(4,2)
    
    % % 
    % %Log 3He/3He at end of 4.5 and U+th CC at end.
    % [helium_temporal, Dm_U, Dm_Th, U_left, Th_left, He3_ratio, helium_temp_0_Gya, helium_temp_4_Gya] = Uranium_thorium_6_chidester(He_Conc_metal, Tadi, Padi, frac, Upart_fit, Thpart_fit);
    % 
    % %[helium_temporal,Dm,U_left, Th_left, He3_ratio, helium_temp_0_Gya, helium_temp_4_Gya] = Uranium_thorium_6(He_Conc_metal, Tadi, Padi, FeOmolfrac(i), Femolfrac(i), Omolfracpick(i), dIW(i), SiO2molfrac(i), frac); 
    % 
    % %[helium_temporal,Dm,U_left, Th_left, He3_ratio, helium_temp_0_Gya, helium_temp_4_Gya] = Uranium_thorium_6_old(He_Conc_metal,  Tadi, dIW(i), SiO2molfrac(i), frac);
    % 
    % helium_temp_0_Gya_logger(i) = helium_temp_0_Gya;
    % helium_temp_4_Gya_logger(i) = helium_temp_4_Gya;
    % U_left_logger(i) = U_left(i);
    % Th_left_logger(i) = Th_left(i);
    % He3_ratio_logger(i) = He3_ratio(i);
    % 
    % DIW_LOG = dIW_logger;
    % 
    % Dm_U_logger = Dm_U;
    % Dm_Th_logger = Dm_Th;
    % 
    % 

    D = [A;B];
    count=count+1
end

figure(1); hold on
plot(r_v_logger, Tb, 'k', 'linewidth', 2)
xlabel('Radius, Km')
ylabel('Temperature (K)')
title('Temp, K')
legend(sprintf('Temperature\nMainMO(.2, 4200, 6000, 50, 7,-2, 1,.25, .4)'))
hold off
ca = gca;
ca.FontSize = 20;

% 
% figure(2)
% plot(r_v_logger, Pressure_bottom_pa/101300, 'linewidth', 2)
% xlabel('Radius, Km')
% ylabel('Pressure Bottom, Atmospheres')
% title('Pressure as radius increases')
% legend(sprintf('Pressure Bottom\nMainMO(.2, 4200, 6000, 50, 7,-2, 1,.25, .4)'))
% hold off
% ca = gca;
% ca.FontSize = 20;

figure(3); hold on
plot(r_v_logger, He_atm_grams_logger, 'b', 'linewidth', 1.5)
plot(r_v_logger, He_Bse_grams_logger, 'r', 'linewidth', 1.5)
plot(r_v_logger, He_Core_grams_logger, 'k', 'linewidth', 1.5)
hold off
xlabel('Radius')
ylabel('Helium Weight Grams')
title('Helium Distribution Grams')
legend('Atmosphere Grams', 'Magma Ocean Grams', sprintf('Core\nMainMO(.2, 4200, 6000, 50, 7,-2, 1,.25, .4)'))
ca = gca;
ca.FontSize = 20;
figure(4)
plot(r_v_logger, He_Core_grams_logger, 'k', 'linewidth', 1.5)
xlabel('Radius')
ylabel('Core Helium Weight Grams')
title('Helium Distribution Grams')
legend(sprintf('Core Grams\nMainMO(.2, 4200, 6000, 50, 7,-2, 1,.25, .4)'))
ca = gca;
ca.FontSize = 20;

% figure(5)
% plot(r_v_logger, He_sol_g_g_logger, 'linewidth', 2)
% xlabel('Radius, Km')
% ylabel('He Solubility g/g')
% title('Solubility as Radius increases')
% legend(sprintf('He solubility\nMainMO(.2, 4200, 6000, 50, 7,-2, 1,.25, .4)'))
% ca = gca;
% ca.FontSize = 20;
% 
% 

% 
% % Add text to all figures
% text(0.2, 0.9, 'Silicon Metal: 0.2')
% text(0.2, 0.8, 'Radius Start: 4000 Km')
% text(0.2, 0.7, 'Radius End: 6000 Km')
% text(0.2, 0.6, 'Time of Accretion: 1e7 Myr')
% text(0.2, 0.5, 'Helium Partitioning Coefficient: -2')
% text(0.2, 0.4, 'Magma Ocean Equilibration Deph: 0.25')
% text(0.2, 0.3, 'Adiabat: 0.4 K/Km')




% 

% columnTitles = {'radius', 'Tadi','Padi', 'Helium_Temp_4_Gya_logger','DIW', 'Dm_logger'};
% IPNB_1e6 = cat(7, Radius, Tadi, Pressure_adiabat, helium_temp_4_Gya_logger,DIW_LOG, Dm_logger);
% dataMatrix = reshape(IPNB_1e6, [], numel(columnTitles));
% 
% % Combine column titles and data into a cell array
% outputData = [columnTitles; num2cell(dataMatrix)];
% 
% % Save the data to a CSV file
% %writecell(outputData, 'C:\Users\Osori\Documents\MATLAB\Masters_code\output.csv');
% writecell(outputData, '/Users/marv/Library/CloudStorage/Box-Box/Masters_code/outputs/output.csv');
% 
% % % Convert data to a table 
% % dataTable = array2table(dataMatrix, 'VariableNames', columnTitles); 
% % % save CSV
% % writetable(dataTable,'C:\Users\Osori\Documents\MATLAB\Current code\ip_im_check_taccr_1e6.csv')

end



