% Rewrite Code there is an error in melt_mols_i

function [n_SiO2, n_MgO, n_FeO, n_Al2O3, n_CaO, n_NiO, n_Fe, n_Ni, n_O_metal, n_Si_metal] = enstatitecomp_2(psi_si)
%Initial Composition
Si_mg_g = 167;
Mg_mg_g = 106;
Al_mg_g = 8.1;
Ca_mg_g = 8.5;
Fe_mg_g = 290;
Ni_mg_g = 17.5;
O_mg_g = 263;

% Calculate the molar amounts of MnO, CrO, P2O5, and Na2O
Mn_mg_g = 2.2;
Cr_mg_g = 3.15;
P_mg_g = 2;
Na_mg_g = 6.8;

% Define molar masses (in g/mol)
M_Si = 28.0855;
M_Mg = 24.305;
M_Al = 26.9815;
M_Ca = 40.078;
M_Fe = 55.845;
M_Ni = 58.6934;
M_O = 15.999;
M_P = 30.9738;
M_Na = 22.9898;
M_Mn = 54.9;
M_Cr = 52;

% Convert from mg/g to mols
Si = Si_mg_g / M_Si;
Mg = Mg_mg_g / M_Mg;
Al = Al_mg_g / M_Al;
Ca = Ca_mg_g / M_Ca;
Fe = Fe_mg_g / M_Fe;
Ni = Ni_mg_g / M_Ni;
O = O_mg_g / M_O;

n_MnO = Mn_mg_g / (M_Mn);
n_CrO = Cr_mg_g / (M_Cr);
n_P2O5 = P_mg_g / (2 * M_P);
n_Na2O = Na_mg_g / (2 * M_Na);

% Subtract the molar amounts of MnO, CrO, P2O5, and Na2O from the total O
n_O_remaining = O - n_MnO - n_CrO - 5 * n_P2O5 - n_Na2O;

% Calculate the molar amounts of MgO, Al2O3, and CaO (using all Mg, Al, and Ca)
n_MgO = Mg;
n_Al2O3 = Al / 2;
n_CaO = Ca;
% Calculate the molar amount of remaining O after forming MgO, Al2O3, and CaO
n_O_remaining = n_O_remaining - n_MgO - 3 * n_Al2O3 - n_CaO;

% Calculate the molar amount of SiO2 using the remaining O
n_SiO2 = min(Si, n_O_remaining / 2); % each SiO2 requires 1 mole of Si and 2 moles of O
% Calculate the molar amount of remaining Si after forming SiO2
n_Si_remaining = Si - n_SiO2;

% The remaining Fe, Ni, Si, and O go into the metal component
n_Fe = Fe;
n_Ni = Ni;
n_Si_metal = max(n_Si_remaining, 0); % ensure that n_Si_metal is non-negative
n_O_metal = n_O_remaining - 2 * n_SiO2;


melt_mols_i = [n_MgO, n_Al2O3, n_CaO, n_SiO2];
melt_mol_frac_i=melt_mols_i./sum(melt_mols_i);
metal_mols_i = [n_Fe, n_Ni, n_O_metal, n_Si_metal];
metal_mol_frac_i=metal_mols_i./sum(metal_mols_i);


% Fixed - Move to windows
si_metal_ratio_check = metal_mols_i(4)/(melt_mols_i(4)+metal_mols_i(4));
% Calculate the observed Si-metal ratio
observed_ratio = psi_si;

% Define the tolerance for accepting the solution
tolerance = 0.001; %1
% Define the range of O_mg_g values to test
O_mg_g_range = linspace(100, 400, 10000);
% Calculate the chi-square values for each O_mg_g value in the range
chi_squared_values = zeros(size(O_mg_g_range));

for i = 1:length(O_mg_g_range)
    % Update O_mg_g
    O_mg_g = O_mg_g_range(i);
    O = O_mg_g / M_O;
    n_O_remaining = O - n_MnO - n_CrO - 5 * n_P2O5 - n_Na2O;
    % Calculate the molar amounts of each component
    n_MgO = Mg;
    n_Al2O3 = Al/ 2;
    n_CaO = Ca;
    n_O_remaining = n_O_remaining - n_MgO - 3 * n_Al2O3 - n_CaO;
    n_SiO2 = min(Si, n_O_remaining / 2);
    n_Si_remaining = Si - n_SiO2;
    n_Si_metal = max(n_Si_remaining, 0);
    n_O_metal = n_O_remaining - 2 * n_SiO2; % each SiO2 requires 2 oxygen atoms
    n_Fe = Fe;
    n_Ni = Ni;
    % Calculate the Si-metal ratio
    melt_mols_i = [n_MgO, n_Al2O3, n_CaO, n_SiO2];
    metal_mols_i = [n_Fe, n_Ni, n_O_metal, n_Si_metal];
    si_metal_ratio = metal_mols_i(4) / (melt_mols_i(4) + metal_mols_i(4));
    expected_ratio = si_metal_ratio;
    % Check if the solution is within the tolerance
    if abs(si_metal_ratio - observed_ratio) / observed_ratio < tolerance
        disp(['Solution found for O_mg_g = ' num2str(O_mg_g) ' with Si-metal ratio = ' num2str(si_metal_ratio)])
        break
    end
    % Calculate the chi-square value
    chi_squared_values(i) = (expected_ratio - observed_ratio)^2;
end

% Find the index of the minimum chi-square value
[min_chi_squared_value, min_chi_squared_index] = min(chi_squared_values);

% Get the corresponding value of O_mg_g
O_mg_g_best = O_mg_g_range(min_chi_squared_index);

n_FeO =0;
n_NiO =0;

if n_O_metal == 0
    n_O_metal  = 1e-10;
else
    n_O_metal = n_O_metal;
end
if n_Si_metal == 0
    n_Si_metal = 1e-10;
else
    n_Si_metal = n_Si_metal;
end

% Report the results
disp(['MgO: ' num2str(n_MgO)])
disp(['Al2O3: ' num2str(n_Al2O3)])
disp(['CaO: ' num2str(n_CaO)])
disp(['SiO2: ' num2str(n_SiO2)])
disp(['Fe: ' num2str(n_Fe)])
disp(['Ni: ' num2str(n_Ni)])
disp(['Si (metal): ' num2str(n_Si_metal)])
disp(['O (metal): ' num2str(n_O_metal)])

%melt_mols_i = [SiO2, MgO, FeO, Al2O3, CaO, NiO];
%metal_mols_i_scaled=[Fe, Ni, O, Si]*frac;

% melt_mols = [4.762, 4.3612, 0, .45031, 0.21209, 0];
% metal_mols_i_scaled = [5.1929, .29816, 1e-10, 1.1901];
end


