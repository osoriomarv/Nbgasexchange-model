%URanium oxygen partitioning 

% Equation 24
function log_du = calculate_log_du(T, P, mol_O_metal, mol_U)
    % Constants
    a = -0.55;
    b = -10133;
    c = 0;
    E_O_i = -31;

    % Calculate Dfe
    Dfe = 1; % You need to define this function separately

    % Calculate Ln_gamma_U
    Ln_gamma_U =  E_O_i * mol_O_metal *(log(1-mol_O_metal)/mol_O_metal) + ...
                    E_O_i * mol_O_metal.^2 * mol_U *(1/(1-mol_O_metal) + mol_U/2*(1-mol_U)^2) ;

    % Calculate log_gamma_u_metal
    log_gamma_u_metal = Ln_gamma_U / 2.303;

    % Calculate Log_K_U
    Log_K_U = a + b/T + c*(P/T);

    % Calculate Log_Du
    log_du = Log_K_U + log(Dfe) - log_gamma_u_metal;
end

%Mas
%first line is 0

%Equation 4 = log(du) = a + b + c + log(Dfe) - log_gamma_uranium_metal

%

%I want to write a code that involves these equations using these equations
%which should be a matlab function
% a = -.55
% b = -10133
% c = 0
% E_o_i = -31;


%T and P will be input while log_du will be outputted;


%Equation 4 from "The solubility of heat-producing elemetns in Earth's
%core"
%Log_K_U = Log(Du/Dfe) + log(gamma_u_metal/gamma_Fe_Metal) = a + b/T +
%c*(P/T)

%find log_gamma_u_metal
    
    %Ln_gamma_U = - sigma_k to N(E_o_i * mol_O_metal *(1 + (ln(1-mol_O_metal)/mol_O_metal) - (1/(1-mol_U)) + 
                %sigma_K to N(E_o_i * mol_O_metal^2 *mol_U *((1/(1-mol_u) +
                %(1/(1-mol_U) + 1/(1-mol_o_metal) + mol_U/2*(1-mol_U)^2 - 1)
    
    %log_gamma_u_metal = ln_gamma_u / 2.303;
    
    % Dfe = x_metal / x_sil

%Log_(Du) = a + b + c + log(Dfe) - log_gamma_u_metal)




