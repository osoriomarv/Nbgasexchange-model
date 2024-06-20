# Set the path to the MATLAB executable
$matlab = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe"

# Set the path to the MATLAB function you want to call
$function = "C:\Users\Osori\Documents\MATLAB\Masters_code\Main_MO.m"

# Set the path to the CSV file containing input parameters
$inputFile = "C:\Users\Osori\Documents\MATLAB\Masters_code\100_check_randomparameters.csv"

# Set the path to the folder to save the output files
$outputFolder = "C:\Users\Osori\Documents\MATLAB\Masters_code\success_params"

# Create the output folder if it doesn't exist
New-Item -ItemType Directory -Force -Path $outputFolder | Out-Null

# Read the input parameters from the CSV file
$inputParams = Import-Csv $inputFile

# Define the titles for the output CSV file
$outputTitles = "SiO2molfrac", "AlO1molfrac", "FeOmolfrac", "MgOmolfrac", "CaOmolfrac", "T_v_K", "Tadi", "Pressure_bottom_pa", "Padi", "U_left_logger", "Th_left_logger", "He3_ratio_logger", "Dm_U_logger", "Dm_Th_logger", "DIW_LOG", "core_grams", "Bse_grams", "He_atm_grams", "helium_temp_0_Gya_logger", "helium_temp_4_Gya_logger", "PPM_TH_core", "PPM_U_core", "He_Conc_metal_logger", "He_Core_grams_logger", "He_Bse_grams_logger"

# Loop through each row and output the corresponding values
foreach ($param in $inputParams) {
    $Set = $param.Set
    $psi_si = $param.psi_si
    $D_he = $param.d_he
    $Taccretion = $param.taccr
    $frac = $param.frac
    $adiabat = $param.adiabat
    $radius = $param.radius
    $radius2 = $param.radius2
    $step = $param.step
    $moequi = $param.moequi
    
    # Set the filename for the output file
    $outputFile = Join-Path $outputFolder "outputs_1000chirandom_$([int]$Set).csv"

    # Build the command to call the MATLAB function and save the output to the CSV file
    $command = "& `"$matlab`" -nodisplay -nosplash -nodesktop -r `"try, cd('C:\Users\Osori\Documents\MATLAB\Masters_code'); [SiO2molfrac, AlO1molfrac, FeOmolfrac, MgOmolfrac, CaOmolfrac, T_v_K, Tadi, Pressure_bottom_pa, Padi, U_left_logger, Th_left_logger, He3_ratio_logger, Dm_U_logger, Dm_Th_logger, DIW_LOG, core_grams, Bse_grams, He_atm_grams, helium_temp_0_Gya_logger, helium_temp_4_Gya_logger, He_Conc_metal_logger, He_Core_grams_logger,He_Bse_grams_logger] = Main_MO($psi_si, $radius, $radius2, $step, $Taccretion,$D_he,$frac,$adiabat,$moequi); 

    results = [$Set;SiO2molfrac; AlO1molfrac; FeOmolfrac; MgOmolfrac; CaOmolfrac; T_v_K; Tadi; Pressure_bottom_pa; Padi; U_left_logger; Th_left_logger; He3_ratio_logger; Dm_U_logger; Dm_Th_logger; DIW_LOG; core_grams; Bse_grams; He_atm_grams; helium_temp_0_Gya_logger; helium_temp_4_Gya_logger;He_Conc_metal_logger;He_Core_grams_logger;He_Bse_grams_logger]; writematrix(results.','${outputFile}'); catch ME, disp(getReport(ME)); end; exit;`""

    # Add titles to the output CSV file
    $outputData = $outputTitles -join ","
    $outputData | Out-File -FilePath $outputFile -Encoding UTF8 -Append

    # Append the results to the output CSV file
    $resultsData = $results -join ","
    $resultsData | Out-File -FilePath $outputFile -Encoding UTF8 -Append
    
    # Execute the command
    Invoke-Expression $command
}
