function [u, p] = FBAlocalsens(model1, param, factor, inc)

% This function calculates the growth rate as a function of all parameters
% fixed in the model
% INPUTS : model1 : a COBRA model
%          param : vector of parameters containing 67 numeric values
%          factor: fold change for the sensitivity analysis
%          inc : for two parameters with an initial value of 0, inc = the
%          upper bound of the parameter value. 0 is the lower bound
% OUTPUT : A matrix with 67 rows (for each parameter) and 2 columns. The
% first column gives the growth rate with a decrease by a given factor of
% the parameter and the second column also gives the growth rate for the parameter in row,
% but for an increase by a given factor.


i = 0;
j = 0;
u = zeros(length(param),2);
param_ini = param;
for j = 1:2
    for i = 1:length(param)

        param = param_ini;

        if param(i) == 0 && j == 1 % i == 30 || i == 33 || i == 56 && j == 1
            param(i) = param(i) ; % the values of param(30) and param(33) and param(56) are zero, no changes
        end

        if param(i) == 0 && j == 2  % i == 30 || i == 33 || i == 56 && j == 2
            param(i) = param(i) + inc ; 
        end
        
        
        if j == 1 
            param(i) = param(i) / factor;
        else
            param(i) = param(i) * factor;
        end


% Update the biomass_pigm_h equation :
% All pigment concentration are from Guérin et al MSc thesis, Alderkamp et
% al 2011 assuming a chlc2:chlc1 ratio of 2.2 (Owens and Wold, 1986)
% Values of pigment concentrations were converted in stoichiometric
% coefficients

% Pigments measured in continuous light (Guérin et al, in prep)
Ccell = param(68); % 7.5E-12;   % g C / cell
DWcell = param(69); % Ccell * 2.;  % Assume a g DW / g C = 2
chla = param(1); %0.14          % pg Chla / cell
chlc = param(2); % 0.03           % pg Chlc / cell
Fucoxan = param(3); % 0.06        % pg Fucoxanthin / cell
B_carot = param(4); %0.004       % pg B-carotene / cell
DD = param(5); % 9/1000        % pg DD / cell   ???? Alderkamp.... ??

% Pigments in percentage of mass per g DW (% g / gDW)
chla = chla * 1E-12 / DWcell * 100.; 
Fucoxan = Fucoxan * 1E-12 / DWcell * 100.;
B_carot = B_carot * 1E-12 / DWcell * 100.;
DD = DD * 1E-12 / DWcell * 100.;
chlc = chlc * 1E-12 / DWcell * 100.;
chlc_rat = 2.2;    % Ratio chlorophyll c2 : chlorophyl c1 in P. tricornutum (Owens & Wold (1986) Plant Physiol 80: 732?738)		
chlc1 = 1/chlc_rat * chlc;
chlc2 = (1-(1/chlc_rat)) * chlc;

% Molar mass of components (g/mol)
molChla =	892.4782;
molFucoxan =	658.9040;
molB_carot =	536.8704;
molDD = 582.8528;
molChlc1 =	607.9166;
molChlc2 =	609.9324;

% Pigment composition (mol g DW-1)
chla = chla / molChla ;
Fucoxan = Fucoxan / molFucoxan;
B_carot = B_carot / molB_carot;
DD = DD / molDD;
chlc1 = chlc1 / molChlc1;
chlc2 = chlc2 / molChlc2;

% Non-normalized biomass weight (g / g DW)
NonNormBiom = (chla * molChla) + (Fucoxan * molFucoxan) + (B_carot * molB_carot) + (DD * molDD) + (chlc1 * molChlc1) + (chlc2 * molChlc2);
Correction = NonNormBiom  ;

% Stoichiometric coefficients (mmol / g DW)
chla = chla / Correction * 1000.;
Fucoxan = Fucoxan / Correction * 1000.;
B_carot = B_carot / Correction * 1000.;
DD = DD / Correction * 1000.;
chlc1 = chlc1 / Correction * 1000.;
chlc2 = chlc2 / Correction * 1000.;

printRxnFormula(model1, 'biomass_pigm_h');
model1 = addReaction(model1, 'biomass_pigm_h', 'metaboliteList', {'cholphya_h',...
    'cholphyc1_h', 'cholphyc2_h', 'caro_h', 'fxanth_h',...
    'diadinx_h', 'biomass_pigm_h'},...
    'stoichCoeffList', [-chla; -chlc1; -chlc2; -B_carot; -Fucoxan;...
    -DD; 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

% Update biomass_DNA_c for the total F. cylindrus genome size and its low
% G:C content (Mock et al 2017)
GC = param(6); %0.398    % Proportion of total nucleotides as GC (Mock et al 2017)
AT = 1-GC;
GenomeSize = param(7); % 61.1E+6   % Total genome size (number of bases)

% Calculating total number of nucleotides per cell (dAMP, dTMP, dGMP, dCMP)
dNucl = [GenomeSize * AT / 2, GenomeSize * AT / 2, GenomeSize * GC / 2, GenomeSize * GC / 2];

% Molar mass of nucleotides (mol/g) (moldAMP, moldTMP, moldGMP, moldCMP)
molMassNucl = [329.2055, 320.1921, 345.2049, 305.1808];

% Nucleotides (g / cell)
Avog = 6.02E+23;     % Avogadro constant in atom / mol
NuclMass = dNucl .* molMassNucl / Avog;

% Relative nucleotide and pyrophosphate abundance (number or mol/gDW)
Rel_dNucl = dNucl ./ sum(dNucl);
Rel_Pyr = 1;

% triphosphate nucleotide (dNTP) and pyrophosphate molar mass (g/mol)
% (dATP, dTTP, dGTP, dCTP, pyrophosphate or ppi)
MolMassdNTP = [487.1495, 478.1361, 503.1489, 463.1248];
MolMassPyr = 174.9513;

% Non-normalized biomass weight (g / gDW = gDW / gDW)
NonNormBiom = (sum(Rel_dNucl .* MolMassdNTP)) - (Rel_Pyr * MolMassPyr);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (mmol / gDW) of dATP, dTTP, dGTP, dCTP (S_dNTP) and Pyr
% (S_Pyr)
dNTP = Rel_dNucl ./ Correction;
Pyr = Rel_Pyr / Correction;

printRxnFormula(model1, 'biomass_DNA_c');
model1 = addReaction(model1, 'biomass_DNA_c', 'metaboliteList', {'dgtp_c',...
    'dttp_c', 'dctp_c', 'datp_c', 'ppi_c', 'biomass_dna_c'},...
    'stoichCoeffList', [-dNTP(3); -dNTP(2); -dNTP(4); -dNTP(1); Pyr; 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

% Update 'biomass_RNA_c' assuming a RNA:DNA ratio of 8 (Lourenco et al (1998) J Phycol 34: 798?811)		
% RNA Nucleotides (g / cell) (dAMP, dTMP, dGMP, dCMP)
RNA_DNA_r = param(8); % 8.
Mass_Nucl = NuclMass * RNA_DNA_r;

% Relative abundance is the same than for DNA
Rel_Nucl = Rel_dNucl;
Rel_Pyr;

% Molar mass of RNA nucleotides (g/mol) ATP, UTP, GTP, CTP, Pyr
MolMassNTP = [503.1489, 480.1090, 519.1483, 479.1242];
MolMassPyr = 174.9513;

% Non-normalized biomass weight
NonNormBiom = (sum(Rel_Nucl .* MolMassNTP)) - (Rel_Pyr * MolMassPyr);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (ATP, UTP, GTP, CTP, Pyr)
NTP = Rel_Nucl ./ Correction;
Pyr = Rel_Pyr / Correction;

printRxnFormula(model1, 'biomass_RNA_c');
model1 = addReaction(model1, 'biomass_RNA_c', 'metaboliteList', {'atp_c',...
    'gtp_c', 'utp_c', 'ctp_c', 'ppi_c', 'biomass_rna_c'},...
    'stoichCoeffList', [-NTP(1); -NTP(3); -NTP(2); -NTP(4); Pyr; 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

% Update biomass_mem_lipids_c equation
% Molar mass of lipids in g/mol ()
molMassLip = [768, 722, 724, 726, 936, 765, 903, 793, 1123, 764.97,...
    780.97, 824.97, 728.97, 803, 782.97, 806.97, 616, 712, 568];

% Molar ratios of lipids (mol/ gDW)
MolRatLip = param(9:19+8);
%[0.0555, 0.1244, 0.0622, 0.0504, 0.0450, 0.1043, 0.0318, 0.0144,...
%    0.0070, 0.0490, 0.0044, 0.0622, 0.0670, 0.0380, 0.0180, 0.0220,...
%    0.0104, 0.0062, 0.0088]

% Mass ratio (g / gDW = gDW / gDW)
MassRatLip = molMassLip .* MolRatLip;

% Non-normalized biomass weight (g / gDW = gDW / gDW)
NonNormBiom = sum(MassRatLip);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients of each lipid (mmol / gDW)
S_Lip = MolRatLip ./ Correction;

printRxnFormula(model1, 'biomass_mem_lipids_c');
model1 = addReaction(model1, 'biomass_mem_lipids_c', 'metaboliteList', {'mgdg205n3164n1_h',...
    'mgdg1619Z163n4_h', 'mgdg1619Z162n4_h', 'mgdg161_h', 'dgdg205n31619Z_h',...
    'sqdg140160_h', 'sqdg1619Z240_c', 'sqdg160_h', 'asq205n31602o205n3_h',...
    'pg205n31613E_h', 'pc182_9_12_c', 'pc205n3_c', 'pc161_c',...
    'dgta205n3_c', 'pe205n3_c', 'pail1619Z160_c', '12dgr182_9_12_c',...
    '12dgr226n3_c', '12dgr160_h', 'biomass_mem_lipids_c'},...
    'stoichCoeffList', [-S_Lip(1); -S_Lip(2); -S_Lip(3); -S_Lip(4); -S_Lip(5);...
    -S_Lip(6); -S_Lip(7); -S_Lip(8); -S_Lip(9); -S_Lip(10); -S_Lip(11); -S_Lip(12);...
    -S_Lip(13); -S_Lip(14); -S_Lip(15); -S_Lip(16); -S_Lip(17); -S_Lip(18); -S_Lip(19); 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

% Update the 'biomass_pro_c' equation

% Weight amino acid percentage (Brown, 1991). The sum = 100%
% Amino acid list : 'alatrna_c, argtrna_c, asntrna_c, asptrna_c, cystrna_c,
% glntrna_c, glutrna_c, glytrna_c, histrna_c, iletrna_c, leutrna_c,
% lystrna_c, mettrna_c, phetrna_c, protrna_c, sertrna_c, thrtrna_c, trptrna_c,
% tyrtrna_c, valtrna_c
MassPercent = param(28:47);
%[7.2, 6.6, 0, 8.6, 0.38, 0, 11.2, 5.8, 1.6, 4.9, 7.7, 5.6,...
%    1.9, 6.6, 6.3, 5.9, 5.4, 1.6, 4.1, 5.6]

% Normalized weight (% g) : the total = 1
NormWeight = MassPercent / sum(MassPercent);
NormWeight(3) = MassPercent(4) / sum(MassPercent) / 2.;
NormWeight(6) = MassPercent(7) / sum(MassPercent) / 2.;

% Amino acid molar mass (mol/g)
MolMassAmino = [73.0935, 159.2089, 116.1182, 116.0951, 105.1585, 130.1447,...
    130.1216, 59.0670, 139.1548, 115.1730, 115.1730, 131.1955, 133.2115,...
    149.1893, 99.1307, 89.0929, 103.1194, 188.2253, 165.1887, 101.1465];

% Molar Mass of other reactants (g/mol) (Water, atp_c, gtp_c)
MolMassReac = [18.0152, 503.1489, 519.1483];

% Molar Mass of other products (g/mol) (Proton, Phosphate, GDP, ADP,
% Unloaded tRNA)
MolMassProd = [1.0079, 95.9793, 440.1763, 424.1769, 1.0079];

% Molar ratio (% mol/gDW)
MolRatAmino = NormWeight ./ MolMassAmino * 100.;

% Molar ratio of other reactants (mol/ gDW) (Water, atp_c, gtp_c)
SumMolRat = sum(MolRatAmino);
MolRatReac = [3, 1, 2] * SumMolRat;

% Molar ratio other product (mol/ gDW) (Proton, Phosphate, GDP, ADP,
% Unloaded tRNA) where unloaded tRNA accounts for the sum of all molar
% ratio of all amino acids
MolRatProd = [4, 3, 2, 1, 1] * SumMolRat;

% Non-normalized biomass weight (g DW / gDW)
NonNormBiom = sum(NormWeight) * 100. + (sum(MolMassReac .* MolRatReac)) - (sum(MolMassProd .* MolRatProd));
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients of amino acids (mmol / gDW)
S_Amino = - MolRatAmino / Correction;
Sum_S_Amino = - sum(S_Amino);

% Stoichiometric coefficients of other reactants (mmol / gDW)
% Water, ATP, GTP
S_Reac = - [3, 1, 2] * Sum_S_Amino;

% Stoichiometric coefficients of other products (mmol / gDW)
S_Prod = [4, 3, 2, 1] * Sum_S_Amino;

% Stoichiometric coefficients (mmol / gDW) of all 'tra-amino acid species' as products
% 'trnaala_c', 'trnaarg_c'
% But those coefficients are positive...
S_Amino_Prod = - S_Amino;

printRxnFormula(model1, 'biomass_pro_c');
model1 = addReaction(model1, 'biomass_pro_c', 'metaboliteList', {'h2o_c',...
    'atp_c', 'gtp_c', 'alatrna_c', 'argtrna_c', 'asntrna_c',...
    'asptrna_c', 'cystrna_c', 'glntrna_c', 'glutrna_c',...
    'glytrna_c', 'histrna_c', 'iletrna_c', 'leutrna_c', 'lystrna_c',...
    'mettrna_c', 'phetrna_c', 'protrna_c', 'sertrna_c', 'thrtrna_c',...
    'trptrna_c', 'tyrtrna_c', 'valtrna_c',...  % End of reactants
    'h_c', 'pi_c', 'adp_c', 'gdp_c', 'trnaala_c', 'trnaarg_c', 'trnaasn_c',...
    'trnaasp_c', 'trnacys_c', 'trnagln_c', 'trnaglu_c', 'trnagly_c',...
    'trnahis_c', 'trnaile_c', 'trnaleu_c', 'trnalys_c', 'trnamet_c',...
    'trnaphe_c', 'trnapro_c', 'trnaser_c', 'trnathr_c', 'trnatrp_c',...
    'trnatyr_c', 'trnaval_c', 'biomass_pro_c'},...
    'stoichCoeffList', [S_Reac(1); S_Reac(2); S_Reac(3); S_Amino(1); S_Amino(2); S_Amino(3); S_Amino(4);...
    S_Amino(5); S_Amino(6); S_Amino(7); S_Amino(8); S_Amino(9); S_Amino(10); S_Amino(11);...
    S_Amino(12); S_Amino(13); S_Amino(14); S_Amino(15); S_Amino(16); S_Amino(17); S_Amino(18);...
    S_Amino(19); S_Amino(20); S_Prod(1); S_Prod(2); S_Prod(3); S_Prod(4); S_Amino_Prod(1); S_Amino_Prod(2);...
    S_Amino_Prod(3); S_Amino_Prod(4); S_Amino_Prod(5); S_Amino_Prod(6);...
    S_Amino_Prod(7); S_Amino_Prod(8); S_Amino_Prod(9); S_Amino_Prod(10);...
    S_Amino_Prod(11); S_Amino_Prod(12); S_Amino_Prod(13); S_Amino_Prod(14);...
    S_Amino_Prod(15); S_Amino_Prod(16); S_Amino_Prod(17); S_Amino_Prod(18);...
    S_Amino_Prod(19); S_Amino_Prod(20); 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

% Add an objective function with explicit storage_carbon
% (bof_c_accumulation_c)
% Here I replace 'biomass_tag' by 'carbon_storage_c'
% I replaced 'biomass_c' by 'biomass_c_acc_c'

total_prot = param(48); %param5(1) %25./100       %     ratio in g/gDW
total_carb = param(49); %20./100       %    ratio in g/DW
total_lipid = param(50); % 34./100      %    ratio in g/ gDW
total_DNA = param(51); %0.5406/100     %    ratio in g/gDW assumed from P. tricornutum
total_RNA = param(52); %4.3244 /100    %    ratio in g/gDW assumed from P. tricornutum
total_pigm = param(53); %1.56 / 100    %   ratio in g/gDW measured in F. cylindrus
SumTot = total_prot + total_carb + total_lipid + total_DNA + total_RNA + total_pigm;

total_prot = total_prot * 100.;       %     ratio in % g/gDW
total_carb = total_carb * 100.;      %    ratio in % g/DW
total_lipid = total_lipid * 100. ;    %    ratio in % g/ gDW
total_DNA = total_DNA * 100. ;       %    ratio in % g/gDW calculated for F. cylindrus
total_RNA = total_RNA * 100. ;     %    ratio in % g/gDW calculated for F. cylindrus
total_pigm = total_pigm * 100. ;   %   ratio in % g/gDW measured in F. cylindrus
SumTot = total_prot + total_carb + total_lipid + total_DNA + total_RNA + total_pigm;

% corection of main component is necessary because of error in measurements
% and because of unknown components (such as osmolyte) in the cell
% characterization
% prot_cor = (1.-total_pigm-total_DNA-total_RNA) * total_prot / (total_prot + total_carb + total_lipid)
% carb_cor = (1.-total_pigm-total_DNA-total_RNA) * total_carb / (total_prot + total_carb + total_lipid)
% lip_cor = (1.-total_pigm-total_DNA-total_RNA) * total_lipid / (total_prot + total_carb + total_lipid)

prot_cor = ( total_prot + ( 100. - SumTot )/ 3. ) / 100.;
carb_cor = ( total_carb + ( 100. - SumTot )/ 3. ) / 100.;
lip_cor = ( total_lipid + ( 100. - SumTot )/ 3. ) / 100.;
total_DNA = total_DNA / 100.;
total_RNA = total_RNA / 100.;
total_pigm = total_pigm / 100.;

Prop_glucan = param(54); %10.      % Proportion of cabohydrate as glucan
Prop_TAG = param(55); %10.          % Proportion of lipid accumulated as TAG
mem_lip = lip_cor - (lip_cor * Prop_TAG / 100.);
carb_cor2 = carb_cor - (carb_cor * Prop_glucan / 100.);
TAG = lip_cor * Prop_TAG / 100.;
Glucan = carb_cor * Prop_glucan / 100.;
C_storage = TAG + Glucan;
Mass_bal = prot_cor + carb_cor2 + mem_lip + TAG + Glucan + total_DNA + total_RNA + total_pigm
if Mass_bal < 1.01 & Mass_bal > 0.99
    disp('Mass balance is ok')
else
    error('Mass of all components is unbalanced')   
end

rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;
model1 = addReaction(model1, 'bof_c_accumulation_c', 'metaboliteList', {'biomass_pro_c',...
    'biomass_pigm_h', 'biomass_mem_lipids_c', 'biomass_dna_c',...
    'biomass_rna_c', 'biomass_carb_c', 'carbon_storage_c', 'biomass_c_acc_c'},...
    'stoichCoeffList', [-prot_cor; -total_pigm; -mem_lip; -total_DNA;...
    -total_RNA; -carb_cor2; -C_storage; 1],...
    'reversible', false);

% Add a reaction 'biomass_c_storage_c'
% Molar ratios (mmol / g DW) (13glucan_c, tag1619Z1619Z160)
% amount accumulate per time step...
%MolRat = [0.117, 0.05]

% Molar mass (g/mol)
Molmass = [419.97, 802];

% Mass ratio (g / gDW)
%MassRat = Molmass .* MolRat / 1000.

% Here we use Glucan and TAG calculated in the mass balance of the
% objective function
MassRat = [Glucan, TAG];

MolRat = MassRat ./ Molmass * 1000; % mmol / gDW

% NonNormalized biomass weight (g / gDW)
NonNormBio = sum(MassRat);
Correction = NonNormBio;

% Stoichiometric coefficients (mmol / g)
S_Stor = MolRat / Correction;

model1 = addReaction(model1, 'biomass_c_storage_c', 'metaboliteList', {'13glucan_c',...
    'tag1619Z1619Z160_c', 'carbon_storage_c'},...
    'stoichCoeffList', [-S_Stor(1); -S_Stor(2); 1],... 
    'reversible', false);

% Update the biomass_carb_c equation

% Molar ratio (mol/ gDW) of the 9 sugar species
% Glucose, galactose, mannose, Xylulose, arabinose, fucose, rhamnose,
% glucuronic acid, mannose-sulfate
sugar_mol = param(56:64);
    
% Molar mass of sugars (g/mol)
MolMassSug = [564.2851 564.2851 603.3244 534.2592 534.2592 587.3250 548.2857 577.2608 682.3797];

% Mass ratio of sugars (g / g)
MassRatSug = sugar_mol .* MolMassSug;

% Product molar ratio (mol / gDW); 'Proton, UDP, GDP'
Proton_mol = sum(sugar_mol);
UDP_mol = sum(sugar_mol([1,2,4,5,7,8]));
GDP_mol = sum(sugar_mol([3,6,9]));
Prod_mol = [Proton_mol, UDP_mol, GDP_mol];

% Molar mass products (g/mol)
MolMassProd = [1.0079, 401.1370, 440.1763];

% Mass ratio products (g / g)
MassRatSug2 = Prod_mol .* MolMassProd;

% Non-normalized biomass weight (g / gDW)
NonNormBiom = sum(MassRatSug) - sum(MassRatSug2);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (mmol / gDW)
Stoich_Sug = [sugar_mol, Prod_mol];
S_Sug = Stoich_Sug / Correction;

printRxnFormula(model1, 'biomass_carb_c');

model1 = addReaction(model1, 'biomass_carb_c', 'metaboliteList', {'udpg_c',...
    'udpgal_c', 'gdpmann_c', 'udpxyl_c',...
    'udparab_c', 'gdpfuc_c', 'udprmn_c', 'udpglcur_c', 'gdpman2s_c',...
    'h_c', 'udp_c', 'gdp_c', 'biomass_carb_c'},...
    'stoichCoeffList', [-S_Sug(1); -S_Sug(2); -S_Sug(3); -S_Sug(4); -S_Sug(5);...
    -S_Sug(6); -S_Sug(7); -S_Sug(8); -S_Sug(9); S_Sug(10); S_Sug(11);...
    S_Sug(12); 1],...
    'reversible', false);
rxns_length1 = length(model1.rxns);
mets_length1 = length(model1.mets) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Objective function
% model1 = changeObjective(model1, 'DM_biomass_c_acc_c')
Ps = param(67);
model1 = changeRxnBounds(model1, 'EX_hco3_e', -Ps, 'b'); % mmol hco3 g DW-1 h-1

C_N = param(65);
N_P = param(66);
fold_change = 1.0;
N_low =  - Ps / (C_N * fold_change) ; % in mmol / gDW / h
N_high = - Ps / (C_N / fold_change);
P_low = - Ps / (C_N * N_P * fold_change);  % in mmol / gDW / h 
P_high = - Ps / (C_N * N_P / fold_change);
model1 = changeRxnBounds(model1,'EX_no3_e', N_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_no3_e', N_low, 'u'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_low, 'u'); % Bd: -1000 / 0

FBAsolution = optimizeCbModel(model1, 'max','one');
u(i,j) = FBAsolution.f * 24;

%u(i,j) = FBAsolution.x(2143) * 24;
% Rxn IDs of biomass function : 2100 to 2105 , 2113 + 2114, 2119, 2140, 2141, 2143, 2144-2146
%model1.rxns( find(contains(model1.rxns(), 'biomass') ) )
%findRxnIDs(model1, ans) % 2100, 2101, 2102, 2110, 2111, 2116, 2143
%p(i,j) = 2
%pigmbio(i, j) = FBAsolution.x(2101)

%u(i,j) = FBAsolution.x(2143);

    end
end