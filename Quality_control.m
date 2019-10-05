%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quality control of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Load the model
run Model_setup

%%

% Testing if any reactions are possible without C and light input
modelNoC = model1
modelNoC = changeRxnBounds(modelNoC, 'EX_hco3_e', 0, 'b');
modelNoC = changeRxnBounds(modelNoC, 'EX_photon_e', 0, 'b');
modelNoC = changeRxnBounds(modelNoC,'ATPM_h', 0, 'b'); 
modelNoC = changeRxnBounds(modelNoC,'ATPM_m', 0, 'b'); 
modelNoC = changeRxnBounds(modelNoC,'ATPM_c', 0, 'b'); 
FBAsolutionNoC = optimizeCbModel(modelNoC, 'max', 'one')

%printFluxVector(modelNoC, FBAsolutionNoC.x);
%printFluxVector(modelNoC, FBAsolutionNoC.x, true);

if any(FBAsolutionNoC.x ~= 0)
    disp('There is a problem! Matter is produced from nothing');
else
    disp('Good, there are no active reactions when C and light input = 0');
end

% Testing whether or not ATP can be produced from nothing and if any
% reactions in the model can be significant without any exchange reactions

modelClosed = model1;
% set exchange reaction boundaries
[selExc, selUpt] = findExcRxns(model1); %, 0, 0)
uptakes = model1.rxns(selUpt) ;
exchanges = model1.rxns(selExc) ;  
EXrxn = exchanges;

EXub = transpose(zeros(1,43)); % [1000;  1000;  x;  1000;  1000;  1000;  0;  1000; 1000; 1000;  1000;  1000; 1000; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
EXlb = transpose(zeros(1, 43)); %[-1000; -1000; x; -1000; -0.17; -1000; 0; -0.34; -1000; -1000; -1000; -1000;   0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
for i = 1:length(EXrxn);
    index = findRxnIDs(modelClosed, EXrxn{i});
    modelClosed.ub(index) = EXub(i);
    modelClosed.lb(index) = EXlb(i);
end

modelClosed = changeObjective(modelClosed, {'ATPM_c'});
modelClosed = changeRxnBounds(modelClosed,'ATPM_c', 0, 'l'); 
modelClosed = changeRxnBounds(modelClosed,'ATPM_c', 1000, 'u'); 
checkObjective(modelClosed) % check biomass function
modelClosed = changeRxnBounds(modelClosed,'DM_biomass_c_acc_c', 0, 'b'); 
printConstraints(modelClosed, -1000, +1000)
modelClosed = changeRxnBounds(modelClosed,'ATPM_h', 0, 'b'); 
modelClosed = changeRxnBounds(modelClosed,'ATPM_m', 0, 'b'); 
modelClosed = changeRxnBounds(modelClosed,'CEF_h', 0, 'b'); 

FBAsolutionClosed = optimizeCbModel(modelClosed,'max', 'one');
printFluxVector(modelClosed, FBAsolutionClosed.x);
printFluxVector(modelClosed, FBAsolutionClosed.x, true);

if any(FBAsolutionClosed.x ~= 0)
    disp('There is a problem! Matter is produced from nothing');
else
    disp('Good, there are no active reactions whenn all exchange reactions are closed');
end

% Performing a Leak test. A metabolite is leaking if the exchange reaction can carry secretion flux in the closed model (no uptake flux through any exchange reactions is permitted).
modelLeak = model1;
EXrxn = EXrxn;
EXub = repmat(1000,43,1) ;
EXlb = transpose(zeros(1,43)) %[0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
for i = 1:length(EXrxn);
    index = findRxnIDs(modelLeak, EXrxn{i});
    modelLeak.ub(index) = EXub(i);
    modelLeak.lb(index) = EXlb(i);
end
printConstraints(modelLeak, -1000, +1000)
modelLeak = changeRxnBounds(modelLeak,'ATPM_c', 0, 'l');
modelLeak = changeRxnBounds(modelLeak,'ATPM_c', 1000, 'u');
modelLeak = changeRxnBounds(modelLeak,'CEF_h', 0, 'b');
modelLeak = changeRxnBounds(modelLeak,'ATPM_h', 0, 'b');
modelLeak = changeRxnBounds(modelLeak,'ATPM_m', 0, 'b');

[LeakRxns,modelTested,LeakRxnsFluxVector] = fastLeakTest(modelLeak, EXrxn,'false');
if ~isempty(LeakRxns) & ~isempty(LeakRxnsFluxVector)
    disp('There is a problem, at least one metabolite is leaking');
else
    disp('Good, no metabolites are leaking');
end

% Test if the model produces energy from water
modelWater = modelLeak;
modelWater = changeObjective(modelWater,'ATPM_c');
modelWater = changeRxnBounds(modelWater,'ATPM_c',0,'l');
modelWater = changeRxnBounds(modelWater,'ATPM_c',1000,'u');
modelWater = changeRxnBounds(modelWater,'EX_h2o_e',-1,'l');
modelWater = changeRxnBounds(modelWater,'EX_h2o_e',1000,'u');
printConstraints(modelWater, -1000, 1000)

FBAWater = optimizeCbModel(modelWater, 'max', 'one');
% Print non-zero fluxes if any
printFluxVector(modelWater, FBAWater.x, true);

if any(FBAWater.x ~= 0)
    disp('There is a problem!');
else
    disp('Good, no active reactions in the presence of water uptake');
end

% Test if the model produces energy from water and oxygen
modelWaterO = modelLeak;
modelWaterO = changeObjective(modelWaterO,'ATPM_c');
modelWaterO = changeRxnBounds(modelWaterO,'ATPM_c',0,'l');
modelWaterO = changeRxnBounds(modelWaterO,'ATPM_c',1000,'u');
modelWaterO = changeRxnBounds(modelWaterO,'EX_h2o_e',-1,'l');
modelWaterO = changeRxnBounds(modelWaterO,'EX_h2o_e',1000,'l');
modelWaterO = changeRxnBounds(modelWaterO,'EX_o2_e',-1,'l');
modelWaterO = changeRxnBounds(modelWaterO,'EX_o2_e',1000,'l');
FBAWaterO = optimizeCbModel(modelWaterO, 'max', 'one');

if FBAWaterO.stat == 0 | all(FBAWaterO.x == 0)
    disp('Good, no matter is produced from oxygen and water only');
else
    disp('Error, Matter is produced when ATP demand is reversed');
end

% Test if the model produces matter when ATP demand (R_NGAM) is reversed
modelRevATP = modelLeak;
modelRevATP = changeObjective(modelRevATP,'ATPM_c');
modelRevATP = changeRxnBounds(modelRevATP,'ATPM_c',-1000,'l');
modelRevATP = changeRxnBounds(modelRevATP,'ATPM_c',0,'u');
printConstraints(modelRevATP, -1000, +1000);
FBARev = optimizeCbModel(modelRevATP, 'max', 'one');
printFluxVector(modelRevATP, FBARev.x, true);

if FBARev.stat == 0 | all(FBARev.x == 0)
    disp('Good, no matter is produced when ATP demand is reversed');
else
    disp('Error, Matter is produced when ATP demand is reversed');
end

%Adding an artificial reaction validating that NADP production did not occur in the
%absence of CO2 and light.
modelNADP = modelLeak
modelNADP = addReaction(modelNADP, 'R_NADP_DUMMY',...
    'metaboliteList', {'nadph[C_c]', 'h[C_c]', 'nadp[C_c]'},...
    'stoichCoeffList', [-1; 1; 1], 'reversible', false);
modelNADP = changeObjective(modelNADP, {'R_NADP_DUMMY'});
checkObjective(modelNADP) % check biomass function
FBANADP = optimizeCbModel(modelNADP,'max', 'one')
% Print non-zero fluxes if any
printFluxVector(modelNADP, FBANADP.x, true);
rxnID_NADP_DUMMY = findRxnIDs(modelNADP, 'R_NADP_DUMMY');
sol_R_NADP_DUMMY = FBANADP.x(rxnID_NADP_DUMMY);

if FBANADP.stat == 0 | all(FBANADP.x == 0)
    disp('Good, no NADP is produced from nothing');
else
    disp('Error, NADP is produced from nothing');
end

% Testing for NADPH synthesis with no inputs
modelNADPH = modelLeak;
modelNADPH = addReaction(modelNADPH, 'R_NADPH_DUMMY',...
    'metaboliteList', {'nadph[C_c]', 'h[C_c]', 'nadp[C_c]'},...
    'stoichCoeffList', [1; -1; -1], 'reversible', false);
modelNADPH = changeObjective(modelNADPH, {'R_NADPH_DUMMY'});
checkObjective(modelNADPH); % check biomass function
FBANADPH = optimizeCbModel(modelNADPH, 'max', 'one');
% Print non-zero fluxes if any
printFluxVector(modelNADPH, FBANADPH.x, true);
rxnID_NADPH_DUMMY = findRxnIDs(modelNADPH, 'R_NADPH_DUMMY');
sol_R_NADPH_DUMMY = FBANADPH.x(rxnID_NADPH_DUMMY);

if FBANADP.stat == 0 | all(FBANADP.x == 0)
    disp('Good, no NADPH is produced from nothing');
else
    disp('Error, NADPH is produced from nothing');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model verification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General verification of model fields
verifyModel(model1, 'simpleCheck', true)
assert(ans == 1)

% Verify charge blance
% verifyModel(model1, 'chargeBalance', true) % This function cannot be run
% with biomass function components, demand reactions and exchange reactions.

% set exchange reaction boundaries
[selExc, selUpt] = findExcRxns(model1) %, 0, 0)
uptakes = model1.rxns(selUpt) 
exchanges = model1.rxns(selExc)   
EXrxn = exchanges

Bio_Toremove = {'biomass_DNA_c', 'biomass_RNA_c',...
'biomass_carb_c', 'biomass_mem_lipids_c', 'biomass_pigm_h', 'biomass_pro_c', 'bof_c_accumulation_c',...
'sink_Asn-X-Ser_Thr_c', 'DM_biomass_c_acc_c', 'biomass_c_storage_c'}
DM_Remove = {'DM_2mop_m', 'DM_5drib_c', 'DM_dmsp_c',...
'DM_fald_m', 'DM_for_c', 'DM_indole_c', 'DM_m2masn_c', 'DM_mhpglu_c',...
'DM_no3_c', 'DM_phyt_c', 'DM_thmppp_c', 'DM_tre_c'}

modelmass1 = removeRxns(model1, EXrxn)
modelmass1 = removeRxns(modelmass1, Bio_Toremove)
modelmass1 = removeRxns(modelmass1, DM_Remove)
modelmass1 = removeRxns(modelmass1, {'NADHOR_m', 'PSI_u'})
% The reactions 'NADHOR_m' and 'PSI_u' are also charged balance. A coupling
% constant is added to those reactions

ChargeBalanceCheck = verifyModel(modelmass1, 'chargeBalance', true)

Ind = find(ChargeBalanceCheck.chargeBalance.imBalanceCharge ~= 0) 
modelmass1.rxns(Ind)
% The reaction 'CBFC_u' is charge balanced but with a minor rounding error of '2E-04'

% Do a mass balance check
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(modelmass1);
length(imBalancedMass)

Ind = find(~isempty(imBalancedMass)) % Find Rxns that are not mass balanced

if isempty(imBalancedMass)
    disp('Good, all masses are balanced')
else
    disp('Error, mass balance not Ok')
end

% Test deadend metabolites
model1.mets(detectDeadEnds(model1));
fprintf('Number of Deadend metabolites : %4.0f\n', length(model1.mets(detectDeadEnds(model1))))

% test for blocked reactions
blocked_Rxns = identifyFastBlockedRxns(model1, model1.rxns(:));
fprintf('Number of blocked reactions : %4.0f\n', length(blocked_Rxns))

% Identify all blocked metabolites (downstream of a gap)
gapFind(model1);
fprintf('Number of blocked metabolites : %4.0f\n', length(gapFind(model1)))

% Check uniqueness of reactions and metabolites name
checkCobraModelUnique(model1);

% Check duplicate reactions
method = 'S'
[modelout, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model1, 'method')
method = 'FR'
[modelout, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model1, 'method')
if isempty(removedRxnInd)
    disp('Good, no duplicate reactions');
else
    disp('Error, there is at least one duplicated reaction');
end

% Check that there are no empty constraints
if isempty(getEmptyConstraints(model1))
    disp('Good, no empty constraints');
else
    disp('Error, there is at least one empty constraint');
end

% Time the code
toc;