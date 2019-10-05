%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shadow prices, reduced costs as well as reaction/gene deletion analysis
%
% For reaction/gene deletion analysis, the file solveCobraQP.m must be
% modified since the MOMA method does not work otherwise. Appears to have a bug with gurobi...
% I commented the lines 600-601-602 in the function solveCobraQP.m (in
% '/Users/mlavoie/cobratoolbox/src/base/solvers')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% load the model
run Model_setup

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating shadow prices, which are sensitivity analysis of the fluxes involving metabolites on the predicted growth rate by FBA
% Shadow prices (.y) and Reduced cost (.w) vectors are only computed with
% standard (zero-norm) FBA. This is because with the quadratic minimization
% (norm=1) there are one objective followed by a secondary objectives and
% hence there would be too sets of shadow prices or reduced costs for the
% FBA. Both results will not take into account the whole FBA analysis with
% quadratic minimization.

FBAsolution = optimizeCbModel(model1, 'max');
% Shoadow prices are in the field "FBAsolution.y" or "FBAsolution.dual"
% all(FBAsolution.y == FBAsolution.dual)

% Shadow Prices : The increase (if positive) or decrease (if negative) in growth rate for an addition of 1 unit

% Representing all non-zeros shadow prices
% NZSP = find(FBAsolution.y ~= 0)
% length(NZSP)
NZSP = find(abs(FBAsolution.y) > 1E-09); % Accuracy = 1E-09
length(NZSP);
FBAsolutionNZSP = FBAsolution.y(NZSP);
ZSP = find(abs(FBAsolution.y) < 1E-09);
length(ZSP);
FBAsolutionZSP = FBAsolution.y(ZSP);
assert(length(ZSP) + length(NZSP) == length(FBAsolution.y));
[shadowPrices,MetID]= sort(FBAsolutionNZSP, 'descend');
a = model1.mets(MetID); % Metabolite names with NZSP
NZSP_Table = [a, num2cell(shadowPrices)];
% create a table
NZSP_table = cell2table(NZSP_Table,'variableNames', {'Met_names','Shadow_Prices'});
% write to file
writetable(NZSP_table,'NZSP_table.xls');

NZSP_low = find(abs(FBAsolution.y) > 1E-09 & abs(FBAsolution.y) < 0.3);

% Compute the proportion of non-zeros and zeros shadow prices
fprintf('The proportion of metabolites with non-zero shadow prices is : %5.1f\n', 100 * length(NZSP) / (length(FBAsolution.y)));
fprintf('The proportion of metabolites with zero shadow prices is : %5.1f %\n', 100 * length(ZSP) / (length(FBAsolution.y)));
fprintf('Shadow prices varied between %5.2f and %5.2f\n', max(FBAsolutionNZSP), min(FBAsolutionNZSP));
fprintf('The proportion of metabolites with shadow prices between -0.3 and 0.3 is : %5.1f\n', 100 * length(NZSP_low) / (length(FBAsolution.y)));

% Subsystem of metabolites with higher shadow prices

% Select only the metabolites with the 50 highest Shadow Prices
a_high = cell(50,1);
%a_high_c = cellfun(@num2str, cellarray, 'UniformOutput', false)
a_high([1:50]) = {a{1:50}};
shadow_high = shadowPrices(1:50);
NZSP_high_table = [a_high, num2cell(shadow_high)];
% create a table
NZSP_high_table = cell2table(NZSP_high_table,'variableNames', {'Met_Names','Shadow_Prices'});
% write to file
writetable(NZSP_high_table,'NZSP_high_table.xls');

% Find subsystems in which the metabolites with the 50 highest shadow
% prices are
NZSP_subSystem = model1.subSystems(findMetIDs(model1, a_high));
vertcat(NZSP_subSystem{:});

% Find reactions associated with those metabolites
modelS = full(model1.S);

i = 0;
shadow_r = [];
shadow_react = [];
for i = 1:length(a_high)
    shadow_r = find(modelS(findMetIDs(model1, a_high(i)),:) ~= 0);
    shadow_react = cat(2, shadow_react, shadow_r); % Grow a vector of reaction IDs with metabolites
end
shadow_react = unique(shadow_react);
NZ_shadow = find(FBAsolution.x(shadow_react) ~= 0);
fprintf('Number of active reactions involving the metabolites : %4.0f \n', length(NZ_shadow));
fprintf('Proportion of active reactions out of all reactions involving the metabolites : %4.0f \n', 100 * length(NZ_shadow)/length(shadow_react));
% Most reactions involving metabolites with high shadow prices are inactives
FBAsolution2 = optimizeCbModel(model1, 'max');
NZ_shadow = find(FBAsolution2.x(shadow_react) ~= 0);
fprintf('Number of active reactions involving the metabolites : %4.0f \n', length(NZ_shadow));
fprintf('Proportion of active reactions out of all reactions involving the metabolites : %4.0f \n', 100 * length(NZ_shadow)/length(shadow_react));

% Select only the metabolites with the 50 LOWEST Shadow Prices
a_low = cell(50,1);
%a_low_c = cellfun(@num2str, cellarray, 'UniformOutput', false)
a_low([1:50]) = {a{length(a)-49:length(a)}};
shadow_low = shadowPrices(length(a)-49:length(a));
NZSP_low_table = [a_low, num2cell(shadow_low)];
% create a table
NZSP_low_table = cell2table(NZSP_low_table,'variableNames', {'Met_Names','Shadow_Prices'});
% write to file
writetable(NZSP_low_table,'NZSP_low_table.xls');

% Find subsystems in which the metabolites with the 50 highest shadow
% prices are
NZSP_subSystem = model1.subSystems(findMetIDs(model1, a_low));
vertcat(NZSP_subSystem{:});

% Find reactions associated with those metabolites
modelS = full(model1.S);

i = 0;
shadow_r = [];
shadow_react = [];
for i = 1:length(a_low)
    shadow_r = find(modelS(findMetIDs(model1, a_low(i)),:) ~= 0);
    shadow_react = cat(2, shadow_react, shadow_r); % Grow a vector of reaction IDs with metabolites
end
shadow_react = unique(shadow_react);
NZ_shadow = find(FBAsolution.x(shadow_react) ~= 0);
fprintf('Number of active reactions involving the metabolites : %4.0f \n', length(NZ_shadow));
fprintf('Proportion of active reactions out of all reactions involving the metabolites : %4.0f \n', 100 * length(NZ_shadow)/length(shadow_react));
% Most reactions involving metabolites with low shadow prices are inactives
FBAsolution2 = optimizeCbModel(model1, 'max');
NZ_shadow = find(FBAsolution2.x(shadow_react) ~= 0);
fprintf('Number of active reactions involving the metabolites : %4.0f \n', length(NZ_shadow));
fprintf('Proportion of active reactions out of all reactions involving the metabolites : %4.0f \n', 100 * length(NZ_shadow)/length(shadow_react));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced cost analysis (sensitivity analysis of reactions)
% Note that by definition the reduced cost are positive. Meaning that a reaction A with a reduced cost of 35 would lead to an increase in the objective value by 35 if the flux through this reaction would be increased by 1 flux unit.
% This analysis works only with the standard FBA as above.
FBAsolution = optimizeCbModel(model1, 'max');
% all(FBAsolution.w == FBAsolution.rcost)

%%%% Representing all non-zeros reduced costs
% NZRC = find(FBAsolution.w ~= 0)
% length(NZRC)
NZRC = find(abs(FBAsolution.w) > 1E-09);
length(NZRC);
FBAsolutionNZRC = FBAsolution.w(NZRC);
ZRC = find(abs(FBAsolution.w) < 1E-09);
FBAsolutionZRC = FBAsolution.w(ZRC);
assert(length(NZRC) + length(ZRC) == length(FBAsolution.w));

[redCost,rxnID] = sort(FBAsolution.w(NZRC), 'descend');
a = model1.rxns(rxnID(1:length(NZRC)));
b = redCost(1:length(NZRC));
c = model1.subSystems([rxnID(1:length(NZRC))], 1);
c = vertcat(c{:});
NZRC_all_table = [a, num2cell(b), c];
varNames = {'Reaction','Reduced_Costs','Subsystems'}; %Rxn = a ; Reduced_Cost = num2cell(b) ; Subsystems = c
NZRC_all_table = cell2table(NZRC_all_table, 'VariableNames',varNames);
% write to file
writetable(NZRC_all_table,'NZRC_all_table.xls');

% Compute the proportion of non-zeros and zeros reduced costs
fprintf('The proportion of reactions with non-zero reduced costs is : %5.1f\n', 100 * length(NZRC) / (length(FBAsolution.w)));
fprintf('The proportion of reactions with zero reduced costs is : %5.1f\n', 100 * length(ZRC) / (length(FBAsolution.w)));
fprintf('Reduced costs vary between %5.2f and %5.2f\n', max(FBAsolutionNZRC), min(FBAsolutionNZRC));

% reduced costs among active reactions
Ind_FBA_Act = find(abs(FBAsolution.x) > 1E-09);
length(Ind_FBA_Act);
FBAsolution.w(Ind_FBA_Act);

% Most reactions with non-zero reduced costs are inactives
NZ_reduced = find(FBAsolution.x(NZRC) ~= 0);
model1.rxns(NZ_reduced); % Reactions for which reduced costs are different from zero
fprintf('Number of reactions with non-zero reduced costs that are actives: %4.0f \n', length(NZ_reduced));
fprintf('Proportion of reactions with non-zero reduced costs that are actives : %4.0f \n', 100 * length(NZ_reduced)/length(NZRC));

FBAsolution2 = optimizeCbModel(model1, 'max', 'one');
NZ_reduced = find(FBAsolution2.x(NZRC) ~= 0);
fprintf('Number of reactions with non-zero reduced costs that are actives: %4.0f \n', length(NZ_reduced));
fprintf('Proportion of reactions with non-zero reduced costs that are actives : %4.0f \n', 100 * length(NZ_reduced)/length(NZRC));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use single-reaction deletion with FBA: Compute FBA when removing one
% reaction at a time.
% This gives an idea of the global robustness of the network to
% gene/reaction deletion 

% Add constraints on the uptake of N and other elements. Otherwise, N exchange is frequently
% at -1000 mmol / gDW / h when a reaction is deleted, which is unrealistic.
% check urea and glycolate production...
% A flux of -1000 mmol NO3- / h-1 / gDW = 15 000 fmol / cell / h (* 2 g DW/gC
% and 7.5 x 10-12 gC/cell), this is much higher than what is expected for
% diatoms (Lomas and Gilbert, 2000) and even for N-limited diatoms (Zevenboom and Mur, 1981)

FBAsolution = optimizeCbModel(model1, 'max', 'one');
printFluxVector(model1, FBAsolution.x, false, true);

%{
% Compute the maximum diffusive uptake rate of nutrients for F. cylindrus
% Diffusion coefficients are from (Li and Gregory, 1974) at 25 °C (in cm^2/
% sec)
DNO3 = 19*1E-06      % NO3-
DH = 93.1E-06        % H+
DHCO3 = 11.8E-06     % HCO3- or 1.18E-09 m^2/sec
% DO = ??             % O2
DP = 6.12E-06        % PO43- diffusion coefficient
DSO4 = 10.7E-06      % SO42-
DMg = 7.05E-06       % Mg2+

% Kariuki and Dewald, 1995
Ccell = 7.5E-12
DWcell = Ccell * 2
Radius = 5 * 1E-04             % Mean mesured radius in cm ?? (Reference )?
Bound = 2 * Radius             % Boundary layer thickness
Cbulk = 10 * 1000              % pmol / cm3 or nM
Csurf = 0 
Area = 4 * pi() * Radius^2     % Cell surface (cm^2 / cell)

Jmax =  ( (4 * pi() * DNO3) * ( (Radius * (Radius+Bound)) / ((Radius+Bound) - Radius) ) * (Cbulk - Csurf) ) / Area    % Jmax in pmol / cm^2 / sec
Jmax = Jmax * 60 * 60 * Area * (DWcell)^-1 / (1E+09) % * 60 sec/min * 60 min/h * cm^2/cell * (gDW)-1 / 1000 mmol/umol / 1000 uml/nmol / 1000 nmol/ pmol : Jmax in mmol (gDW)-1 h-1

[selExc, selUpt] = findExcRxns(model1); %, 0, 0)
uptakes = model1.rxns(selUpt) ;
exchanges = model1.rxns(selExc) ;  
EXrxn = exchanges;
%}


%{
% Try to force the model around a given C:N value (7.5 or 5.5) and a
measured N:P value of 12.5 (Garcia et al 2018).
% This suggests that another P species is stored inside the cell since a DM
reaction (DM_phyt_c) is activated and most P entering the cell is stored
there.

N_P = 13
fold_change = 1.1
N_low =  - Ps / (C_N * fold_change) ; % in mmol / gDW / h % 6.6
N_high = - Ps / (C_N / fold_change);
P_low = - Ps / (C_N * N_P * fold_change);  % in mmol / gDW / h % 106
P_high = - Ps / (C_N*N_P / fold_change); % not 2
Pho_low = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) * fold_change ; % not 5
Pho_high = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) / fold_change ;
%}

fold_change = 1.5;
%N_low =  - Ps / (C_N * fold_change) ; % in mmol / gDW / h % 6.6
%N_high = - Ps / (C_N / fold_change);
%P_low = - Ps / (C_N * N_P * fold_change);  % in mmol / gDW / h % 106
%P_high = - Ps / (C_N*N_P / fold_change); % not 2
Pho_low = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) * fold_change ; % not 5
Pho_high = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) / fold_change ;

%model1 = changeRxnBounds(model1,'EX_no3_e', N_high, 'l'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_no3_e', N_low, 'u'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_pi_e', P_high, 'l'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_pi_e', P_low, 'u'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_photon_e', Pho_low, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_photon_e', P_low, 'u'); % Bd: -1000 / 0

FBAsolution = optimizeCbModel(model1, 'max', 'one');
% printFluxVector(model1, FBAsolution.x, false, true);

%{
% printFluxVector(model1, fluxSolution(:,2), false, true);
fluxSolution(findRxnIDs(model1, 'EX_photon_e'), IDMore)

% Adding constraints on N and P uptake does not change much the conclusion
% : Still, around 30% of reaction deletion does not affect the growth rate

% Adding a constraint on photon uptake (5-fold higher or lower than the
% control), does not change much the conclusion, again 30% of reaction
% deletion does not affect the growth rate. However, several robust
% reactions absorbs more photon , so need more energy to grow optimally

model1.rxns(Ind_Comp) % last one, 187e, no CEF, increase photon uptake
surfNet(model1, model1.rxns(Ind_Comp))
model1.rxns(IDMore)


%{
EX_photon_e         	       -10.5
EX_hco3_e           	     -0.7778
EX_no3_e            	    -0.06346
EX_h_e              	     -0.8422
EX_h2o_e            	      0.1738
EX_o2_e             	       1.011
EX_pi_e             	   -0.005693
EX_so4_e            	   -0.002729
EX_mg2_e            	  -0.0001841
DM_biomass_c_acc_c  	     0.01563
%}
%}

% Reaction deletion analysis with the function singleRxnDeletion() computes
% by default a quadratic minimization FBA (See the source code)

% Try the 'FBA' method first, then the 'MOMA' method
% The MOMA method does not work well. Appears to have a bug with gurobi...
% I commented the lines 600-601-602 in the function solveCobraQP.m (in
% '/Users/mlavoie/cobratoolbox/src/base/solvers')

% If needed, adjust the tolerance to 1e-06 rather than 1e-09. Maximum tolerance is 0.01
% global CBT_QP_PARAMS
% CBT_QP_PARAMS.feasTol=1e-09

methodDel = {'FBA', 'MOMA'};
k = 0;
for k = 1:length(methodDel)
if k == 1
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('First, compute the effect of single reaction deletion using the FBA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleRxnDeletion(model1, methodDel{k});
else
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('Second, compute the effect of single reaction deletion using the MOMA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleRxnDeletion(model1, methodDel{k});
end

%%%%%%%%%%%%%
% Calculating the effect of reaction deletion on growth rate
FBAsolution = optimizeCbModel(model1, 'max', 'one');
IndActive = find(FBAsolution.x ~= 0); % 567 indices
%IndActive = find(abs(FBAsolution.x) > 1E-09); % same result
IndInac = find(FBAsolution.x == 0); % Useful for the last part looking at inactive fluxes
assert(length(IndActive) + length(IndInac) == length(FBAsolution.x));
grRatioAct = grRatio(IndActive);
fluxSolutionAct = fluxSolution(:, IndActive);
size(fluxSolutionAct);

NaN_ratio = find(isnan(grRatioAct));
grRatioAct(NaN_ratio) = 0; % Assume that NaN values in grRatio or 'Infeasible' FBA yields a growth rate = 0

length(grRatioAct);
length(grRatioAct >= 0);
assert(length(grRatioAct) == length(grRatioAct >= 0)); % Make sure no negatiev values are there
assert(length(grRatioAct) == size(fluxSolutionAct,2));

model1.rxns(IndActive(NaN_ratio)) ; % find the reaction names of reaction deletion yielding infeasible flux distribution

% Summary of grRatio for each reaction names
horzcat(model1.rxns(IndActive), num2cell(grRatioAct));
fprintf('Number of active reactions: %4.2f \n', length(IndActive));
fprintf('Proportion of active reactions: %4.2f \n', 100 * length(IndActive) / length(FBAsolution.x));
fprintf('Number of inactive reactions: %4.2f \n', length(IndInac));
fprintf('Number of reaction with infeasible reaction deletions: %4.2f \n', length(NaN_ratio));

i = 0;
z = 0;
L1 = 0;
Lactive = 0;
L = 0;
NumberRobustRxn = 0;
NumberActiveRxn = 0;
Percentage_RobustRxn = 0;
for i = 1:1000
j = i / 1000;
IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
Lactive(i) = length(IndActive);
%    for z = 1:length(IndFlux)
%        Flux_Vector = fluxSolution(:,IndFlux(z)); % Extract a Flux vector of 2156 values for each relevant value of grRatio (z)
%        Ind_Act = find(Flux_Vector ~= 0);
%        L(z) = length(Ind_Act);
%    end
%Lactive(i) = sum(L) / length(L)
NumberRobustRxn(i) = L1(i); % Number of robust reactions among all active reactions
NumberActiveRxn(i) = Lactive(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberActiveRxn(i);
end
figure (1);
Y = Percentage_RobustRxn;
length(Y);
X = 1/1000:1/1000:1;
X = [0,X];
Y = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + Y(1),Y]  ;  
Y = fliplr(Y);
Y(end) = 100 ;  % This correct for the 3.1 % of reactions for which deletion abolishes the growth rate, but yield negative grRatioAct values.
X = X * 100 ;% Makes X-axis in percentage
length(X);
plot(X, Y);
xlim([-2 102]);
ylim([0 103]);
xlabel('Maximum growth rate inibition (%)', 'FontSize', 12); % 1 - (KO / WT) x 100
ylabel('Cumulative percentage of reactions out of all active reactions', 'FontSize', 12);
ax = gca;
ax.FontSize = 10; 
% title('Robustness of the FBA result to reaction deletion')
% This plot gives the the cumulative proportion of reactions out of all
% active reactions, which when deleted one by one yield a given growth rate
% inhibition or less.
if k == 1
saveas(gcf,'RxnDeletion_FBA.pdf');
else
saveas(gcf,'RxnDeletion_MOMA.pdf');
end

% How the proportion of active reactions change for different reaction deletions?
%{
i = 0
L = 0
for i = 1:size(fluxSolutionAct, 2)
    Flux_Vector = fluxSolutionAct(:,i); 
    Ind_Act = find(Flux_Vector ~= 0);
    L(i) = length(Ind_Act);
end
Mean_L = sum(L) / length(L)
Std = std(L)
CV = ( Std / Mean_L ) * 100
Median = median(L)
IQR = iqr(L)
% This means that the number of active reactions can become smaller or
higher than in the Wild-type (595 reactions)
%}
%%%%%%%%%%%%%

% Extracting fluxes of most robust reactions (Flux_mat: fluxes in row and rxn ID (ID_Rob) in column)
ID_Rob = find(grRatioAct >= 0.99 & grRatioAct <= 1.01);
% ID_other = find(grRatioAct < 0.99)
fprintf('Number of more robust reactions: %4.2f \n', length(ID_Rob))
fprintf('Proportion of more robust reactions relative to active reactions: %4.2f \n', 100 * length(ID_Rob) / length(IndActive))

L_r = length(ID_Rob);

Flux_mat = zeros(size(fluxSolutionAct,1),L_r);
z = 0;
for z = 1:L_r
        Flux_mat(:,z) = fluxSolutionAct(:,ID_Rob(z));
end

% Listing all reactions of most robust reactions
rxn_list_red = model1.rxns(IndActive); % Reduced rxn list of active reactions
modelRob = rxn_list_red(ID_Rob);

% Looking for the proportion of reactions requiring ATP or NADPH that are
% MORE robusts
IDMore = findRxnIDs(model1,modelRob); % Reaction IDs of more robust reactions
modelS = full(model1.S);

IDatpc = find(modelS(findMetIDs(model1, 'atp_c'),:) ~= 0);
IDatph = find(modelS(findMetIDs(model1, 'atp_h'),:) ~= 0);
IDatpm = find(modelS(findMetIDs(model1, 'atp_m'),:) ~= 0);
IDatpx = find(modelS(findMetIDs(model1, 'atp_x'),:) ~= 0); % No 'atp_u'

RxnID_atpc = find(ismember(IDMore, IDatpc') ~= 0); %ismember : % find if IDatpc is in IDLess
RxnID_atph = find(ismember(IDMore, IDatph') ~= 0);
RxnID_atpm = find(ismember(IDMore, IDatpm') ~= 0);
RxnID_atpx = find(ismember(IDMore, IDatpx') ~= 0);

fprintf('Proportion of more robust reactions using atp : %4.2f \n', 100 *  ( length(RxnID_atpc) + length(RxnID_atph) + length(RxnID_atpm) + length(RxnID_atpx) ) / length(IDMore))

% How the fluxes changes through all biomass reactions for robust reaction
% deletions?
IDbio = find(contains(model1.rxns(), 'biomass'));
i = 0;
Mean_bio = [];
Range_bio = [];
for i = 1:length(IDbio)
    Mean_bio(i) = mean( fluxSolution(IDbio(i), IDMore) );
    Range_bio(i) = max( fluxSolution(IDbio(i), IDMore) ) - min ( fluxSolution(IDbio(i), IDMore) );
    
end
Mean_bio; % vector of mean  of fluxes through the 8th biomass reactions for all robust reactions
Range_bio;% vector of ranges of fluxes through the 1-8th biomass reactions
Range_bio ./ Mean_bio;
fprintf('Error in percent on each main biomass component: %4.3f \n', 100 * Range_bio ./ Mean_bio)
% Since range_bio is much smaller than Mean_bio, 
% this loop shows that the fluxes do not vary significantly for all 8
% biomass equation when a given robust reaction is removed
% Therefore, in all cases, not only the growth rte is robust to reaction deletion, but
% also, the rate of synthesis of all major cell components is robust.

% How the proportion of active reactions change for different reaction deletions for ROBUST reactions?

i = 0;
L = 0;
Compensatory = 0;
for i = 1:size(Flux_mat, 2)
    Flux_Vector = Flux_mat(:,i); 
    Ind_Act = find(Flux_Vector ~= 0);
    %Ind_Act = find(abs(Flux_Vector) > 1E-09);
    L(i) = length(Ind_Act); % Vector of number of active reactions for each rxn deletion
    if L(i) > length(IndActive) 
        Compensatory(i) = i; % vector reporting the rxn deletion leading to a compensatory response
    end
end
% first reaction : 606-576 = 30 compensatory active reactions
% second reaction : 595-576 = 19 compensatory active reactions

Mean_L = sum(L) / length(L);
Std = std(L);
CV = ( Std / Mean_L ) * 100;
Median = median(L);
IQR = iqr(L);
% Proportion of reaction deletion with number of fluxes higher than 567
% (control WT has 567 active reactons)
length(find(L > length(IndActive) ));
length(find(L < length(IndActive) ));
fprintf('Number of reaction deletion (KO strain) with number of non-zero fluxes > number of fluxes in control : %4.2f \n', length(find(L > length(IndActive))) )
fprintf('Number of reaction deletion (KO strain) with number of non-zero fluxes < number of fluxes in control : %4.2f \n', length(find(L < length(IndActive))) )

% What are the 100% of reaction deletion that can activate more reactions than in the control?
length(Compensatory);
Compensatory(Compensatory == 0) = []; % Remove zeros from Compensatory, which represents L < 565
length( Compensatory );
modelRob(Compensatory); % Here this is equivalent ot modelRob
% So, deleting this set of reaction leads to a compensatory response
% : in all the cases, activation of more reactions than WT is observed. 
% This suggests that the resilience of the model to perturbation is due to
% compensatory activation of other reactions.

% Find reaction ID that are activated in KO (i.e., row in flux mat with non-zero
% fluxes, but that are not in IndActive or in control. Ncomp
% Find also reactions that are deactivated in KO relative to WT. Ncomp2
i = 0;
Ind_Act = 0;
Flux_Vector = 0;
Nact = [];
Ncomp = [];
Ncomp2 = [];
rangeComp = [];
maxRange = [];
j = 0;
for i = 1:size(Flux_mat, 2)
    Flux_Vector = Flux_mat(:,i); 
    length(Flux_Vector);
    Ind_Act = find(Flux_Vector ~= 0);
    Nact(i) = length(Ind_Act);
    %Ind_Act = find(abs(Flux_Vector) > 1E-09);
    %length(Ind_Act)
    Ind_Comp = setdiff(Ind_Act, IndActive); % find newly active reaction in KO relative to WT, i.e., findindices in Ind_Act (KO active reaction IDs) that are not in IndActive (control active reaction IDs)
    rangeComp = [];
    j = 0;
    for j = 1:length(Ind_Comp)
        rangeComp(j) = range(Flux_mat(Ind_Comp(j),:));
    end
    maxRange(i) = max(rangeComp);
    Ncomp(i) = length(Ind_Comp);
    Ind_Comp2 = setdiff(IndActive, Ind_Act); % find reactions in WT that was active , but are not active anymore
    Ncomp2(i) = length(Ind_Comp2);
    assert( length(IndActive) - length(Ind_Comp2) + length(Ind_Comp) == length(Ind_Act)) % Total active reaction in control - reactions becoming non-active + new active reaction = total active reactions in KO
    if i == 1
        Ind_Comp_Store = Ind_Comp;
    end
end

% Calculating range in fluxes for compensatory reactions
maxRange; % vector of maximum differences among fluxes of a row

% not good, indices are only for the last iteration : Ind_New_Comp = Ind_Comp % Store Indices of compensatory reactions in 'Ind_New_Comp' for subsequent use.
% plot(Nact, Ncomp, 'o')
mean(Ncomp);
std(Ncomp);
CV =  std(Ncomp) / mean(Ncomp) * 100;
median(Ncomp);
IQR = iqr(Ncomp);
range(Ncomp);
min(Ncomp);
max(Ncomp);
% plot(Nact, Ncomp2, 'o')
mean(Ncomp2);
std(Ncomp2);
CV = std(Ncomp2) / mean(Ncomp2) * 100;
median(Ncomp2);
IQR = iqr(Ncomp2);
min(Ncomp2);
max(Ncomp2);
% In general (based on mean and median), slightly more activation than
% inacivation of reactions occured in KO compared to WT.
% Therefore, net number of activated reactions in KO is higher than total
% number of activated reactions in WT.

%{
% Plot all fluxes for KO and WT for each reaction deletion
% Deleting Rxn 1 : ATPS_m
% Ind_Comp_Store
% model1.rxns(Ind_Comp_Store)  
% surfNet( model1, model1.rxns(Ind_Comp_Store))
% Flux_mat(Ind_Comp_Store(1), :)
% Rxn #1551  CYOR_m, Bd: 0 / 1000, Ubiquinol-cytochrome c oxidoreductase Complex III
% Rxn #1552  CYOO_m, Bd: 0 / 1000, Cytochrome c oxidase Complex IV
% Rxn #1897  ETFQO_m, Bd: 0 / 1000, Electron transfer flavoprotein-ubiquinone oxidoreductas
% Rxn #2103  CEF_h, produces ATP

% In the control, 2pg_m comes from ENO_m (Rxn 89 : h2o_m + pep_m <=> 2pg_m ) 
% in the knock-out, 2 pg_m, which is ultimatly used for PPP, is produced
% from pyr_m and ser__L_m

% Rxn #2138  SPT_m, Bd: -1000 / 1000, SPT_m
fluxSolution(2138,1) % + 0.068 production of hpyr_m and ala__L_m from pyr_m and ser__L_m
% Rxn #2139  HYPRRx_m, Bd: -1000 / 1000, HYPRRx_m
fluxSolution(2139,1) % + 0.068 production of glyc__R_m and nad_m from nadh and hpyr_m
% Rxn #2140  GLYCK2_m, Bd: 0 / 1000, GLYCK2_m
fluxSolution(2140,1) % + 0.068 production of adp_m and 2pg_m from glyc__R_m and atp

% Here production of 3pg_m comes from transport of 2pg_m as in the control
% #277  PGAM_m, 3pg_m <=> 2pg_m 
fluxSolution(277,1)  % - 0.068 production of 3pg_m from 2pg_m
%  #89  ENO_m, h2o_m + pep_m <=> 2pg_m 
fluxSolution(89,1) % zero flux
% #384  PGK_m, atp_m + 3pg_m <=> adp_m + 13dpg_m  
fluxSolution(384,1) % + 0.068 production of 13dpg_m (glycolysis in reverse)
% #413  SERA_m, nad_m + 3pg_m <=> h_m + nadh_m + 3php_m  
fluxSolution(413,1) % zero
% #87  GAPDH_m, h_m + nadh_m + 13dpg_m <=> pi_m + nad_m + g3p_m  
fluxSolution(87,1) % + 0.068
% #383  TPI_m, g3p_m <=> dhap_m  
fluxSolution(383,1) % zero
%  #388  TPTPt_m, pi_m + g3p_c <=> pi_c + g3p_m  
fluxSolution(388,1) % -0.068 TRANSFER of grp from mito to cytosol
%  #394  EDA_m, 2ddg6p_m <=> pyr_m + g3p_m  
fluxSolution(394, 1) % zero
% #111  PYDXS_c, gln__L_c + g3p_c + ru5p__D_c -> h_c + 3 h2o_c + pi_c + glu__L_c + pydx5p_c  
fluxSolution(111, 1) % zero
% #272  TPTPti_h, pi_h + g3p_c -> pi_c + g3p_h 
fluxSolution(272, 1) % zero
% #313  TPI_c, g3p_c <=> dhap_c  
fluxSolution(313, 1) % 0.036 production of dhap_c
% #328  TKL2_c, f6p_c + g3p_c <=> e4p_c + xu5p__D_c  

% What fluxes produce pyr_m in KO  and control ?
FBAsolution.x([17, 23, 52, 59, 95, 104, 390, 394, 1439, 1530, 1867])
fluxSolution([17, 23, 52, 59, 95, 104, 390, 394, 1439, 1530, 1867],1)

% In the KO, Rxn #1867  ALATA_L_m, produces pyr_m, but not in the control

% Consuming pyr_m
FBAsolution.x([42, 71, 1469, 1471, 1476, 1744, 2138])
fluxSolution([42, 71, 1469, 1471, 1476, 1744, 2138],1)

model1.rxns([17,390,1867])
surfNet(model1, model1.rxns([17,390,1867]))
% In WT, pyr_m is consumed though PDH and PYC, while in KO, pyr_m is
% redirected to upward glycolysis with SPT_m

% draw_by_met(model1, {'2pg_m'}, 'true', 1, 'struc', {''}, fluxSolution(:,1))

% alanine transaminase : ALATA_L_c (in mitochondria and cytosol)
% aspartate transaminase (ASPTA_L_c) (in mitochondria and cytosol)
% several transport reactions
%  #2140  GLYCK2_m and #2139  HYPRRx_m and #2138  SPT_m: 
% pyr_m -> hpyr_m -> 2pg_m
% Energy production : CEF_h, CYOR_m
% P and N balance : DM_NO3, PYROHt_c : h2o_c + ppi_c -> h_e + 2 pi_c  
%}
    
% Large change in fluxes for those exchange reactions when deleting robust
% reactions
KO_photon = fluxSolution(findRxnIDs(model1, 'EX_photon_e'), IDMore);
WT_photon = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e'));
PropPhoton = 100 * length( find( abs(KO_photon) ./ abs( abs(WT_photon) ) > 1 ) ) / length(IDMore);
fprintf('Proportion of robust reaction deletion with up-regulated photon uptake : %4.2f \n', PropPhoton);

fluxSolution(findRxnIDs(model1, 'EX_no3_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_no3_e'));
fluxSolution(findRxnIDs(model1, 'DM_no3_c'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'DM_no3_c'));

% Change a little bit :
fluxSolution(findRxnIDs(model1, 'EX_o2_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_o2_e'));
fluxSolution(findRxnIDs(model1, 'EX_h2o_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_h2o_e'));
fluxSolution(findRxnIDs(model1, 'EX_h_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_h_e'));

% Does NOT CHANGE AT ALL :
fluxSolution(findRxnIDs(model1, 'EX_mg2_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_mg2_e'));
fluxSolution(findRxnIDs(model1, 'EX_so4_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_so4_e'));
fluxSolution(findRxnIDs(model1, 'EX_pi_e'), IDMore);
FBAsolution.x(findRxnIDs(model1, 'EX_pi_e'));

% Key reactions, large changes :
KO_CEF = fluxSolution(findRxnIDs(model1, 'CEF_h'), IDMore); % CEF_h is often activated for robust reactions
WT_CEF = FBAsolution.x(findRxnIDs(model1, 'CEF_h'));
PropCEF = 100 * length( find( KO_CEF > 0) )  / length(IDMore);
fprintf('Proportion of robust reaction deletion with up-regulated cyclic electron transport : %4.2f \n', PropCEF);

fluxSolution(findRxnIDs(model1, 'ATPS_h'), IDMore) ;
FBAsolution.x(findRxnIDs(model1, 'ATPS_h'));
KO_ATP_h = fluxSolution(findRxnIDs(model1, 'ATPS_h'), IDMore);
WT_ATP_h = FBAsolution.x(findRxnIDs(model1, 'ATPS_h'));
PropATPh = 100 * length( find( KO_ATP_h ./ WT_ATP_h > 1) )  / length(IDMore);
fprintf('Proportion of robust reaction deletion with up-regulated ATPS_m : %4.2f \n', PropATPh);

fluxSolution(findRxnIDs(model1, 'ATPS_m'), IDMore) ;
FBAsolution.x(findRxnIDs(model1, 'ATPS_m'));
KO_ATP_m = fluxSolution(findRxnIDs(model1, 'ATPS_m'), IDMore);
WT_ATP_m = FBAsolution.x(findRxnIDs(model1, 'ATPS_m'));
PropATPm = 100 * length( find( KO_ATP_m ./ WT_ATP_m > 1) )  / length(IDMore);
fprintf('Proportion of robust reaction deletion with up-regulated ATPS_m : %4.2f \n', PropATPm);

% Energy dissipation
fluxSolution(findRxnIDs(model1, 'EX_glyclt_e'), IDMore); % 0
fluxSolution(findRxnIDs(model1, 'PTOX_h'), IDMore); % 0
fluxSolution(findRxnIDs(model1, 'MEHLER_h'), IDMore); % 0
fluxSolution(findRxnIDs(model1, 'RUBISO_h'), IDMore); % photorespiration = 0
fluxSolution(findRxnIDs(model1, 'PGLYCP_h'), IDMore); % photorespiration = 0
fluxSolution(findRxnIDs(model1, 'EX_urea_e'), IDMore); 

%{
% Other tests
fluxSolution(findRxnIDs(model1, 'HYPRRx_m'), IDMore) ;
FBAsolution.x(findRxnIDs(model1, 'HYPRRx_m'));
fluxSolution(findRxnIDs(model1, 'SPT_m'), IDMore) ;
FBAsolution.x(findRxnIDs(model1, 'SPT_m'));
fluxSolution(findRxnIDs(model1, 'GLYCK2_m'), IDMore) ;
FBAsolution.x(findRxnIDs(model1, 'GLYCK2_m' ));
% Enzymes involved in redox state homeostasis
% Rxn #2047  FRNDPR_c, Bd: -1000 / 1000, Ferredoxin:NADP reductase
fluxSolution(findRxnIDs(model1, 'FRNDPR_c'), IDMore)  % reduction of fdxox_c to fdxrd_c using nadph_c : Needed for NOR_c (nitrite reductase cytosolic, using fdxrd_c)
FBAsolution.x(findRxnIDs(model1, 'FRNDPR_c' )) 
% Rxn #1756  THD1_h, 
fluxSolution(findRxnIDs(model1, 'THD1_h'), IDMore)  % negative or positive fluxes...
FBAsolution.x(findRxnIDs(model1, 'THD1_h' )) 
% Rxn #1939  NTRC_h, Bd: 0 / 1000, NADPH-dependent thioredoxin reductase 3
fluxSolution(findRxnIDs(model1, 'NTRC_h'), IDMore)  % Sometimes equal to 0.0021, Synthesis of trdrd_h from trdox_h and nadph_h
FBAsolution.x(findRxnIDs(model1, 'NTRC_h' )) 
% Rxn #2020  ACALDt_m, acald_c <=> acald_m , acetaldehyde transport
fluxSolution(findRxnIDs(model1, 'ACALDt_m'), IDMore)  % Usually at -1000 mmol / gDW / h
FBAsolution.x(findRxnIDs(model1, 'ACALDt_m' ))  
% Rxn #74  THRA_m, (-1000) threonine aldolase, gly_m + acald_m <=> thr__L_m 
fluxSolution(findRxnIDs(model1, 'THRA_m'), IDMore)  % Usually at -1000
FBAsolution.x(findRxnIDs(model1, 'THRA_m' ))  % 0
%%%%% #2020  ACALDt_m (-1000), acald_c <=> acald_m 
fluxSolution(findRxnIDs(model1, 'ACALDt_m'), IDMore)  % Usually at -1000
FBAsolution.x(findRxnIDs(model1, 'ACALDt_m' ))  % 0
% #74  THRA_c (1000), gly_c + acald_c <=> thr__L_c  
fluxSolution(findRxnIDs(model1, 'THRA_c'), IDMore)  % Usually at 1000
FBAsolution.x(findRxnIDs(model1, 'THRA_c' ))  % 0
% #264  THRTL_c (0.0024), Bd: 0 / 1000, Threonine-tRNA ligase
% atp_c + thr__L_c + trnathr_c -> ppi_c + amp_c + thrtrna_c  
fluxSolution(findRxnIDs(model1, 'THRTL_c'), IDMore)  % 0.0024 NO CHANGE !
FBAsolution.x(findRxnIDs(model1, 'THRTL_c' ))  % 0.0024 NO CHANGE !
%   #352  THRD_c (0.01406), Bd: 0 / 1000, L-threonine deaminase
% thr__L_c -> nh4_c + 2obut_c 
fluxSolution(findRxnIDs(model1, 'THRD_c'), IDMore)  % USUALLY IMCREASE !
FBAsolution.x(findRxnIDs(model1, 'THRD_c' ))  % 0.002
%  #1693  THRNA1t_m (-1000), Bd: -1000 / 1000, Threonine:Na+ symporter, mitochondrial
% na1_m + thr__L_m <=> na1_c + thr__L_c  
fluxSolution(findRxnIDs(model1, 'THRNA1t_m'), IDMore)  % USUALLY -1000
FBAsolution.x(findRxnIDs(model1, 'THRNA1t_m' ))  % 0
% #420  THRS_c (0.01647), Bd: 0 / 1000, threonine synthase
% h2o_c + phom_c -> pi_c + thr__L_c  
fluxSolution(findRxnIDs(model1, 'THRS_c'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'THRS_c' ))  % 0
% #2026  PHSERt_h (0.01647), phom_h -> phom_c  
fluxSolution(findRxnIDs(model1, 'PHSERt_h'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'PHSERt_h' ))  % 0
% #419  HSK_h (0.01647), atp_h + hom__L_h -> h_h + adp_h + phom_h  
fluxSolution(findRxnIDs(model1, 'HSK_h'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'HSK_h' ))  % 0.0044
% #418  HSDH_h (0.01647), h_h + nadh_h + aspsa_h <=> nad_h + hom__L_h  
fluxSolution(findRxnIDs(model1, 'HSDH_h'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'HSDH_h' ))  % 0.0044
% #185  ASADH_h (0.01908), h_h + nadph_h + 4pasp_h -> pi_h + nadp_h + aspsa_h  
fluxSolution(findRxnIDs(model1, 'ASADH_h'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'ASADH_h' ))  % 0.007
% #186  ASPK_h (0.01908), atp_h + asp__L_h -> adp_h + 4pasp_h  
fluxSolution(findRxnIDs(model1, 'ASPK_h'), IDMore)  % USUALLY INCREASE
FBAsolution.x(findRxnIDs(model1, 'ASPK_h' ))  % 0.007
% #12  ASPTA_L_h (-0.02706), akg_h + asp__L_h <=> glu__L_h + oaa_h  
fluxSolution(findRxnIDs(model1, 'ASPTA_L_h'), IDMore)  % USUALLY decrease
FBAsolution.x(findRxnIDs(model1, 'ASPTA_L_h' ))  % -0.0507

% gly_m

% #228  PPRGL_c (0.00049), Bd: 0 / 1000, Phosphoribosylamine-glycine ligase
% atp_c + gly_c + pram_c -> h_c + pi_c + adp_c + gar_c  
fluxSolution(findRxnIDs(model1, 'PPRGL_c'), IDMore)  % DOes not change
FBAsolution.x(findRxnIDs(model1, 'PPRGL_c' ))  % 4.8E-04

% #311  GHMT_c (-0.00512), Bd: -1000 / 1000, h2o_c + gly_c + mlthf_c <=> thf_c + ser__L_c 
fluxSolution(findRxnIDs(model1, 'GHMT_c'), IDMore)  % INCREASE (more negative)
FBAsolution.x(findRxnIDs(model1, 'GHMT_c' ))  % -0.0057



%  #320  CBMK_m (0.00239), Bd: 0 / 1000, Carbamate kinase
% atp_m + nh4_m + co2_m -> 2 h_m + adp_m + cbp_m 
fluxSolution(findRxnIDs(model1, 'CBMK_m'), IDMore)  % Some activated.
FBAsolution.x(findRxnIDs(model1, 'CBMK_m' ))  % 0
% cbp_m is converted in citrulline, which is converted to fumarate and
% malate..., oaa_m, oaa_m is converted to pep with PEPCK_m
% oaa_m can also be converted to aspartate and other amino acids

 
fluxSolution(findRxnIDs(model1, 'PEPCK_m'), IDMore)  % Some activated.
FBAsolution.x(findRxnIDs(model1, 'PEPCK_m' ))  % 0

fluxSolution(findRxnIDs(model1, 'CA_h'), IDMore)  % Variable.
FBAsolution.x(findRxnIDs(model1, 'CA_h' ))  % 0.59
 

% #1532  CA_m (0.13398), h_m + hco3_m <=> h2o_m + co2_m  
fluxSolution(findRxnIDs(model1, 'CA_m'), IDMore)  % Several increase (becomes positive)
FBAsolution.x(findRxnIDs(model1, 'CA_m' ))  % -0.04

% #1535  CO2t_m (-0.16317),co2_c <=> co2_m  
fluxSolution(findRxnIDs(model1, 'CO2t_m'), IDMore)  % USUALLY Higher
FBAsolution.x(findRxnIDs(model1, 'CO2t_m' ))  % -0.05

% #1536  CO2t_h (0.17584),, co2_c <=> co2_h
fluxSolution(findRxnIDs(model1, 'CO2t_h'), IDMore)  % Higher
FBAsolution.x(findRxnIDs(model1, 'CO2t_h' ))  % 0.06 

% #395  RUBISC_h (0.94699), h2o_h + co2_h + rb15bp_h -> 2 h_h + 2 3pg_h 
fluxSolution(findRxnIDs(model1, 'RUBISC_h'), IDMore)  % sometimes Higher
FBAsolution.x(findRxnIDs(model1, 'RUBISC_h' ))  % 0.9470





% Rxn #46  PEPC_m, Phosphoenolpyruvate carboxylase, h2o_m + co2_m + pep_m -> h_m + pi_m + oaa_m  
fluxSolution(findRxnIDs(model1, 'PEPC_m'), IDMore)  % Some activated.
FBAsolution.x(findRxnIDs(model1, 'PEPC_m' ))  %
% #22  MALOAAt_m (-1000), oaa_c + mal__L_m <=> oaa_m + mal__L_c  
fluxSolution(findRxnIDs(model1, 'MALOAAt_m'), IDMore)  % 
FBAsolution.x(findRxnIDs(model1, 'MALOAAt_m' ))  %
% #35  OAACITt_m (-1000), oaa_c + cit_m <=> oaa_m + cit_c  
fluxSolution(findRxnIDs(model1, 'OAACITt_m'), IDMore)  % 
FBAsolution.x(findRxnIDs(model1, 'MALOAAt_m' ))  %
% #13  OAAAKGt_m (1000), akg_m + oaa_c <=> akg_c + oaa_m  

% #14  ASPTA_L_m (999.985), akg_m + asp__L_m <=> glu__L_m + oaa_m  

% #11  ASPTA_L_c (-1000), akg_c + asp__L_c <=> glu__L_c + oaa_c  

%  #1686  ASPNA1t_m (-999.99), na1_m + asp__L_m <=> na1_c + asp__L_c  

% #14  ASPTA_L_m (999.985)
% #2141  ASPGLU2_m (0.00531), glu__L_c + asp__L_m <=> glu__L_m + asp__L_c  

% #1851  NH4t_m (0.28722),h_c + nh4_c <=> h_m + nh4_m  
%  #402  GLUDH2_m (0.3841), h2o_m + nad_m + glu__L_m <=> h_m + nadh_m + nh4_m + akg_m  
% 




% #118  MTHFO1_c, 5-methyltetrahydrofolate:NAD+ oxidoreductase
fluxSolution(findRxnIDs(model1, 'MTHFO1_c'), IDMore)  % Usually activated at 0.0051
FBAsolution.x(findRxnIDs(model1, 'MTHFO1_c' ))  %
% MTHFO1_c (nad dependant) replaces MTHF02_c (5-methyltetrahydrofolate:NADP+
% oxidoreductase) with the same flux !

% G6PDH_c (Glucose 6-phosphate dehydrogenase), PGDH_c (phosphogluconate dehydrogenase)
fluxSolution(findRxnIDs(model1, 'G6PDH_c'), IDMore)
fluxSolution(findRxnIDs(model1, 'PGDH_c'), IDMore)
% Key enzymes of PPP not often activated

% #114  PGLYCP_h, Bd: 0 / 1000, Phosphoglycolate phosphatase, chloroplast
fluxSolution(findRxnIDs(model1, 'PGLYCP_h'), IDMore)

% #1444  RUBISO_h, Bd: 0 / 1000, Ribulose-1,5-bisphosphate oxygenase
fluxSolution(findRxnIDs(model1, 'RUBISO_h'), IDMore)

% #409  UREASE_m, Bd: 0 / 1000, Urease, mitochondria

fluxSolution(findRxnIDs(model1, 'UREASE_m'), IDMore) % In one case (7e), it is activated
fluxSolution(findRxnIDs(model1, 'AGMT_m'), IDMore)
fluxSolution(findRxnIDs(model1, 'ARG_m'), IDMore) % 7e

% What about NO3 uptake and metabolism ?
surfNet(model1, 'no3_c', 0, FBAsolution.x)
surfNet(model1, 'no3_c', 0, fluxSolution(:,IDMore(186)))

surfNet(model1, 'no2_c', 0, FBAsolution.x) % transport of no2 in chloroplast
surfNet(model1, 'no2_c', 0, fluxSolution(:,IDMore(186))) % Use of nitrite reductase in cytosol

surfNet(model1, 'no2_h', 0, FBAsolution.x) % Use only nitrite reductase in chloroplast, not cytosolic
surfNet(model1, 'no2_h', 0, fluxSolution(:,IDMore(186))) % Use of nitrite reductase in chloro

% #1548  NTRIR_h :  5 h_h + 3 nadph_h + no2_h -> 2 h2o_h + 3 nadp_h + nh4_h 
% #332  NOR_c : 8 h_c + 6 fdxrd_c + no2_c -> 2 h2o_c + 6 fdxox_c + nh4_c  

% TPI_c (Rxn 313) conversion of xx to dhap
fluxSolution(findRxnIDs(model1, 'TPI_c'), IDMore)
fluxSolution(313, IDMore(186))
FBAsolution.x(313)

fluxSolution(467, IDMore(186))
FBAsolution.x(467)

fluxSolution(335, IDMore(186))
FBAsolution.x(335)


model1.rxns(Ind_Comp) % Newly activated reactions when deleting IDMore(186) or ATPM_m : CEF increase, and EX_photon_w increase, 
surfNet(model1, model1.rxns(Ind_Comp))

fluxSolution(findRxnIDs(model1, 'CEF_h'), IDMore(186)) % CEF_h : 1.8870, is often activated for robust reactions
FBAsolution.x(findRxnIDs(model1, 'CEF_h'))
fluxSolution(findRxnIDs(model1, 'ATPM_c'), IDMore(186)) %ATPM_c : 0.03
surfNet(model1, 'CEF_h', 0, FBAsolution.x)
surfNet(model1, 'CEF_h', 0, fluxSolution(:, IDMore(186)) )
surfNet(model1, 'pqh2_u', 0, FBAsolution.x)
surfNet(model1, 'pqh2_u', 0, fluxSolution(:, IDMore(186)) )
surfNet(model1, 'ATPS_h', 0, FBAsolution.x)
surfNet(model1, 'ATPS_h', 0, fluxSolution(:, IDMore(186)) )
surfNet(model1, 'atp_h', 0, FBAsolution.x)
surfNet(model1, 'atp_h', 0, fluxSolution(:, IDMore(186)) )

% Use of light-dependant enzyme for chlorophyll synthesis instead of using
% nadph_h
surfNet(model1, 'DVPCHLDOR_h', 0, FBAsolution.x)
surfNet(model1, 'DVPCHLDOR_h', 0, fluxSolution(:, IDMore(186)) )

% Use of nadph to synthesize chla with (DVCHLDR_h) divinyl chlorophyllide vinyl-reductase
% instead of using PCHLDOR_h
surfNet(model1, 'DVCHLDR_h', 0, FBAsolution.x)
surfNet(model1, 'DVCHLDR_h', 0, fluxSolution(:, IDMore(186)) )

surfNet(model1, 'CHLPAS_h', 0, FBAsolution.x)
surfNet(model1, 'CHLPAS_h', 0, fluxSolution(:, IDMore(186)) )


surfNet(model1, 'DVPCHLDOR_h')
%}

% Determining the subSystems of most robust reactions
Sub_r = model1.subSystems(findRxnIDs(model1, rxn_list_red(ID_Rob)));
Sub_r = [Sub_r{:}];
R_r = unique(Sub_r);
Total_Sub_r = length(R_r);
Tot_Sub = getModelSubSystems(model1);
Prop_Sub_r = 100 * ( length(R_r)/ length(Tot_Sub) );

% Listing all reactions of most robust reactions WITHOUT BIOMASS components
Bio_Toremove = {'biomass_DNA_c', 'biomass_RNA_c',...
'biomass_carb_c', 'biomass_mem_lipids_c', 'biomass_pigm_h', 'biomass_pro_c', 'bof_c_accumulation_c',...
'DM_biomass_c_acc_c', 'biomass_c_storage_c'};

Index = find(contains(modelRob, Bio_Toremove));
modelRobBio = modelRob;
modelRobBio(Index) = [];

% Writing a csv file without BIOMASS components
modelRobBio = strcat('R_', modelRobBio);
T = cell2table(modelRobBio);
if k == 1
writetable(T,'Rxns_Rob_bio.csv');
else
writetable(T,'Rxns_Rob_bio_MOMA.csv');
end
    
% Writing a .csv file listing the most robust reactions
% Convert cell to a table and use first row as variable names
modelRob_R = strcat('R_', modelRob);
T = cell2table(modelRob_R);
if k == 1
writetable(T,'Rxns_Robust.csv');
else
writetable(T,'Rxns_Robust_MOMA.csv');
end

% Isolating all metabolites involved in robust reactions and Writing a .csv
% file listing them
MetRob = modelS(:,findRxnIDs(model1, modelRob)); % Matrix of all metabolites for robust reactions
i = 0;
MetRobVec = zeros(size(MetRob, 1), 1); % 1707 * 1
MetRobMat = []; % empty matrix of undefined;
for i = 1:size(MetRob,2)
    MetRobVec = find(MetRob(:,i) ~= 0); % Find a list of all metabolite IDs for each reaction
    [nbrow nbcol] = size(MetRobVec);
    MetRobVecNew = padarray(MetRobVec,[size(MetRob,1)-nbrow,0],0,'post');  
    MetRobMat = horzcat(MetRobMat, MetRobVecNew) ;  
end
MetRobList = unique(MetRobMat(:));
% Create a list of metabolite IDs that are involved in robust reactions
MetRobList(MetRobList == 0) = []; % Remove first zero

% Writing a .csv file listing the most robust metabolites
% Convert cell to a table and use first row as variable names
modelRobMet = strcat('M_', model1.mets(MetRobList));
T = cell2table(modelRobMet);
if k == 1
writetable(T,'Mets_Robust.csv');
else
writetable(T,'Mets_Robust_MOMA.csv'); 
end

%{
% Calculating mean, min, max and range for fluxes of all 2144 reactions for
% ROBUST reaction deletion
tol = 1E-06
% Two values, u and v, are within tolerance if abs(u-v) <= tol*max(abs(A(:)))
% Tolerance is relative by default
i = 0
max1_r = 0
min1_r = 0
range1_r = 0
U = 0
Length_equal = []
Prop_equal = []
matU = zeros(size(Flux_mat, 1), size(Flux_mat,2));
LU_r = []
for i = 1:size(Flux_mat,1)
max1_r(i) = max(Flux_mat(i,:));
min1_r(i) = min(Flux_mat(i,:));
range1_r(i) = max1_r(i) - min1_r(1);
U = uniquetol(Flux_mat(i,:), tol); % Return unique element in a row
U = U(~isnan(U));
Prop_equal(i) = length(U) / length(Flux_mat(i,:));
    %for z = 1:length(U)
    %    I = find(Flux_mat(i,:) == U(z));
    %    matU(i, 1:length(I)) = I; % A matrix of Index of deletion for each reaction
    %end
matU(i, 1:length(U)) = U;
LU_r(i) = length(U);
end
Mean_equal_r = sum(Prop_equal) / length(Prop_equal) % Mean proportion of equal fluxes within 0.1 unit for most robust reactions
SD_r = std(Prop_equal)

% To see the unqiue flux values for each reaction
Unique_Flux = sparse(matU)

% Study the variability of fluxes for reactions that are inactive for the control
% Remove the rows of Flux_mat (with MORE robust reactions only in columns) corresponding to Rxn ID of inactive reactions
Flux_mat_inac = Flux_mat(IndInac,:)
size(Flux_mat_inac)

% find the rxnID (rows number) that have fluxes different from zero
i = 0;
j = 0;
ID = zeros(size(Flux_mat_inac, 1), 1);
IDnew = 0;
for i = 1:size(Flux_mat_inac, 2)
    IDnew = find(Flux_mat_inac(:,i) ~= 0);
    [nbrow nbcol]= size(IDnew);
    IDnew = padarray(IDnew,[size(Flux_mat_inac,1)-nbrow,0],0,'post');
    ID = horzcat(ID, IDnew) ; % This is a matrix of rxn ID, for which fluxes are different from 0   
end

% Count the number of rxns with fluxes different than 0 for each column (or
% gene deletion)
ID(:,1) = []; % Delete the first column of ID with several zeros

NumbRxnNonZero = sum(ID ~= 0,1);
mean(NumbRxnNonZero)
std(NumbRxnNonZero)
range = max(NumbRxnNonZero) - min(NumbRxnNonZero)
median(NumbRxnNonZero)
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting fluxes of less robust reactions (Flux_mat)
ID_Sens = find(grRatioAct >= -0.01 & grRatioAct <= 0.01); % Fluxes of less robust reactions
length(ID_Sens);
% Remove IDs of reaction deletion yielding infeasible problem
ID_Infeasible = find(ismember(ID_Sens, NaN_ratio) == 1); % indices of active reactions that yield infeasible problem
length(ID_Infeasible);
ID_Sens(ID_Infeasible); % ID in IndSens that yield infeasible solution
IndActive(ID_Sens(ID_Infeasible)); % reaction IDs in IndActive yielding an infeasible solution
model1.rxns( IndActive(ID_Sens(ID_Infeasible)) ); % Reaction names, for which deletion is infeasible
fprintf('Number of reactions for which reaction deletion is infeasible: %4.2f \n', length(model1.rxns( IndActive(ID_Sens(ID_Infeasible)))) );
fprintf('Proportion of active reactions for which reaction deletion is infeasible: %4.2f \n',100 * length(model1.rxns( IndActive(ID_Sens(ID_Infeasible)))) / length(IndActive) );
ID_Sens(ID_Infeasible) = [] ; % Removing infeasible reaction deeltion from ID_Sens
length(ID_Sens);
fprintf('Number of less robust reactions : %4.2f \n', length(ID_Sens));
fprintf('Proportion of less robust reactions out of all active reactions: %4.2f \n', 100 * length(ID_Sens)/length(IndActive));

L_l = length(ID_Sens); 

Flux_mat = zeros(size(fluxSolutionAct,1),L_l);
z = 0;
for z = 1:L_l
        Flux_mat(:,z) = fluxSolutionAct(:,ID_Sens(z));
end

% Listing all reactions of LESS robust reactions
rxn_list_red = model1.rxns(IndActive);
modelSens = rxn_list_red(ID_Sens);

% Looking for the proportion of reactions requiring ATP or NADPH that are
% LESS robusts
IDLess = findRxnIDs(model1,modelSens); % Reaction IDs of less robust reactions
modelS = full(model1.S);

IDatpc = find(modelS(findMetIDs(model1, 'atp_c'),:) ~= 0);
IDatph = find(modelS(findMetIDs(model1, 'atp_h'),:) ~= 0);
IDatpm = find(modelS(findMetIDs(model1, 'atp_m'),:) ~= 0);
IDatpx = find(modelS(findMetIDs(model1, 'atp_x'),:) ~= 0); % No 'atp_u'

RxnID_atpc = find(ismember(IDLess, IDatpc') ~= 0); %ismember : % find if IDatpc is in IDLess
RxnID_atph = find(ismember(IDLess, IDatph') ~= 0);
RxnID_atpm = find(ismember(IDLess, IDatpm') ~= 0);
RxnID_atpx = find(ismember(IDLess, IDatpx') ~= 0);

fprintf('Proportion of less robust reactions using atp : %4.2f \n', 100 *  ( length(RxnID_atpc) + length(RxnID_atph) + length(RxnID_atpm) + length(RxnID_atpx) ) / length(IDLess));

% How the proportion of active reactions change for different reaction deletions for LESS robust reactions?
i = 0;
L = 0;
Compensatory = 0;
for i = 1:size(Flux_mat, 2)
    Flux_Vector = Flux_mat(:,i); 
    %Ind_Act = find(Flux_Vector ~= 0);
    Ind_Act = find(abs(Flux_Vector) > 1E-03); % if I use > 1E-03, no sensitive reaction deletion with more active reactions than control
    L(i) = length(Ind_Act);
    if L(i) > length(IndActive) 
        Compensatory(i) = i;
    end
end
Mean_L = sum(L) / length(L);
Std = std(L);
CV = ( Std / Mean_L ) * 100;
Median = median(L);
IQR = iqr(L);
% Proportion of reaction deletion with number of fluxes higher than in the
% control (WT)
% (WT)
length(find(L > length(IndActive) ));
length(find(L < length(IndActive) ));
fprintf('Number of reaction deletion (KO strain) with number of non-zero fluxes > number of active fluxes in control : %4.2f \n', length(find(L > length(IndActive) )));
fprintf('Number of reaction deletion (KO strain) with number of non-zero fluxes < number of active fluxes in control : %4.2f \n', length(find(L < length(IndActive) )));

% What are the only reaction deletion that can activate more reactions than in the control?
% they are exchange reaction yielding an infeasible problem, a NaN grRatio and thus an
% abberrant flux distribution
Compensatory(Compensatory == 0) = [];
modelSens(Compensatory);
printRxnFormula(model1, modelSens(Compensatory)); % One compensatory reaction (PIt_e, ID_Sens(302))
ID_Sens(Compensatory);

% So, deleting this set of reaction leads to a compensatory response
% This suggests that deletion of less robust reactions is
% coupled to a lower compensatory activation of other reactions (compared to the robust reaction case).
% Those strange distribution fluxes are probably erroneous since they are obtained for reaction deletions yielding an infeasible flux distribution...        

% Find reaction ID that are activated in KO (i.e., row in flux mat with non-zero
% fluxes, but that are not in IndActive or in control. Ncomp
% Find also reactions that are deactivated in KO relative to WT. Ncomp2
i = 0;
Ind_Act = 0;
Flux_Vector = 0;
Nact = [];
Ncomp = [];
Ncomp2 = [];
maxRange = [];
for i = 1:size(Flux_mat, 2)
    %if ismember(i, Compensatory) == 1 % Remove Compensatory reactions that are artefacts because the problem is infeasible... (grRatio = NaN)
    %%disp('infeasible')
    %else
    Flux_Vector = Flux_mat(:,i); 
    length(Flux_Vector);
    Ind_Act = find(Flux_Vector ~= 0);
    Nact(i) = length(Ind_Act);
    %Ind_Act = find(abs(Flux_Vector) > 1E-09);
    %length(Ind_Act)
    Ind_Comp = setdiff(Ind_Act, IndActive); % find newly active reaction in KO relative to WT, i.e., findindices in Ind_Act (KO active reaction IDs) that are not in IndActive (control active reaction IDs)
    rangeComp = [];
    j = 0;
    for j = 1:length(Ind_Comp)
        rangeComp(j) = range(Flux_mat(Ind_Comp(j),:));
    end
    maxRange(i) = max(rangeComp);
    Ncomp(i) = length(Ind_Comp);
    Ind_Comp2 = setdiff(IndActive, Ind_Act); % find reactions in WT that was active , but are not active anymore
    Ncomp2(i) = length(Ind_Comp2);
    assert( length(IndActive) - length(Ind_Comp2) + length(Ind_Comp) == length(Ind_Act)) % Total active reaction in control - reactions becoming non-active + new active reaction = total active reactions in KO
    %end
end
% Maximum range in fluxes for all newly activated reaction in each reaction
% deletion
maxRange;

% plot(Nact, Ncomp, 'o')
mean(Ncomp); % Newly active reactions
std(Ncomp);
CV =  std(Ncomp) / mean(Ncomp) * 100;
median(Ncomp);
IQR = iqr(Ncomp);
% plot(Nact, Ncomp2, 'o')
mean(Ncomp2); % Newly inactive reactions
std(Ncomp2);
CV = std(Ncomp2) / mean(Ncomp2) * 100;
median(Ncomp2);
IQR = iqr(Ncomp2);
% Much more reactions are inactivated than activated in KO compared to WT,
% leading to a net total number of active reactions in KO smaller than
% total number of active reaction in WT

Flux_mat(Ind_Comp, :);

%{
% Determine the reactions that are activated for both robust and sensitive
% reactions
intersect(Ind_New_Comp, Ind_Comp) % NOT good !!! Here we just compare two reaction deletion test...
length(ans)
model1.rxns( intersect(Ind_New_Comp, Ind_Comp) )
surfNet(model1, model1.rxns( intersect(Ind_New_Comp, Ind_Comp) ) )
% alanine transaminase : ALATA_L_c (in mitochondria and cytosol)
% aspartate transaminase (ASPTA_L_c) (in mitochondria and cytosol)
% several transport reactions
%  #2140  GLYCK2_m and #2139  HYPRRx_m and #2138  SPT_m: 
% pyr_m -> hpyr_m -> 2pg_m
% Energy production : CEF_h, CYOR_m
% P and N balance : DM_NO3, PYROHt_c : h2o_c + ppi_c -> h_e + 2 pi_c  
%}
    
% Large change in fluxes for those exchange reactions when deleting robust
% reactions
KO_photon = fluxSolution(findRxnIDs(model1, 'EX_photon_e'), IDLess);
WT_photon = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e'));
PropPhoton = 100 * length( find( abs(KO_photon) ./ abs( abs(WT_photon) ) > 1 ) ) / length(IDLess);
fprintf('Proportion of sensitive reaction deletion with up-regulated photon uptake : %4.2f \n', PropPhoton);

fluxSolution(findRxnIDs(model1, 'EX_no3_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_no3_e'));
fluxSolution(findRxnIDs(model1, 'DM_no3_c'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'DM_no3_c'));

% Change a little bit :
fluxSolution(findRxnIDs(model1, 'EX_o2_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_o2_e'));
fluxSolution(findRxnIDs(model1, 'EX_h2o_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_h2o_e'));
fluxSolution(findRxnIDs(model1, 'EX_h_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_h_e'));

% Does NOT CHANGE AT ALL :
fluxSolution(findRxnIDs(model1, 'EX_mg2_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_mg2_e'));
fluxSolution(findRxnIDs(model1, 'EX_so4_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_so4_e'));
fluxSolution(findRxnIDs(model1, 'EX_pi_e'), IDLess);
FBAsolution.x(findRxnIDs(model1, 'EX_pi_e'));

% Key reactions, large changes :
KO_CEF = fluxSolution(findRxnIDs(model1, 'CEF_h'), IDLess); % CEF_h is often activated for robust reactions
WT_CEF = FBAsolution.x(findRxnIDs(model1, 'CEF_h'));
PropCEF = 100 * length( find( KO_CEF > 0) )  / length(IDLess);
fprintf('Proportion of sensitive reaction deletion with up-regulated cyclic electron transport : %4.2f \n', PropCEF);

% CEF between around 0.002 and one time 2

fluxSolution(findRxnIDs(model1, 'ATPS_h'), IDLess) ;
FBAsolution.x(findRxnIDs(model1, 'ATPS_h'));
KO_ATP_h = fluxSolution(findRxnIDs(model1, 'ATPS_h'), IDLess);
WT_ATP_h = FBAsolution.x(findRxnIDs(model1, 'ATPS_h'));
PropATPh = 100 * length( find( KO_ATP_h ./ WT_ATP_h > 1) )  / length(IDLess);
fprintf('Proportion of sensitive reaction deletion with up-regulated ATPS_m : %4.2f \n', PropATPh);

fluxSolution(findRxnIDs(model1, 'ATPS_m'), IDLess) ;
FBAsolution.x(findRxnIDs(model1, 'ATPS_m'));
KO_ATP_m = fluxSolution(findRxnIDs(model1, 'ATPS_m'), IDLess);
WT_ATP_m = FBAsolution.x(findRxnIDs(model1, 'ATPS_m'));
PropATPm = 100 * length( find( KO_ATP_m ./ WT_ATP_m > 1) )  / length(IDLess);
fprintf('Proportion of sensitive reaction deletion with up-regulated ATPS_m : %4.2f \n', PropATPm);

% MOMA : 0.86

% Determining the subSystems of less robust reactions
Sub_l = model1.subSystems(findRxnIDs(model1, rxn_list_red(ID_Sens)));
Sub_l = [Sub_l{:}];
R_l = unique(Sub_l);
Total_Sub_l = length(R_l);
Tot_Sub_l = getModelSubSystems(model1);
Prop_Sub_l = 100 * ( length(R_l)/ length(Tot_Sub_l) );

% Listing all reactions of less robust reactions without biomass equations
Bio_Toremove = {'biomass_DNA_c', 'biomass_RNA_c',...
'biomass_carb_c', 'biomass_mem_lipids_c', 'biomass_pigm_h', 'biomass_pro_c', 'bof_c_accumulation_c',...
'DM_biomass_c_acc_c', 'biomass_c_storage_c'};
Index = find(contains(modelSens, Bio_Toremove));
modelSensBio = modelSens;
modelSensBio(Index) = [];

% Writing a csv file without BIOMASS components
modelSensBio = strcat('R_', modelSensBio);
T = cell2table(modelSensBio);
if k == 1
writetable(T,'Rxns_Sens_bio.csv');
else
writetable(T,'Rxns_Sens_bio_MOMA.csv');  
end

% Writing a .csv file listing the most robust reactions
% Convert cell to a table and use first row as variable names
modelSens_R = strcat('R_', modelSens);
T = cell2table(modelSens_R);
if k == 1
writetable(T,'Rxns_Sens.csv');
else
writetable(T,'Rxns_Sens_MOMA.csv');
end

% Isolating all metabolites involved in LESS robust reactions and Writing a .csv
% file listing them
MetLess = modelS(:,findRxnIDs(model1, modelSens)); % Matrix of all metabolites for robust reactions
i = 0;
MetLessVec = zeros(size(MetLess, 1), 1); % 1707 * 1
MetLessMat = []; % empty matrix of undefined;
for i = 1:size(MetLess,2)
    MetLessVec = find(MetLess(:,i) ~= 0); % Find a list of all metabolite IDs for each reaction
    [nbrow nbcol] = size(MetLessVec);
    MetLessVecNew = padarray(MetLessVec,[size(MetLess,1)-nbrow,0],0,'post');  
    MetLessMat = horzcat(MetLessMat, MetLessVecNew) ;  
end
MetLessList = unique(MetLessMat(:));
% Create a list of metabolite IDs that are involved in robust reactions
MetLessList(MetLessList == 0) = []; % Remove first zero

% Writing a .csv file listing the most robust metabolites
% Convert cell to a table and use first row as variable names
modelLessMet = strcat('M_', model1.mets(MetLessList));
T = cell2table(modelLessMet);
if k == 1
writetable(T,'Mets_Less_Robust.csv')
else
writetable(T,'Mets_Less_Robust_MOMA.csv')
end

% Display the subsystems of most and less robust reactions
% display the number of reactions per subsystems
% header = {'Subsys_Most_Robust', 'Subsys_Less_Robust'};
%R_r_vert = vertcat(cell((50-42),1), R_r')
%cell2table([R_r_vert, R_l'],'variableNames', {'Subsys_Most_Robust', 'Subsys_Less_Robust'})

% Comparison of subsystems involved in both sets of reactions
intersect(R_r, R_l, 'stable');
setdiff(R_r, R_l); % Values in R_r not in R_l
setdiff(R_l, R_r); % Values in R_l not in R_r

fprintf('Number of robust reactions : %4.2f \n', L_r);
fprintf('Number of not robust reactions : %4.2f \n', L_l);
fprintf('There are %4.2f times less robust reactions than non-robust ones \n', L_l / L_r);
fprintf('There are %4.2f times less subsystems for robust reactions than for non-robust reactions \n', length(R_l) / length(R_r));
% 33 out of 37 subsystems for the most robust reaction deletions are also
% present for the less robust reaction deletion

%{
% Calculating mean, min, max and range for all 2144 reactions (see
% proportion of Less Robust reaction that are more variable)
tol = 1E-06
i = 0
max1 = 0
min1 = 0
range1 = 0
U = 0
Length_equal = []
Prop_equal = []
LU = 0
matU = zeros(size(Flux_mat, 1), size(Flux_mat,2));
for i = 1:size(Flux_mat,1)
max1(i) = max(Flux_mat(i,:));
min1(i) = min(Flux_mat(i,:));
range1(i) = max1(i) - min1(1);
U = uniquetol(Flux_mat(i,:),tol);
U = U(~isnan(U));
Length_equal(i) = length(Flux_mat(i,:)) - length(U) ;
Prop_equal(i) = Length_equal(i) / length(Flux_mat(i,:));
matU(i, 1:length(U)) = U;
LU(i) = length(U);
end
Mean_equal_l = sum(Prop_equal) / length(Prop_equal) % Mean proportion of equal fluxes within 0.1 unit for less robust reactions
SD_l = std(Prop_equal)

sum(LU_r) % Sum of the number of unique fluxes for robust reaction deletion
sum(LU) % Sum of the number of unique fluxes for less robust reaction deletion
% A little bit more unique values for the robust set of reactions for a
% given tolerance
fprintf('%4.2f-fold more reactions have unique fluxes for the robust set of reactions at this tolerance', (sum(LU_r) / sum(LU) ))
% 1.84-fold more at a tolerance of 1E-06
%}

%{
figure (1)
plot(min1_r)
figure (2)
plot(min1)

plot(max1_r)
plot(max1)

figure (1)
plot(range1_r)
figure (2)
plot(range1)
% Much more variability for some reactions of the sensitive set of reactions
Mean_range1_r = sum(range1_r) / length(range1_r) 
SD_r = std(range1_r)
Mean_range1 = sum(range1) / length(range1) 
SD1 = std(range1)
% See Rxn 1390 and 1401 and 1408, those reactions fluxes become negative..
%model1.rxns(1390) 
%model1.rxns(1401) 
%model1.rxns(1408)
%}

%{
% Study the variability of fluxes for reactions that are inactive for the control
% Remove the rows of Fluxmat (with less robust reactions only in columns) corresponding to Rxn ID of inactive reactions
Flux_mat_inac = Flux_mat(IndInac,:)
size(Flux_mat_inac)

% find the rxnID (rows number) that have fluxes different from zero
i = 0;
j = 0;
ID = zeros(size(Flux_mat_inac, 1), 1);
IDnew = 0;
for i = 1:size(Flux_mat_inac, 2)
    IDnew = find(Flux_mat_inac(:,i) ~= 0);
    [nbrow nbcol]= size(IDnew);
    IDnew = padarray(IDnew,[size(Flux_mat_inac,1)-nbrow,0],0,'post');
    ID = horzcat(ID, IDnew) ; % This is a matrix of rxn ID, for which fluxes are different from 0   
end


% Count the number of rxns with fluxes different than 0 for each column (or
% gene deletion)
ID(:,1) = []; % Delete the first column of ID with several zeros

NumbRxnNonZero = sum(ID ~= 0,1);
mean(NumbRxnNonZero)
std(NumbRxnNonZero)
range = max(NumbRxnNonZero) - min(NumbRxnNonZero)
median(NumbRxnNonZero)
%}

end


% Time the code
toc;












%{
% How many compensatory Rxns are common for all reaction analysis?

for i = 2:size(ID, 2)
    %num = setdiff(ID(:,1), ID(:,i))
    num = intersect(ID(:,1), ID(:,i)); % find common values between col 1 and col2
    share = ismember(ID(:,1), ID(:,i), 'rows') ;% find elements that are shared between both rows
    notshared = find(share == 0); % find elements in row one that are not shared.
    ID(find(share == 1),1) = 0;  % produces a first row with non-unique elements
    % repeat this process and check at the end the unique elements in the
    % first row of the matrix
    % This shows that all elements in the first row are shared with other
    % rows.
    
end


for i = 3:size(ID, 2)
    %num = setdiff(ID(:,1), ID(:,i))
    num = intersect(ID(:,2), ID(:,i)); % find common values between col 1 and col2
    share = ismember(ID(:,2), ID(:,i), 'rows') ;% find elements that are shared between both rows
    notshared = find(share == 0); % find elements in row one that are not shared.
    ID(find(share == 1),2) = 0;  % produces a first row with non-unique elements
    % repeat this process and check at the end the unique elements in the
    % first row of the matrix
    % This shows that all elements in the second row are shared with other
    % rows.
    
end


for i = 4:size(ID, 2)
    %num = setdiff(ID(:,1), ID(:,i))
    num = intersect(ID(:,3), ID(:,i)); % find common values between col 1 and col2
    share = ismember(ID(:,3), ID(:,i), 'rows') ;% find elements that are shared between both rows
    notshared = find(share == 0); % find elements in row one that are not shared.
    ID(find(share == 1),3) = 0;  % produces a first row with non-unique elements
    % repeat this process and check at the end the unique elements in the
    % first row of the matrix
    % This shows that all elements in the second row are shared with other
    % rows.
    
end

% count unique values
[uv,~,idx] = unique(ID);
n = accumarray(idx(:),1);

uv = unique(x);
n  = histc(x,uv);




%%%%%%%%%%%%%
% Calculate at each growth rate inhibition, the percentage of reaction
% Should I use only FBAsolution.x ~= 0 ?
NaN_ratio = find(isnan(grRatio))
grRatio(NaN_ratio) = 0 % Assume that NaN values in grRatio or 'Infeasible' FBA yields a growth rate = 0

length(grRatio)

i = 0
z = 0
L1 = 0
Lactive = 0
L = 0
for i = 1:1000
j = i / 1000;
IndFlux = find(grRatio > j*0.9 | grRatio< j*1.1) % & grRatioAct < j*1.1); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
Lactive(i) = length(IndActive)
    for z = 1:length(IndFlux)
        Flux_Vector = fluxSolution(:,IndFlux(z)); % Extract a Flux vector of 2156 values for each relevant value of grRatio (z)
        Ind_Act = find(Flux_Vector ~= 0);
       L(z) = length(Ind_Act);
    end
MeanL = sum(L) / length(L)
NumberRobustRxn(i) = L1(i); % Number of robust reactions among all active reactions
NumberActiveRxn(i) = Lactive(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberActiveRxn(i);
end
figure (1)
Y = Percentage_RobustRxn;
length(Y)
X = 0.001:0.001:1;
X = [0,X]
Y = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + Y(1),Y]    
Y = fliplr(Y)
length(X)
plot(X, Y)
xlim([-0.02 1.02])
ylim([0 103])
xlabel('Growth rate inibition (KO / WT)')
ylabel('Percentage of reactions out of all active reactions')
title('Robustness of the FBA result to gene deletion')


























% Which gene deletion would lead to a lower growth rate?
plot(grRateKO)
%
FBAsolution = optimizeCbModel(model1,'max', 'one')
FBAsolutionNZ.x = find(FBAsolution.x ~= 0) % Find the number of non-zeros flux in the vector FBAsolution.x
FBAsolutionZ.x = find(FBAsolution.x == 0) % Find the number of zeros flux in the vector FBAsolution.x
length(FBAsolutionNZ.x) % number of non-zeros fluxes
length(FBAsolutionZ.x) % number of zeros fluxes

NaN_ratio = find(isnan(grRatio))
length(NaN_ratio)
length(grRatio)
length(FBAsolutionNZ.x) % number of non-zeros fluxes
length(FBAsolutionZ.x) % number of zeros fluxes
assert(length(grRatio) == length(FBAsolutionNZ.x) + length(FBAsolutionZ.x))

% There are multiple reactions with a grRatio close to 1 (ex.: 1.00002...)
% There are one reaction with a grRatio > 1.01
grRat_high = find(grRatio > 1.01)
grRatio(grRat_high)
model.rxns(grRat_high) % if we remove DM_dmsp_c, the growth rate increases by 1.2-fold

% Testing the robustness of the modeled FBA growth rate to deletion of
% reactions.
% If grRatio2 is used, the script assume that when the modeled growth rate is undefined (NaN),
% i.e., when FBA has no solutions, the cell does not grow

% Note that this analysis also assumes that the active reaction and
% inactive reaction does not change when deleting one reaction....

NaN_ratio = find(isnan(grRatio))
grRatio2 = grRatio
grRatio2(NaN_ratio) = 0 % Assuming that undefined FBA means no growth

IndFlux = find(grRatio >= 0)
IndNeg = find(grRatio < 0)
length(IndFlux) % 
assert(length(model1.rxns) == length(IndFlux) + length(IndNeg) + length(NaN_ratio)) 

i = 0
L1 = 0
LZ = 0
LNZ = 0
NumberRobustRxn = 0
NumberNonZeroRxn = 0
Percentage_RobustRxn = 0
for i = 1:1000
j = i / 1000
IndFlux = find(grRatio > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
LZ(i) = length(FBAsolutionZ.x); % Number of fluxes that are equal to zeros and hence does not affect the growth rate
LNZ(i) = length(FBAsolutionNZ.x);
NumberRobustRxn(i) = L1(i) - LZ(i); % Number of robust reactions among all active reactions
NumberNonZeroRxn(i) = LNZ(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberNonZeroRxn(i);
end
figure (6)
Y = Percentage_RobustRxn;
Y(1000)=[] % We must remove the last data in 'Percentage_RobustRxn' because at i = 1000, j = 1 and grRatio cannot be greater than 1.
length(Y)
X = 0.001:0.001:1;
X(1000)=[]
length(X)
%X = X * 100 % grRatio * 100, or the ratio of the KO strain growth rate relative to the WT strain
xlim([0 1]) %xlim([0 99.9999])
xticks([0 1]) %xticks([0 99.999])
plot(X, Y)
xlabel('Minimum ratio in growth rate of the KO strain relative to WT')
ylabel('Percentage of reactions robust to gene deletion out of all active reactions')
title('Robustness of the FBA result to gene deletion')
text(3, 37.4, 'The 1561 inactive reactions (73 % of all rxns) are not considered')
text(3, 37, 'There are 595 active reactions with non-zero fluxes')
text(3, 35.9, 'This means that a single reaction deletion for around 33% of all active rxns does not affect growth')
text(3, 35, 'For around 60% of active reactions, single reaction deletion abolishs growth rate')
% Around 100-37=63% of all active reaction when deleted (one by one)
% abolish the growth of our model cell.
% This is impressive and suggest that an important number of all quantitiatively important reactions
% can be deleted without affecting the growth rate appreciably. This is
% surprising given the fact that the network is a simplification, does not
% include metabolic feedback. This is also a prediction of reaction
% deletion and, hence, the network will be even more robust to deletion of
% one isozyme than to one reaction with more than one isozymes. This
% highlights the connectivity (resilience) among reactions of the network,
% and hence, by extension of a real cell, which possess more regulatory
% possibility and uses more reactions than those included in our
% genome-scale model.


% Computing gene deletion robustness analysis taking into account NaN
% Here grRatio with NaN are remove (using grRatio2) and hence the total number of reaction
% can be accounted for

i = 0
for i = 1:1000
j = i / 1000
IndFlux = find(grRatio2 > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
LZ(i) = length(FBAsolutionZ.x); % Number of fluxes that are equal to zeros and hence does not affect the growth rate
LNZ(i) = length(FBAsolutionNZ.x);
NumberRobustRxn(i) = L1(i) - LZ(i); % Number of robust reactions among all active reactions
NumberNonZeroRxn(i) = LNZ(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberNonZeroRxn(i);
end
figure (7)
Y = Percentage_RobustRxn;
Y(1000)=[] % We must remove the last data in 'Percentage_RobustRxn' because at i = 1000, j = 1 and grRatio cannot be greater than 1.
length(Y)
X = 0.001:0.001:1;
X(1000)=[]
length(X)
%X = X * 100 % grRatio * 100, or the ratio of the KO strain growth rate relative to the WT strain
xlim([0 1])
xticks([0 1])
plot(X, Y)
xlabel('Minimum ratio in growth rate of the KO strain relative to WT')
ylabel('Percentage of reactions robust to gene deletion out of all active reactions')
title('Robustness of the FBA result to gene deletion')
text(3, 37.4, 'The 1561 inactive reactions (73 % of all rxns) are not considered')
text(3, 37, 'There are 595 active reactions with non-zero fluxes')
text(3, 35.9, 'This means that a single reaction deletion for around 33% of all active rxns does not affect growth')
text(3, 35, 'For around 60% of active reactions, single reaction deletion abolishs growth rate')


% Repeating the computations , but considering all reactions, not only
% active reactions.
i = 0
for i = 1:1000
j = i / 1000
IndFlux = find(grRatio2 > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
LZ(i) = length(FBAsolutionZ.x); % Number of fluxes that are equal to zeros and hence does not affect the growth rate
LNZ(i) = length(FBAsolutionNZ.x);
NumberRobustRxn(i) = L1(i); % Number of robust reactions among all active reactions
NumberAllRxn(i) = LZ(i) + LNZ(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberAllRxn(i);
end
figure (6)
Y = Percentage_RobustRxn;
Y(1000)=[] % We must remove the last data in 'Percentage_RobustRxn' because at i = 1000, j = 1 and grRatio cannot be greater than 1.
length(Y)
X = 0.001:0.001:1;
X(1000)=[]
length(X)
%X = X * 100 % grRatio * 100, or the ratio of the KO strain growth rate relative to the WT strain
xlim([0 1])
xticks([0 1])
plot(X, Y)
xlabel('Minimum ratio in growth rate of the KO strain relative to WT')
ylabel('Percentage of reactions robust to gene deletion out of all active reactions')
title('Robustness of the FBA result to gene deletion')

%%%%%%%%%%
% repeating the computation again, but here considering only grRatio that
% are link to active reactions
IndActive = find(FBAsolution.x ~= 0); % 595 indices
grRatioAct = grRatio(IndActive)
i = 0
for i = 1:1000
j = i / 1000
IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
Lactive(i) = length(IndActive)
NumberRobustRxn(i) = L1(i); % Number of robust reactions among all active reactions
NumberActiveRxn(i) = Lactive(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberActiveRxn(i);
end
figure (6)
Y = Percentage_RobustRxn;
Y(1000)=[] % We must remove the last data in 'Percentage_RobustRxn' because at i = 1000, j = 1 and grRatio cannot be greater than 1.
length(Y)
X = 0.001:0.001:1;
X(1000)=[]
length(X)
%X = X * 100 % grRatio * 100, or the ratio of the KO strain growth rate relative to the WT strain
xlim([0 1])
xticks([0 1])
plot(X, Y)
xlabel('Minimum ratio in growth rate of the KO strain relative to WT')
ylabel('Percentage of reactions robust to gene deletion out of all active reactions')
title('Robustness of the FBA result to gene deletion')



%%%%%%%%%%%%%
% Testing reaction deletion robustness considering active and inactive reaction
% for each grRatio; i.e., considering each flux distribution (fluxSolution)
% obtained with the command 'singleGeneDeletion'

z = 0;
i = 0;
for i = 1:1000
j = i / 1000;
IndFlux = find(grRatio > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain

for z = 1:length(IndFlux)
    Flux_Vector = fluxSolution(:,IndFlux(z)); % Extract a Flux vector of 2156 values for each relevant value of grRatio (z)
    Ind_FBA_Active = find(Flux_Vector ~= 0);
    L1(i, z) = length(Ind_FBA_Active);
    Ind_FBA_Inactive = find(Flux_Vector == 0);
    L2(i,z) = length(Ind_FBA_Inactive);
    L3(i,z) = length(IndFlux) - length(Ind_FBA_Inactive); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted. This means : Number of robust reactions among all active reactions
    RobustRxn(i,z) = 100 * L3(i,z) / L1(i,z); % Percentage of robust reactions for each of the grRatio (> j)
end

Mean_RobustRxn(i) = sum(RobustRxn(i,:)) / length(RobustRxn(i,:));
Std(i) = std(RobustRxn(i,:));
CV(i) = ( Std(i) / Mean_RobustRxn(i) ) * 100;
end
% The error around the mean for each calculation is relatively small
range(CV(1:end-1)) % Gives the range (max - min of CV)

figure (6)
Y = Mean_RobustRxn;
Y(1000)=[] % We must remove the last data in 'Percentage_RobustRxn' because at i = 1000, j = 1 and grRatio cannot be greater than 1.
length(Y)
X = 0.001:0.001:1;
X(1000)=[]
length(X)
%X = X * 100 % grRatio * 100, or the ratio of the KO strain growth rate relative to the WT strain
xlim([0.001 1]) %xlim([0 99.9999])
xticks([0.001 1]) %xticks([0 99.999])
plot(X, Y)
xlabel('Minimum ratio in growth rate of the KO strain relative to WT')
ylabel('Mean percentage of reactions with single reaction deletion')
title('Robustness of the FBA result to gene deletion')
% This plot tells us the mean percentage of reactions out of all active reactions for which a single
% reaction deletion leads to a given range in growth rate ratio : 
% from the minimum ratio in growth rate (KO/WT) to the maximum 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
% Testing reaction deletion robustness considering active and inactive reaction
% for each grRatio; i.e., considering each flux distribution (fluxSolution)
% obtained with the command 'singleGeneDeletion'
% HERE WE compute percentage of robust reactions for each growth rate ratio
% and not the cumulative proportion...

% produce a histogram of frequency with all reactions

figure (2)
histogram(grRatio)

%%%%%%%%%%%%%%
% Plotting only active reactions
Ind_FBA_A = find(FBAsolution.x ~= 0)
grRatio3 = grRatio(Ind_FBA_A)

for i = 1:length(grRatio3)
    if grRatio3(i) > 1 & grRatio3 < 1.1 % This round 
        grRatio3(i) = 1
    elseif grRatio3(i) < 0
        grRatio3(i) = 0
    end
end

% Or avoid a loop : grRatio3(grRatio3 < 0) = 0
% grRatio3(grRatio3  < 1.1 & grRatio3 > 1) = 1

figure (3)
histogram(grRatio3, 'BinWidth', 0.1, 'BinLimits', [0,1.4] ) % 'NumBins', 5,
xlabel('Category of growth rate ratio , KO/WT')
ylabel('Frequency of each category')

% Extract frequencies and bin limits
[N,edges] = histcounts(grRatio3)
rel_freq = N ./ sum(N)
rel_freq = rel_freq * 100
bar(rel_freq)
xlabel('category of growth rate ratio (KO / WT)')
ylabel('Proportion of reactions out of all active reactions')
intervals={'0-0.2'; '0.2-0.4'; '0.4-0.6';'0.6-0.8';'0.8-1.0';'1.0-1.2';'1.2-1.4'};
set(gca,'xticklabel',intervals)

% Histogram of relative frequency
% Find numIntervals
intervalWidth = 0.2
numIntervals = 1 / ( intervalWidth / (1.4 - min(grRatio3)) )

grRatio3
numIntervals = numIntervals;
x = 0:intervalWidth:intervalWidth*7

ncount = histc(grRatio3,x);
relativefreq = ncount/length(grRatio3);

bar(x-intervalWidth/2, relativefreq,1)
xlim([-0.1 1.4])
set(gca, 'xtick', x)

% Try 2:
intervalWidth = 0.2
numIntervals = 1 / ( intervalWidth / (1.4 - min(grRatio3)) )

grRatio3
numIntervals = numIntervals;
x = -0.1:intervalWidth:intervalWidth*7

ncount = histc(grRatio3,x);
relativefreq = ncount/length(grRatio3);

bar(x-intervalWidth/2, relativefreq,1)
xlim([-0.1 1.4])
set(gca, 'xtick', x)

sum(grRatio3) / length(grRatio3)

%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gene deletion analysis 

[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model1, 'FBA')

model1 = generateRules(model1)
i = 1;
for i = 1:length(model.genes)
    disp(i)
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model1, 'FBA', model1.genes(i));
end
% The function stop working at model1.genes(11) = 149227 and model1.genes(139) =
% 184309 and model1.genes(236) = 195514
% find(model1.genes{:} == '149227')
Index = strfind(model1.genes, '149227');
find(not(cellfun('isempty',Index)))
Index = find(contains(model1.rules,'x(11)'));
Index
model1.rules(Index)

Index = strfind(model1.grRules, '149227');
Index = find(not(cellfun('isempty',Index)))
model1.grRules(Index)


Index = strfind(model1.genes, '195514');
find(not(cellfun('isempty',Index)))
Index = find(contains(model1.rules,'x(236)'));
Index

Index = strfind(model1.grRules, '195514');
find(not(cellfun('isempty',Index)))

deleteModelGenes(model,geneList{i})


%%%%%%%%%%%%%
% Calculating the effect of reaction deletion on growth rate
FBAsolution = optimizeCbModel(model1, 'max', 'one');
IndActive = find(FBAsolution.x ~= 0); % 567 indices
%IndActive = find(abs(FBAsolution.x) > 1E-09); % same result
IndInac = find(FBAsolution.x == 0); % Useful for the last part looking at inactive fluxes
assert(length(IndActive) + length(IndInac) == length(FBAsolution.x));
grRatioAct = grRatio(IndActive);
fluxSolutionAct = fluxSolution(:, IndActive);
size(fluxSolutionAct);

NaN_ratio = find(isnan(grRatioAct));
grRatioAct(NaN_ratio) = 0; % Assume that NaN values in grRatio or 'Infeasible' FBA yields a growth rate = 0

length(grRatioAct);
length(grRatioAct >= 0);
assert(length(grRatioAct) == length(grRatioAct >= 0)); % Make sure no negatiev values are there
assert(length(grRatioAct) == size(fluxSolutionAct,2));

model1.rxns(IndActive(NaN_ratio)) ; % find the reaction names of reaction deletion yielding infeasible flux distribution

% Summary of grRatio for each reaction names
horzcat(model1.rxns(IndActive), num2cell(grRatioAct));
fprintf('Number of active reactions: %4.2f \n', length(IndActive));
fprintf('Proportion of active reactions: %4.2f \n', 100 * length(IndActive) / length(FBAsolution.x));
fprintf('Number of inactive reactions: %4.2f \n', length(IndInac));
fprintf('Number of reaction with infeasible reaction deletions: %4.2f \n', length(NaN_ratio));

i = 0;
z = 0;
L1 = 0;
Lactive = 0;
L = 0;
for i = 1:1000
j = i / 1000;
IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
L1(i) = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
Lactive(i) = length(IndActive);
%    for z = 1:length(IndFlux)
%        Flux_Vector = fluxSolution(:,IndFlux(z)); % Extract a Flux vector of 2156 values for each relevant value of grRatio (z)
%        Ind_Act = find(Flux_Vector ~= 0);
%        L(z) = length(Ind_Act);
%    end
%Lactive(i) = sum(L) / length(L)
NumberRobustRxn(i) = L1(i); % Number of robust reactions among all active reactions
NumberActiveRxn(i) = Lactive(i);
Percentage_RobustRxn(i) = 100 * NumberRobustRxn(i) / NumberActiveRxn(i);
end
figure (1);
Y = Percentage_RobustRxn;
length(Y);
X = 0.001:0.001:1;
X = [0,X];
Y = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + Y(1),Y]  ;  
Y = fliplr(Y);
X = X * 100 ;% Makes X-axis in percentage
length(X);
plot(X, Y);
xlim([-2 102]);
ylim([0 103]);
xlabel('Maximum growth rate inibition (%)'); % 1 - (KO / WT) x 100
ylabel('Cumulative percentage of reactions out of all active reactions');
% title('Robustness of the FBA result to reaction deletion')
% This plot gives the the cumulative proportion of reactions out of all
% active reactions, which when deleted one by one yield a given growth rate
% inhibition or less.
%}
