%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shadow prices, reduced costs as well as reaction/gene deletion analysis
%
% For reaction/gene deletion analysis, the file solveCobraQP.m must be
% modified since the MOMA method does not work otherwise. Appears to have a bug with gurobi...
% I commented the lines 600-601-602 in the function solveCobraQP.m (in
% '/Users/mlavoie/cobratoolbox/src/base/solvers')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
tic;

% load the model
run Model_setup

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use single-reaction deletion with FBA: Compute FBA when removing one
% reaction at a time.
% This gives an idea of the global robustness of the network to
% reaction deletion 

% The model.setup.m file considers constraints on the uptake of N and other elements. 
% This is needed. Otherwise, N exchange is frequently
% at -1000 mmol / gDW / h when a reaction is deleted, which is unrealistic.
% check urea and glycolate production...
% A flux of -1000 mmol NO3- / h-1 / gDW = 15 000 fmol / cell / h (* 2 g DW/gC
% and 7.5 x 10-12 gC/cell), this is much higher than what is expected for
% diatoms (Lomas and Gilbert, 2000) and even for N-limited diatoms (Zevenboom and Mur, 1981)

FBAsolution = optimizeCbModel(model1, 'max', 'one');
printFluxVector(model1, FBAsolution.x, false, true);
printConstraints(model1, -1000, 1000)

% Constraints the photon uptake rate so that no unrealistic photon uptake
% occur.
fold_change = 1.5;
Pho_low = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) * fold_change ; % not 5
Pho_high = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) / fold_change ;

model1 = changeRxnBounds(model1,'EX_photon_e', Pho_low, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_photon_e', P_low, 'u'); % Bd: -1000 / 0

FBAsolution = optimizeCbModel(model1, 'max', 'one');
% printFluxVector(model1, FBAsolution.x, false, true);

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
assert(length(grRatioAct) == length(grRatioAct >= 0)); % Make sure no negative values are there
assert(length(grRatioAct) == size(fluxSolutionAct,2));

model1.rxns(IndActive(NaN_ratio)) ; % find the reaction names of reaction deletion yielding infeasible flux distribution

grRatioAct(find(grRatioAct < 0)) = 0
grRatioAct(find(grRatioAct > 1)) = 1
grRatioAct(find(grRatioAct < 0.001 & grRatioAct > 0)) = 0

% Summary of grRatio for each reaction names
horzcat(model1.rxns(IndActive), num2cell(grRatioAct));
fprintf('Number of active reactions: %4.2f \n', length(IndActive));
fprintf('Proportion of active reactions: %4.2f \n', 100 * length(IndActive) / length(FBAsolution.x));
fprintf('Number of inactive reactions: %4.2f \n', length(IndInac));
fprintf('Number of reaction with infeasible reaction deletions: %4.2f \n', length(NaN_ratio));

i = 0;
j = 0;
L1 = 0;
Lactive = 0;
Percentage_RobustRxn = 0;
for i = 1:1000
    j = i ./ 1000;
    IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
    L1 = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
    Lactive = length(IndActive); 
    Percentage_RobustRxn(i) = 100 * L1 ./ Lactive;
end

if k == 1
    figure (1);
end

X = Percentage_RobustRxn;
length(X);
Y = 1/1000:1/1000:1;
Y = [0,Y];
X = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + X(1),X]  ;  
X = fliplr(X);
X(end) = 100 ;  % This correct for the 3.1 % of reactions for which deletion abolishes the growth rate, but yield negative grRatioAct values.
Y = Y * 100 ; % Makes X-axis in percentage
length(Y);
plot(X, Y, 'LineWidth', 1.5);
ylim([-2 102]);
xlim([0 103]);
ylabel('Maximum growth rate inibition (%)', 'FontSize', 12); % 1 - (KO / WT) x 100
xlabel('Cumulative percentage of reactions out of all active reactions (%)', 'FontSize', 12);
ax = gca;
ax.FontSize = 10; 
% title('Robustness of the FBA result to reaction deletion')
% This plot gives the the cumulative proportion of reactions out of all
% active reactions, which when deleted one by one yield a given growth rate
% inhibition or less.

%if k == 1
%saveas(gcf,'RxnDeletion_FBA.pdf');
%else
%saveas(gcf,'RxnDeletion_MOMA.pdf');
%end

if k == 1
    hold on
else
    hold off
    legend({'FBA','MOMA'}, 'FontSize', 11, 'Location','northwest')
    saveas(gcf,'RxnDeletion_FBA_MOMA.pdf');
end

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
T2 = [modelRob, model1.rxnNames(findRxnIDs(model1, modelRob)), model1.grRules(findRxnIDs(model1, modelRob))]
T2 = cell2table(T2, 'variableNames', {'Reaction_Abbr', 'Reaction_Name', 'ProteinID'})
if k == 1
writetable(T,'Rxns_Robust.csv');
writetable(T2,'Rxns_RobustFBA.xls');
else
writetable(T,'Rxns_Robust_MOMA.csv');
writetable(T2,'Rxns_RobustMOMA.xls');
end

% Isolating all metabolites involved in robust reactions and Writing a .csv
% file listing them
MetRob = model1.S(:,findRxnIDs(model1, modelRob)); % Matrix of all metabolites for robust reactions
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Writing a .csv file listing the most sensitive reactions
% Convert cell to a table and use first row as variable names
modelSens_R = strcat('R_', modelSens);
T = cell2table(modelSens_R);
T2 = [modelSens, model1.rxnNames(findRxnIDs(model1, modelSens)), model1.grRules(findRxnIDs(model1, modelSens))]
T2 = cell2table(T2, 'variableNames', {'Reaction_Abbr', 'Reaction_Name', 'ProteinID'})
if k == 1
writetable(T,'Rxns_Sens.csv');
writetable(T2,'Rxns_SensFBA.xls');
else
writetable(T,'Rxns_Sens_MOMA.csv');
writetable(T2,'Rxns_SensMOMA.xls');
end

% Isolating all metabolites involved in LESS robust reactions and Writing a .csv
% file listing them
MetLess = model1.S(:,findRxnIDs(model1, modelSens)); % Matrix of all metabolites for robust reactions
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


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gene deletion analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute gene robustness without (model1) or with (model2) allelic variants

% First compute the cases without considering allelic variants (model1)

methodDel = {'FBA', 'MOMA'};
k = 0;
for k = 1:length(methodDel)
if k == 1
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('First, compute the effect of single gene deletion using the FBA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model1, methodDel{k});
else
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('Second, compute the effect of single gene deletion using the MOMA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model1, methodDel{k});
end

%%%%%%%%%%%%%
% plotting the effect of gene deletion on growth rate
FBAsolution = optimizeCbModel(model1, 'max', 'one');
IndActive = find(FBAsolution.x ~= 0); % 567 indices
%IndActive = find(abs(FBAsolution.x) > 1E-09); % same result

% find a list of genes associated with active reactions
i = 0;
GeneActList = {};
for i = 1:length(IndActive)
    model1.grRules(IndActive(i));
    expression = '\d*\d*';
    matchNumber = regexp(model1.grRules(IndActive(i)),expression,'match');
    GeneActListNew = [matchNumber{:}];
    GeneActList = {GeneActList, GeneActListNew};
    GeneActList = [GeneActList{:}];
    %flatCellArray = flatCellArray'  
end
GeneActList = unique(GeneActList);
GeneActList = GeneActList';

% Find the reactions associated with active genes
i = 0;
ActRxns = [];
for i = 1:length(GeneActList)
    expression = GeneActList(i);
    matchNumber = regexp(model1.grRules,expression,'match');
    Ind = find(~cellfun(@isempty,matchNumber));
    %matchNumber{Ind};
    ActRxns = [ActRxns, Ind'];
end
ActRxns = unique(ActRxns) % This gives a vector of active Rxns IDs.
% Same than IndActive above , but with more reactions : Some active genes
% are with active and inactive reactions...

% find grRatio elements that are associated with active genes
GeneTot = cell2mat(model1.genes);
GeneTot = str2num(GeneTot);
GeneActList = cell2mat(GeneActList);
GeneActList = str2num(GeneActList);
Ind_grRatioAct = find(ismember(GeneTot, GeneActList) == 1);
length(Ind_grRatioAct);
assert(length(Ind_grRatioAct) == length(GeneActList));

grRatioAct = grRatio(Ind_grRatioAct);

fluxSolutionAct = fluxSolution(:, Ind_grRatioAct);
size(fluxSolutionAct);

NaN_ratio = find(isnan(grRatioAct));
length(NaN_ratio)
grRatioAct(NaN_ratio) = 0; % Assume that NaN values in grRatio (due to missing gene) yields a growth rate = 0

length(grRatioAct);
length(grRatioAct >= 0);
assert(length(grRatioAct) == length(grRatioAct >= 0)); % Make sure no negative values are there
assert(length(grRatioAct) == size(fluxSolutionAct,2));

% if grRatioAct > 1, grRation == 1
grRatioAct(find(grRatioAct > 1)) = 1
grRatioAct(find(grRatioAct < 0)) = 0
grRatioAct(find(grRatioAct < 0.001 & grRatioAct > 0)) = 0

% Summary of grRatio for each reaction names
fprintf('Number of active genes: %4.2f \n', length(grRatioAct));
fprintf('Proportion of active genes: %4.2f \n', 100 * length(grRatioAct) / length(grRatio));
fprintf('Number of inactive genes: %4.2f \n', length(grRatio) - length(grRatioAct));
fprintf('Number of gene with infeasible solutions: %4.2f \n', length(NaN_ratio));

i = 0;
j = 0;
L1 = 0;
Lactive = 0;
Percentage_Robust = [];
for i = 1:1000
    j = i ./ 1000;
    IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
    L1 = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
    Lactive = length(grRatioAct);
    Percentage_Robust(i) = 100 * L1 ./ Lactive;
end

%if k == 1
%    subplot(2,1,1)
%else
%    subplot(2,1,2)
%end

if k == 1
    figure (1);
end

X = Percentage_Robust;
length(X);
Y = 0.001:0.001:1;
Y = [0,Y]; % Add a zero at the right
X = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + X(1),X]  ;  
X = fliplr(X);
Y = Y * 100 ;% Makes Growth inhibition in percentage
length(Y);
plot(X, Y, 'LineWidth',1.5);
ylim([-3 103]);
xlim([0 103]);
ylabel('Maximum growth rate inibition (%)', 'FontSize', 12); % 1 - (KO / WT) x 100
xlabel('Cumulative percentage of genes out of all active genes (%)', 'FontSize', 12);
% title('Robustness of the FBA result to reaction deletion')
% This plot gives the the cumulative proportion of genes out of all
% genes involved in active reactions, which when deleted one by one yield a given growth rate
% inhibition or less.
ax = gca;
ax.FontSize = 10; 

%if k == 1
%saveas(gcf,'GeneDeletion_FBA.pdf');
%else
%saveas(gcf,'GeneDeletion_MOMA.pdf');
%end

if k == 1
    text(30,95, 'A', 'FontWeight', 'Bold', 'Fontsize', 16);
    hold on
else
    hold off
    legend({'FBA','MOMA'}, 'FontSize', 11, 'Location','northwest')
    saveas(gcf,'GeneDeletion_FBA_MOMA_Without.pdf');
end

%if k == 1
%    text(5,90, 'A: FBA', 'FontWeight', 'Bold', 'Fontsize', 16);
%else
%    text(5,90, 'B: MOMA', 'FontWeight', 'Bold', 'Fontsize', 16);
%    saveas(gcf,'GeneDeletion_FBA_MOMA_Without.pdf');
%end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute gene robustness with allelic variants (model2)
FBAsolution = optimizeCbModel(model2, 'max', 'one');
printFluxVector(model2, FBAsolution.x, false, true);

fold_change = 1.5;
Pho_low = FBAsolution.x(findRxnIDs(model2, 'EX_photon_e')) * fold_change ; % not 5
Pho_high = FBAsolution.x(findRxnIDs(model2, 'EX_photon_e')) / fold_change ;

model2 = changeRxnBounds(model2,'EX_photon_e', Pho_low, 'l'); % Bd: -1000 / 0
model2 = changeRxnBounds(model2,'EX_photon_e', P_low, 'u'); % Bd: -1000 / 0

FBAsolution = optimizeCbModel(model2, 'max', 'one');
% printFluxVector(model1, FBAsolution.x, false, true);

methodDel = {'FBA', 'MOMA'};
k = 0;
for k = 1:length(methodDel)
if k == 1
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('First, compute the effect of single gene deletion using the FBA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model2, methodDel{k});
else
    fprintf('------------------------------------------------------------------------------------- \n')
    fprintf('Second, compute the effect of single gene deletion using the MOMA method \n')
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model2, methodDel{k});
end

%%%%%%%%%%%%%
% plotting the effect of gene deletion on growth rate
FBAsolution = optimizeCbModel(model2, 'max', 'one');
IndActive = find(FBAsolution.x ~= 0); % 567 indices
%IndActive = find(abs(FBAsolution.x) > 1E-09); % same result

% find a list of genes associated with active reactions
i = 0;
GeneActList = {};
for i = 1:length(IndActive)
    model2.grRules(IndActive(i));
    expression = '\d*\d*';
    matchNumber = regexp(model2.grRules(IndActive(i)),expression,'match');
    GeneActListNew = [matchNumber{:}];
    GeneActList = {GeneActList, GeneActListNew};
    GeneActList = [GeneActList{:}];
    %flatCellArray = flatCellArray'  
end
GeneActList = unique(GeneActList);
GeneActList = GeneActList';

% Find the reactions associated with active genes
i = 0
ActRxns = []
for i = 1:length(GeneActList)
    expression = GeneActList(i);
    matchNumber = regexp(model2.grRules,expression,'match');
    Ind = find(~cellfun(@isempty,matchNumber));
    %matchNumber{Ind};
    ActRxns = [ActRxns, Ind'];
end
ActRxns = unique(ActRxns) % This gives a vector of active Rxns IDs.
% Same than IndActive above , but with more reactions : Some active genes
% are with active and inactive reactions...

% find grRatio elements that are associated with active genes
GeneTot = cell2mat(model2.genes);
GeneTot = str2num(GeneTot);
GeneActList = cell2mat(GeneActList);
GeneActList = str2num(GeneActList);
Ind_grRatioAct = find(ismember(GeneTot, GeneActList) == 1);
length(Ind_grRatioAct);
assert(length(Ind_grRatioAct) == length(GeneActList));

grRatioAct = grRatio(Ind_grRatioAct);

fluxSolutionAct = fluxSolution(:, Ind_grRatioAct);
size(fluxSolutionAct);

NaN_ratio = find(isnan(grRatioAct));
length(NaN_ratio)
grRatioAct(NaN_ratio) = 0; % Assume that NaN values in grRatio (due to missing gene) yields a growth rate = 0

length(grRatioAct);
length(grRatioAct >= 0);
assert(length(grRatioAct) == length(grRatioAct >= 0)); % Make sure no negative values are there
assert(length(grRatioAct) == size(fluxSolutionAct,2));

% if grRatioAct > 1, grRation == 1
grRatioAct(find(grRatioAct > 1)) = 1
grRatioAct(find(grRatioAct < 0)) = 0
grRatioAct(find(grRatioAct < 0.001 & grRatioAct > 0)) = 0

% Summary of grRatio for each reaction names
fprintf('Number of active genes: %4.2f \n', length(grRatioAct));
fprintf('Proportion of active genes: %4.2f \n', 100 * length(grRatioAct) / length(grRatio));
fprintf('Number of inactive genes: %4.2f \n', length(grRatio) - length(grRatioAct));
fprintf('Number of gene with infeasible solutions: %4.2f \n', length(NaN_ratio));

i = 0;
j = 0;
L1 = 0;
Lactive = 0;
Percentage_Robust = [];
for i = 1:1000
    j = i ./ 1000;
    IndFlux = find(grRatioAct > j); % Find the indices in grRatio, grRatio is the ratio of growth rate for the KO strain relative to the growth rate of the WT strain
    L1 = length(IndFlux); % Number of fluxes that affect the modeled growth rate by less than j * 100 % when deleted
    Lactive = length(grRatioAct);
    Percentage_Robust(i) = 100 * L1 ./ Lactive;
end

if k == 1
    figure (1);
end

X = Percentage_Robust;
length(X);
Y = 0.001:0.001:1;
Y = [0,Y];
X = [ (100. * length(find(grRatioAct == 0)) ./ length(grRatioAct)) + X(1),X]  ;  
X = fliplr(X);
Y = Y * 100 ;% Makes X-axis in percentage
length(Y);
plot(X, Y, 'LineWidth', 1.5);
ylim([-3 103]);
xlim([0 103]);
ylabel('Maximum growth rate inibition (%)', 'FontSize', 12); % 1 - (KO / WT) x 100
xlabel('Cumulative percentage of genes out of all active genes (%)', 'FontSize', 12);
% title('Robustness of the FBA or MOMA result to gene deletion')
ax = gca;
ax.FontSize = 10; 

%if k == 1
%saveas(gcf,'GeneDeletion_Variant_FBA.pdf');
%else
%saveas(gcf,'GeneDeletion_Variant_MOMA.pdf');
%end
%saveas(gcf,'GeneDeletion_Variant_FBA_MOMA.pdf');

if k == 1
    text(30,95, 'B', 'FontWeight', 'Bold', 'Fontsize', 16);
    hold on
else
    hold off
    legend({'FBA','MOMA'}, 'FontSize', 11, 'Location','northwest')
    saveas(gcf,'GeneDeletion_Variant_FBA_MOMA.pdf');
end

end

% Time the code
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
