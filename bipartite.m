%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating a list of metabolite and reaction nodes for creating a weighted bipartite network %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code produces a text file in which we can find a table of three
% columns
% Node 1, Node 2, Flux

% For instance, for a reaction R1 : 1 A + 2 B <-> 3 C, and for calculate reaction flux (F). 
% We have several production or consumption flux (F1 to F6) for each node-reaction
% pair:
% Node 1, Node 2, Flux
% A, R1, F1 = F * 1
% B, R1, F2 = F * 2
% R1, C, F3 = F * 3
% C, R1, F4  = F * 3
% R1, B, F5 = F * 2
% R1, A, F6 = F * 1
 
% If reaction R1 is backward (right to left), F1, F2 and F3 = 0 and flux F4, F5 et F6 have positive fluxes. 
% There are no negative fluxes
% If reaction R1 is forward (left to right), F1, F2 et F3 have positive
% flux. But, F4, F5 and F6 have fluxes = 0.

% So Node 1 or Node 2 are either metabolite or reaction name. Reaction
% names are in capital letters except sink_ and DM_xxxx reactions, which
% contains lowercase letters.


%%
run Model_setup.m

%%

%{
% Be careful, use ', instead of " in functions of the toolbox and respect
% the spaces

% Initiate Cobra Toolbox
% initCobraToolbox()

changeCobraSolver('gurobi', 'all')

% Do not read The Initial Network in .xml as it introduces error in the
% cobra model. e.g., No subsystems are present
%model = readCbModel('iLB1025.xml')

% Either, load directly the model as a .mat file
load('iLB1025bon.mat')

% model = readCbModel('iLB1027_lipid.xml')

% Since there are errors when loading the model with readCbModel()
% We must run this function. Otherwise the function getModelSubSystems does
% not work.
pti = convertOldStyleModel(pti)
model = pti
model.id = 'Adapted from iLB1025'
model.description = 'Metabolic network reconstruction of Fragilariopsis cylindrus'

%cont = contains(model.rxns(), 'EX')
%Ind = find(cont ~= 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a FBA

% bof_c	0.343 biomass_pro_c + 0.023 biomass_pigm_h + 0.325 biomass_mem_lipids_c + 0.052 biomass_tag_c + 0.005 biomass_dna_c + 0.037 biomass_rna_c + 0.215 biomass_carb_c 	->	biomass_c 
printRxnFormula(model, 'bof_c')
% Update the bof_c equation : The cellular mass ratio of pigment concentrations is now
% 0.0156 g pigments / g DW (Guérin Msc, Alderkamp et al 2011)
rxns_length1 = length(model.rxns)
mets_length1 = length(model.mets) 
model = addReaction(model, 'bof_c', 'metaboliteList', {'biomass_pro_c',...
    'biomass_pigm_h', 'biomass_mem_lipids_c', 'biomass_tag_c', 'biomass_dna_c',...
    'biomass_rna_c', 'biomass_carb_c', 'biomass_c'},...
    'stoichCoeffList', [-0.1470; -0.0156; -0.325; -0.052; -0.005;...
    -0.037; -0.043; 1],...  % Prot: -0.343 : -0.1470   . Carb -0.1964: -0.043
    'reversible', false);

% Update the biomass_pigm_h equation :
% All pigment concentration are from Guérin et al MSc thesis, Alderkamp et
% al 2011 assuming a chlc2:chlc1 ratio of 2.2 (Owens and Wold, 1986)
printRxnFormula(model, 'biomass_pigm_h')
model = addReaction(model, 'biomass_pigm_h', 'metaboliteList', {'cholphya_h',...
    'cholphyc1_h', 'cholphyc2_h', 'caro_h', 'fxanth_h',...
    'diadinx_h', 'biomass_pigm_h'},...
    'stoichCoeffList', [-0.6473; -0.1110; -0.0929; -0.3771; -0.0533;...
    -0.0347; 1],...
    'reversible', false);
rxns_length1 = length(model.rxns)
mets_length1 = length(model.mets) 

% All other biomass equations are from Levering et al (2016) assuming that
% the amino acid, DNA, RNA, lipid, sugar composition is the same.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change constraints and perform a FBA
checkObjective(model)
changeObjective(model, 'DM_biomass_c')
printUptakeBound(model);
printConstraints(model, -1000, 1000)
printFluxBounds(model);

%%%%%%%%%
% Modeling growth rate
CO2_fixed = -((0.5 / 12 * 1000) * 0.3 / 24) % in mmol / gDW / h
N_fixed =  CO2_fixed / (6.63)  % in mmol / gDW / h
P_fixed =  CO2_fixed / 106  % in mmol / gDW / h
PhotonAbs = -6.79 % (Alderkamp et al 2011) gives ratio of 10 mol C / mol photon in F. cylindrus
model = changeRxnBounds(model,'EX_photon_e', PhotonAbs, 'b');
model = changeRxnBounds(model,'EX_co2_e', CO2_fixed, 'l'); 
model = changeRxnBounds(model,'EX_co2_e', 0, 'u'); 
model = changeRxnBounds(model,'EX_no3_e', N_fixed, 'b'); 
model = changeRxnBounds(model,'EX_nh4_e', 0, 'b'); 
model = changeRxnBounds(model,'EX_pi_e', P_fixed, 'b'); 
model = changeRxnBounds(model,'EX_so4_e', -28.8, 'l');
model = changeRxnBounds(model,'EX_so4_e', 0, 'u');
model = changeRxnBounds(model,'CEF_h', 0.001, 'b');
DMSP = 1.666 * 0.56 / 24 % 1.666 mmol DMSP / g DW
model = changeRxnBounds(model,'DM_dmsp_c', DMSP, 'b'); % 
model = changeRxnBounds(model,'ATPM_c', 0.03, 'b');
model = changeRxnBounds(model,'ATPM_h', 0.03, 'b');
model = changeRxnBounds(model,'ATPM_m', 0.03, 'b');
printUptakeBound(model);
printConstraints(model, -1000, 1000)
FBAsolution = optimizeCbModel(model, 'max', 'one')
% Printing the predicted growth rate at 4 °C
fprintf('The predicted growth rate at 4 °C is %2.4f d-1 \n', FBAsolution.f * 24)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we replace lowercase letters in biomass reactions by uppercase
% letter to be consistent with all other reactions which are capitalized.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = startsWith(model.rxns, 'bio')
j = find(i ~= 0)
model.rxns(j)
% Replace 'biomass' by 'BIOMASS'
model.rxns(j(1)) = {'BIOMASS_PIGM_h'}
model.rxns(j(2)) = {'BIOMASS_MEM_LIPIDS_c'}
model.rxns(j(3)) = {'BIOMASS_TAG_c'}
model.rxns(j(4)) = {'BIOMASS_CARB_c'}
model.rxns(j(5)) = {'BIOMASS_DNA_c'}
model.rxns(j(6)) = {'BIOMASS_RNA_c'}
model.rxns(j(7)) = {'BIOMASS_PRO_c'}
%}

model = model1

% I add 'R_' before each reaction and 'M_' before each
% metabolite, so that the table generated below is clearer
model.rxns = strcat('R_', model.rxns());
model.mets = strcat('M_', model.mets());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is a loop producing the .txt table

% Calculate consumption and production fluxes for each link in the
% bipartite network

% But, sometimes, if the flux is reversible, then the reactants may become
% the products (and vice versa)... This is seen when the flux column is
% negative.

modelS = full(model.S);
% First, initialize the table of nodes and fluxes for the first reaction
i = 1
ReacInd = find(modelS(:,i) < 0) % Extract reactant rows
Reac = model.mets(ReacInd)
RxnName = cell(1,length(Reac))
RxnName(:) = model.rxns(i)
RxnName = RxnName'
Val1 = [Reac, RxnName]

z = 0;
FluxR = zeros(1, length(ReacInd));
FluxRrev = zeros(1, length(ReacInd));
for z = 1:length(ReacInd)
    if FBAsolution.x(i) > 0 %= 0 && model.lb(i) == 0 && model.ub(i) == 1000
    FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
    elseif FBAsolution.x(i) < 0 && model.lb(i) == -1000 && model.ub(i) == 1000
    FluxR(z) = 0
    FluxRrev(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
    elseif FBAsolution.x(i) == 0 && model.lb(i) == 0 && model.ub(i) == 1000
    FluxR(z) =0
    elseif FBAsolution.x(i) == 0 && model.lb(i) == -1000 && model.ub(i) == 1000
    FluxR(z) = 0  
    end
end

ProdInd = find(modelS(:,i) > 0) % Extract product rows
if isempty(ProdInd) == 1    % This make the loop begin again for EX, DM, SINK reactions with no products.
    T = table(Val1, FluxR');   % Only consider the reactants for those reactions
    T.Properties.VariableNames = {'Nodes' 'Flux'};
else 
Prod = model.mets(ProdInd)
RxnName = cell(1,length(Prod))
RxnName(:) = model.rxns(i)
RxnName = RxnName'
Val2 = [RxnName, Prod]
Val3 = vertcat(Val1, Val2)

z = 0;
FluxP = zeros(1, length(ProdInd));
FluxPrev = zeros(1, length(ProdInd));
for z = 1:length(ProdInd)
    if FBAsolution.x(i) > 0 %= 0 && model.lb(i) == 0 && model.ub(i) == 1000
    FluxP(z) = abs(FBAsolution.x(i)) * abs(modelS(ProdInd(z),i));
    elseif FBAsolution.x(i) < 0 && model.lb(i) == -1000 && model.ub(i) == 1000
    FluxP(z) = 0
    FluxPrev(z) = abs(FBAsolution.x(i)) * abs(modelS(ProdInd(z),i));
    elseif FBAsolution.x(i) == 0 && model.lb(i) == 0 && model.ub(i) == 1000
    FluxP(z) = 0
    elseif FBAsolution.x(i) == 0 && model.lb(i) == -1000 && model.ub(i) == 1000
    FluxP(z) = 0
    end
end
Flux = [FluxR FluxP];
Flux = Flux';
end

if model.lb(1) == -1000 && model.ub(1) == 1000 % This takes into account reversible reaction
   Val4 = cell(size(Val1, 1), size(Val1, 2)) % Here we swap the column of Val1 so that we produce Val4
   Val4(:,2) = Val1(:,1)
   Val4(:,1) = Val1(:,2)
   % FluxRrev
   
   Val5 = cell(size(Val2, 1), size(Val2, 2))
   Val5(:,2) = Val2(:,1)
   Val5(:,1) = Val2(:,2)
   % FluxPrev is the vector of flux 
   
   Val6 = vertcat(Val3, Val4, Val5)
   
   FluxRev = [FluxRrev FluxPrev];
   FluxRev = FluxRev';

   Flux = [Flux FluxRev]
   Flux = Flux(:)
   T = table(Val6, Flux);
   T.Properties.VariableNames = {'Nodes' 'Flux'};
   
else % for reactions that are not reversibles
   T = table(Val3, Flux);
   T.Properties.VariableNames = {'Nodes' 'Flux'};
end 

for i = 2:size(modelS,2)

ReacInd = find(modelS(:,i) < 0); % Extract reactant rows
Reac = model.mets(ReacInd);
RxnName = cell(1,length(Reac));
RxnName(:) = model.rxns(i);
RxnName = RxnName';
Val1 = [Reac, RxnName];

    z = 0;
    FluxR = zeros(1, length(ReacInd));
    FluxRrev = zeros(1, length(ReacInd));
    for z = 1:length(ReacInd)
        if FBAsolution.x(i) > 0 && model.lb(i) == 0 && model.ub(i) == 1000
        FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
        elseif FBAsolution.x(i) > 0 && model.lb(i) == -1000 && model.ub(i) == 1000
        FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
        elseif FBAsolution.x(i) < 0 && model.lb(i) == -1000 && model.ub(i) == 1000
        FluxR(z) = 0;
        FluxRrev(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
        elseif FBAsolution.x(i) == 0 && model.lb(i) == 0 && model.ub(i) == 1000
        FluxR(z) = 0;
        elseif FBAsolution.x(i) == 0 && model.lb(i) == -1000 && model.ub(i) == 1000
        FluxR(z) = 0 ; 
        elseif FBAsolution.x(i) < 0 && strncmp(string(model.rxns(i)), 'EX_', 3) == true % This filters Exchange reactions
        FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));
        elseif strncmp(string(model.rxns(i)), 'DM_', 3) == true
        FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));  
        elseif ( strcmp(model.rxns(i), 'ATPM_c') == 1 || strcmp(model.rxns(i), 'ATPM_h') == 1 || strcmp(model.rxns(i), 'ATPM_m') == 1 ) % || strcmp(model.rxns(i), 'CEF_h') == 1 % non-reversible reaction reaction with constraints
        FluxR(z) = abs(FBAsolution.x(i)) * abs(modelS(ReacInd(z),i));    
        end
    end
    
ProdInd = find(modelS(:,i) > 0); % Extract product rows
if isempty(ProdInd) == 1    % This make the loop begin again for EX, DM, SINK reactions with no products.
    Tnew = table(Val1, FluxR');   % Only consider the reactants for those reactions
    Tnew.Properties.VariableNames = {'Nodes' 'Flux'};
    T = [T ; Tnew];

elseif ( model.lb(i) == 0 && model.ub(i) == 1000 ) || ( strcmp(model.rxns(i), 'R_ATPM_c') == 1 || strcmp(model.rxns(i), 'R_ATPM_h') == 1 || strcmp(model.rxns(i), 'R_ATPM_m') == 1 ) % || strcmp(model.rxns(i), 'R_CEF_h') == 1 )  % For non-reversible reactions , including non-reversible reactions with constraints
  
Prod = model.mets(ProdInd);
RxnName = cell(1,length(Prod));
RxnName(:) = model.rxns(i);
RxnName = RxnName';
Val2 = [RxnName, Prod];
Val3 = vertcat(Val1, Val2);

    z = 0;
    FluxP = zeros(1, length(ProdInd));
    for z = 1:length(ProdInd)
        if FBAsolution.x(i) > 0
        FluxP(z) = abs(FBAsolution.x(i)) * abs(modelS(ProdInd(z),i));
        elseif FBAsolution.x(i) == 0 
        FluxP(z) = 0;
        end
    end

Flux = [FluxR FluxP];
Flux = Flux';

Tnew = table(Val3, Flux);
Tnew.Properties.VariableNames = {'Nodes' 'Flux'};
T = [T ; Tnew];


elseif model.lb(i) == -1000 && model.ub(i) == 1000 % This takes into account reversible reaction

Prod = model.mets(ProdInd); % This considers production of products in the forward sens
RxnName = cell(1,length(Prod));
RxnName(:) = model.rxns(i);
RxnName = RxnName';
Val2 = [RxnName, Prod];
Val3 = vertcat(Val1, Val2);
    
    z = 0;
    FluxP = zeros(1, length(ProdInd));
    FluxPrev = zeros(1, length(ProdInd));
    for z = 1:length(ProdInd)
        if FBAsolution.x(i) > 0 %&& model.lb(i) == 0 && model.ub(i) == 1000
        FluxP(z) = abs(FBAsolution.x(i)) * abs(modelS(ProdInd(z),i));
        FluxPrev(z) = 0
        elseif FBAsolution.x(i) < 0 %&& model.lb(i) == -1000 && model.ub(i) == 1000
        FluxP(z) = 0;
        FluxPrev(z) = abs(FBAsolution.x(i)) * abs(modelS(ProdInd(z),i));
        elseif FBAsolution.x(i) == 0 %&& model.lb(i) == 0 && model.ub(i) == 1000
        FluxP(z) = 0;
        FluxPrev(z) = 0;
        %elseif FBAsolution.x(i) == 0 %&& model.lb(i) == -1000 && model.ub(i) == 1000
        %FluxP(z) = 0
        end
    end
    
    Val4 = cell(size(Val3, 1), size(Val3, 2)); % This makes the set of nodes (Val4) for the backward reaction
    Val4(:,1) = Val3(:,2);
    Val4(:,2) = Val3(:,1);
   Val5 = vertcat(Val3, Val4); % Concatenate nodes for forward (Val3) and backward (Val4) reactions
    
   FluxRev = [FluxR FluxP]; % Using FluxR above and the new FluxP
   FluxTot = [FluxRev FluxRrev FluxPrev];
   FluxTot = FluxTot';
   Tnew = table(Val5, FluxTot);
   Tnew.Properties.VariableNames = {'Nodes' 'Flux'};
   T = [T ; Tnew];
end

end

% Save the table to a text file
writetable(T, 'Table_bipartite.txt')
% See the table : type('Table.txt')

%%%%%% Must add '#' in front of the first line of 'Table_bipartite.txt', so that the Python code works.
Str = fileread('Table_bipartite.txt');
fid = fopen('Table_bipartite.txt', 'w');
if fid == -1
  error('Cannot open file for writing.');
end
fprintf(fid,'%s','#');
fprintf(fid, '%s', Str);  % Or faster: fwrite(fid, Str, 'char');
fclose(fid);

% Time the code
toc;








%{

tic
i = 1
ReacInd = find(modelS(:,i) < 0) % Extract reactant rows
Reac = model.mets(ReacInd)
RxnName = cell(1,length(Reac))
RxnName(:) = model.rxns(i)
RxnName = RxnName'
Val1 = [Reac, RxnName]

%Table = table(Val1)
%writetable(Table, 'Table.txt')
%type('Table.txt')

ProdInd = find(modelS(:,i) > 0) % Extract product rows
Prod = model.mets(ProdInd)
RxnName = cell(1,length(Prod))
RxnName(:) = model.rxns(i)
RxnName = RxnName'
Val2 = [RxnName, Prod]
Val3 = vertcat(Val1, Val2)

T = table(Val3)
for i = 2:4%size(modelS,2)
ReacInd = find(modelS(:,i) < 0); % Extract reactant rows
Reac = model.mets(ReacInd);
RxnName = cell(1,length(Reac));
RxnName(:) = model.rxns(i);
RxnName = RxnName';
Val1 = [Reac, RxnName];

ProdInd = find(modelS(:,i) > 0); % Extract product rows
Prod = model.mets(ProdInd);
RxnName = cell(1,length(Prod));
RxnName(:) = model.rxns(i);
RxnName = RxnName';
Val2 = [RxnName, Prod];
Val3 = vertcat(Val1, Val2);

Tnew = table(Val3); %, 'RowNames',{['P',num2str(i)]})
T = [T ; Tnew];

end
toc

% Save the table to a text file
writetable(T, 'Table.txt');
% See the table : type('Table.txt')



% I need to add a column with the flux associated with each reaction...
height(T) % Number of rows in the table
flux = [1:height(T)]
flux = flux'
Table_flux = table(flux)
T = [T Table_flux]

% Use of fprintf ... Does not work well for cell array
fileID = fopen('myfile.txt','w');
fprintf(fileID, '%s\n %s\n', Reac{:}, RxnName{:}); 
fclose(fid);
type('myfile.txt')

% Write column names to text file
column_names = {'Column_1', 'Column_2'};
%v1 = rand(50,1); v2 = sin((1:50).'); v3 = [1:50]' % 0 + ['a':'z, 'A':'X'].';
fid = fopen('myfile.txt','wt');
fprintf(fid, '%s ', column_names{:});
fprintf(fid, '\n');


% Write a text file
%CT = Cell_in';
fprintf(fid,'%s\n', Reac{:}, Prod{:});
fclose(fid)

type('myfile.txt')

A = magic(2); % A matrix with 2 columns

fileID = fopen('myfile.txt','w'); % The 'w' input specifies write access
nbytes = fprintf(fileID,'%5d %5d\n',A) % Write two columns
nbytes = fprintf(fileID,'%5d \n',rowcell) 

fclose(fileID);

type('myfile.txt') % View the content of the file

Cell_in = {'a', 'b'}
fid = fopen('myfile.txt','w');
CT = Cell_in';
fprintf(fid,'%s\n', CT{:});
fclose(fid)


column_names = {'First', 'Asparagus_Quotient', 'Tuba_Percentile'};
v1 = rand(50,1); v2 = sin((1:50).'); v3 = [1:50]' % 0 + ['a':'z, 'A':'X'].';
fid = fopen('abc.txt','wt');
fprintf(fid, '%s ', column_names{:});
fprintf(fid, '\n');
block_of_data = [v1, v2, v3];
fmt = repmat('%15g ', 1, 3);
fmt(end:end+1) = '\n';
fprintf(fid, fmt, block_of_data.');   %transpose is needed to get everything right
fclose(fid);

type('abc.txt')

strs = char(randi([double('A') double('Z')], 5, 3));            % Create Variable
vTbl = table(strs, 'VariableNames',{'Strings'});                % Create ?vTbl?
todayDate = datestr(now,'mm/dd/yyyy');                          % Create Date String
datesvct = repmat(todayDate, size(vTbl,1), 1);                  % Create ?Vector? Of Date String
%vTbl = [ table(datesvct, 'VariableNames', {'VDate'})  vTbl];    % Concatenate Dates In Table

vTblo = [ table(datesvct, 'VariableNames', {'VDate'})]
vTbln = [vTblo vTbl];    % Concatenate Dates In Table

%}