%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes a local sensitivity analysis of the model output to
% changes in the 65 parameterss
% This script uses the function 'FBAlocalsens'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Load the model
run Model_setup
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define initial parameter values
% param1 vector in pg/cell (chla, chlc, Fucoxan, B_carot, DD)
param1 = [chla_P, chlc_P, Fucoxan_P, B_carot_P, DD_P];

% param2 vector of GC proportion (%), GenomeSize (number of bases) and
% RNA:DNA ratio
param2 = [GC_P, GenomeSize_P, RNA_DNA_r_P];

% param3 vector of lipid ratios in mol / gDW.
% 19 lipid species
param3 = MolRatLip_P;

% param4 vector of weight amino acid percentage (% g of a given amino acid
% out of the mass of all amino acids in gram) (% g / g)
param4 = MassPercent_P;

% param5 vector of mass ratio of bof components (total_prot, total_carb,
% total_lipid, total_DNA, total_RNA, total_pigm) (g / gDW)
param5 = [total_prot_P, total_carb_P, total_lipid_P, total_DNA_P, total_RNA_P, total_pigm_P];

% param6 vector of the 
% 1) Proportion of cabohydrate as glucan
% 2) Proportion of lipid accumulated as TAG
param6 = [Prop_glucan, Prop_TAG];

% param7 vector of the 9 sugar species (mol/ gDW)
% Glucose, galactose, mannose, Xylulose, arabinose, fucose, rhamnose,
% glucuronic acid, mannose-sulfate
param7 = sugar_mol;

% param8 : Ps value (productivity in 'mmol hco3 g DW-1 h-1')
% C:N ratio in mol/mol nd N:P ratio in mol/mol
param8 = [C_N, N_P, Ps];

% param9 : vector of carbon per cell and dry weight per cell (g / cell)
param9 = [Ccell, DWcell];

% Bridge all parameter values in one 'param' vector
param = [param1, param2, param3, param4, param5, param6, param7, param8, param9];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call the function to perform a local sensitivity
% analysis of the 67 parameters on the modeled growth rate with FBA.

u = FBAlocalsens(model1, param, 1.4, 0.4)

param_names = {'chla', 'chlc', 'Fucoxan', 'BCarot', 'DD', 'GCprop',...
    'Genomesize', 'RNADNA', 'lip1', 'lip2', 'lip3', 'lip4', 'lip5',...
    'lip6', 'lip7', 'lip8', 'lip9', 'Lip10', 'lip11', 'Lip12',...
    'Lip13', 'Lip14', 'Lip15', 'lip16', 'Lip17', 'Lip18' 'Lip19',...
    'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'hist', 'ile',...
    'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val',...
    'totalProt', 'totalCarb', 'totalLip', 'totalDNA', 'totalRNA',...
    'totalPigm', 'PropGlu', 'PropTGA', 'glu', 'galac', 'mann', 'xylu',...
    'arab', 'fuco', 'rhamn', 'glucur', 'mannSul', 'CN', 'NP', 'Ps', 'Ccell', 'DWcell'};


% This tells us that the growth rate is not sensitive to changes in any
% parameters (taken one by one) except C uptake rate.
% The growth rate is more sensitive to changes in main biomass components,
% although only for very high decrease by 100 fold for parameters
% (total_RNA. total_pigm) that are not very high and hence, when a large
% decrease in those parameters occur, the growth rate may decrease or be
% abolished !!

% Also, the sensitivity of the model to all parameters (except Ps) slightly decrease if Ps
% decreases (and the growth rate decreases). The model becomes even more sensitive to changes in Ps when Ps
% is small !


% Calculating the total sensitivity range between the results of the
% sensitivity for a 1.4-fold increase and the sensitivity for a 1.4-fold
% decreasse.

ur = abs( u(:,2) - u(:,1) ) ./ uc % relative change between the upper and lower estimate, normalized to control growth rate
urel = ur * 100.   % relative change in percent
min(urel)
max(urel)

ur = abs(u(:,2)) ./ abs(u(:,1))

%{
% In absolute
ur = abs( u(:,2) - uc ) ./ uc % Relative change for the 1.4-fold increase with respect to control growth rate
urelhigh = ur * 100. 

ur = abs( u(:,1) - uc ) ./ uc % Relative change for the 1.4-fold decrease with respect to control growth rate
urellow = ur * 100. % Relative change in percent
%}

% ATTENTION
%uc = FBAsolution.x(2143)

ur = ( u(:,2) - uc ) ./ uc % Relative change for the 1.4-fold increase with respect to control growth rate
urelhigh = ur * 100. 
minHigh = min(urelhigh)
maxHigh = max(urelhigh)

ur = ( u(:,1) - uc ) ./ uc % Relative change for the 1.4-fold decrease with respect to control growth rate
urellow = ur * 100. % Relative change in percent
minLow = min(urellow)
maxLow = max(urelhigh)

fprintf('The maximum relative decrease in growth rate in percent for the local sensitivity analysis is : %4.2f \n', min(vertcat(urelhigh, urellow)))
fprintf('The maximum relative increase in growth rate in percent for the local sensitivity analysis is : %4.2f \n', max(vertcat(urelhigh, urellow)))

% Legend as a table
Table_u = [param_names', num2cell(urellow), num2cell(urelhigh)]
% create a table
Table_u = cell2table(Table_u,'variableNames', {'Components', 'Relative_change_decrease', 'Relative_change_increase'})
% write to file
writetable(Table_u,'LocalSens_Table.xls')

% plot the growth rate sensitivity (in %) for a 1.4-fold increase in each paramater
figure(1)
yt = urelhigh;
yt(yt>0) = log10(1+yt(yt>0));
yt(yt<0) = -log10(1-yt(yt<0));

bar1 = bar(yt, 'r');
%ylim([-5 5]);
x = []; % create a vectyor of x thick values
i = 0;
for i = 1:length(param_names)
    x(i)  = i;
end
xticks(x);
xticklabels(param_names);
xtickangle(90);
ax = gca;
set(ax, 'YTick', [-2 -1 0 1 2], ...
    'YTickLabel', {'-100', '-10', '0', '10', '100'});
ax.FontSize = 8; 
xlabel('Parameters', 'FontSize', 12);
ylabel('Relative change in growth rate (%)', 'FontSize', 12); % Relative change in growth rate in response to a rise in each parameter by 1.4-fold
% set(gca,'FontSize',6)
%saveas(gcf,'LocalSens_increase.pdf');
hold on

% plot the effect of a 1.4-fold decrease in each parameter on growth rate
% (in %)
%figure(2)
yt = urellow;
yt(yt>0) = log10(1+yt(yt>0));
yt(yt<0) = -log10(1-yt(yt<0));
bar2 = bar(yt, 'b');
% ylim([-5 50]);
%set(gca,'FontSize',6)
%saveas(gcf,'LocalSens_decrease.pdf');
hold off
legend({'Increase','Decrease'}, 'FontSize', 11, 'Location','northwest')
ax = gca;
set(ax, 'YTick', [-2 -1 0 1 2], ...
    'YTickLabel', {'-100', '-10', '0', '10', '100'});
ax.FontSize = 8; 
xlabel('Parameters', 'FontSize', 12);
ylabel('Relative change in growth rate (%)', 'FontSize', 12); % Relative change in growth rate in response to a rise in each parameter by 1.4-fold
yrule = ax.YAxis;
yrule.FontSize = 12
saveas(gcf,'LocalSens.pdf');


% Repeat the figure, but for only  y-axis < 5%







% find parameters with sensitivity higher than a given threshold
% Ind = find(urel > 0.5); % most sensitive parameter indices

% Find the 9th most sensitive parameters
MaxVal = maxk(urel, 9); % Max relative differences
Ind = find(ismember(urel, MaxVal) ~= 0);
urel(Ind); % sensitivity of most sensitive parameters
param_names{Ind};

%{
c1 = combnk(Ind,2)
c2 = combnk(Ind,3)
c3 = combnk(Ind,4)
c4 = combnk(Ind,5)
c5 = combnk(Ind,6)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing a serie of FBA with different set of random parameters
% increasing or decreasing by 1.4-fold
% Testing whether or not synergistic interactions occur among parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intiialize the random number generator to make the results repeatable
rng(0,'twister');

% Setting up a matrix of low and high values for each parameters
factor = 1.4;
inc = 0.4;
param_Control = param    ;
param_low = param_Control / factor;
param_low( find(param_Control == 0) ) = 0;
param_high = param_Control * factor;
param_high( find(param_Control == 0) ) = inc;
param_mat = zeros(2,length(param_names));
param_mat(1,:) = param_low;
param_mat(2,:) = param_high;


%{
% Create a random vector of 1 or 2 (i.e., low or high values)
i = 0
arrayA = []
for i = 1:65
    arrayA(i) = randi(2);
end
arrayA

% Create a random combination of parameters (increased or decreased by a
% given factor)
i = 0
param_rand = []
for i = 1:length(arrayA)
    param_rand(i) = param_mat(arrayA(i),i)
end
%}

% perform multiple FBA and test whether or not the growth rate vary by more
% than xx %
u_control = FBAlight(param, model1)

tic
i = 0;
u = [];
% Here we can perform 100 (32 sec) or more FBA (1000 : 330 sec or 5-6
% min)...
for i = 1:100000 
    % Create a random vector of 1 or 2 (i.e., low or high values)
    j = 0;
    arrayA = [];
    for j = 1:length(param_names)
        arrayA(j) = randi(2);
    end
    arrayA;

    % Create a random combination of parameters (increased or decreased by a
    % given factor)
    k = 0;
    param_rand = [];
    for k = 1:length(arrayA)
        param_rand(k) = param_mat(arrayA(k),k);
    end
    
    % Optional, if param(65) remains constant...
    %param_rand(65) = param(65);
    
    % Calculate multiple FBA
    u(i) = FBAlight(param_rand,model1);  
end
toc

% Compute standard deviation and mean
std(u);
mean(u);
% compute coefficient of variation
CV = std(u) / mean(u) * 100;

% compute relative error (in percent)
urelrand = 100 * (u - u_control) / u_control;
urelrandmin = min(urelrand);
urelrandmax = max(urelrand);

% plot the growth rate sensitivity (in %) for a random 1.4-fold increase or decrease in each paramater
figure(3)
bar(urelrand);
% ylim([-5 50]);
xlabel('Number of sensitivity analyses'); 
ylabel('Relative change in growth rate (%)'); % Relative change in growth rate in response to a random rise or decrease in each parameter by 1.4-fold

fprintf('The maximum relative decrease in growth rate for multiple random local sensitivity analysis is : %4.2f percent\n', urelrandmin)
fprintf('The maximum relative increase in growth rate for multiple random local sensitivity analysis is : %4.2f percent\n', urelrandmax)

% The above analysis shows that the error in the predicted growth rate
% ranges from 28.6% (min(urel)) to 56.2% (max(urel)). The sensitivity of
% growth rate appears to be  bit higher when multiple parameters are varied in
% combination than when parameters are varied one by one... However, it
% appears that there are no large synergistic interactions among parameters
% affecting growth rate.
% It also appears that much of the variation in growth rate occurs because of changes in C uptake rate.
% If we let constant "param(65)" or C uptake rate, the growth rate varied
% little, i.e., between 0.0053% and 7.75% for a factor of 1.5 !

% Calculating the number of combinations per binary vectors of x elements
 c = [0 1];
 n = 9; % number of parameters
 vec = repelem(2,n);
 cc = c(fullfact(vec)); % with 0 and 1
 cc2 = fullfact(vec); % with 1 and 2
 size(cc,1); % number of binary combinations
% there are more than 33 millions combinations for a binary vector of 25
% elements (n = 25)...
fprintf('There are %4.0f combinations for a binary vector of %1.0f elements \n', size(cc,1), n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the sensitivity of a few key parameters in combinations
% Here we consider all possible combinations for the most sensitive
% parameters excluding CO2 uptake rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up a matrix of low and high values for each parameter (param_mat
% with 2 rows and 65 columns)
factor = 1.4;
inc = 0.4;
param_Control = param;
param_low = param_Control / factor;
param_low( find(param_Control == 0) ) = 0;
param_high = param_Control * factor;
param_high( find(param_Control == 0) ) = inc;
param_mat = zeros(2,length(param_names));
param_mat(1,:) = param_low;
param_mat(2,:) = param_high;

% perform multiple FBA and test whether or not the growth rate vary by more
% than xx %
u_control = FBAlight(param, model1);

% Create a matrix of 1 or 2 (i.e., low or high values) for all
% combinations of most sensitive parameters
n = length(Ind); % number of parameters 
vec = repelem(2,n);
cc = fullfact(vec); % binary matrix with 1 and 2 (512 rows by 9 columns)
size(cc,1) % number of binary combinations 

tic
i = 0;
u = [];
param_cmb = [];
% Here we can perform multiple FBA
for i = 1:size(cc,1)
  param_cmb = cc(i,:); % vectors of 9 most sensitive parameters. combination of parameters between parameter 48 and 53
        
        param_cmbAll = param;
        
        k = 0;
        for k = 1:length(Ind)-1
            param_cmbAll(Ind(k)) = param_mat(param_cmb(k), Ind(k));
        end
        % param_cmbAll(65) automatically = param(65) because param_cmbAll = param at the start
        % CO2 uptake is constant !
        
        % Calculate one FBA for each parameter combination
  u(i) = FBAlight(param_cmbAll,model1);  
end
toc

% Compute standard deviation and mean
std(u);
mean(u);
% compute coefficient of variation
CV = std(u) / mean(u) * 100;

% compute relative error (in percent)
urel = 100 * (u - u_control) / u_control;
min(urel);
max(urel);

fprintf('The maximum relative decrease in growth rate in percent for all combinations of most sensitive parameters except C uptake rate : %4.2f \n', min(urel))
fprintf('The maximum relative increase in growth rate in percent for all combinations of most sensitive parameters except C uptake rate : %4.2f \n', max(urel))

% Here we see that the error remains low and that the interaction among
% sensitive parameters is low (error between 0.2899 and 7.28 %)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the sensitivity of all comibinations of most sensitive parameters
% including CO2 uptake rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a matrix of 1 or 2 (i.e., low or high values) for all
% combinations of parameters param(48) to param(53) and param(65)
n = length(Ind); % number of parameters ( from pram(48) to param(53) )
vec = repelem(2,n);
cc = fullfact(vec); % binary matrix with 1 and 2
size(cc,1); % number of binary combinations 

tic
i = 0;
u = [];
param_cmb = [];
% Here we can perform multiple FBA
for i = 1:size(cc,1)
  param_cmb = cc(i,:); % combination of parameters among parameters in 'Ind'
        
        param_cmbAll = param;
        
        k = 0;
        for k = 1:length(Ind) % Last indices should be 65
            param_cmbAll(Ind(k)) = param_mat(param_cmb(k), Ind(k));
        end
        % Here CO2 uptake IS NOT constant
        
        % Calculate one FBA for each parameter combination
  u(i) = FBAlight(param_cmbAll,model1);  
end
toc

% Compute standard deviation and mean
std(u);
mean(u);
% compute coefficient of variation
CV = std(u) / mean(u) * 100;

% compute relative error (in percent)
urel = 100 * (u - u_control) / u_control;
min(urel);
max(urel);

fprintf('The maximum relative decrease in growth rate in percent for all combinations of most sensitive parameters : %4.2f \n', min(urel))
fprintf('The maximum relative increase in growth rate in percent for all combinations of most sensitive parameters : %4.2f \n', max(urel))

% Here we see that the interactions among parameters is small, the error on
% the predicted growth rate varied by 28% to 60.9%, which is close to the
% error associated with C uptake rate alone.

% Time the code
toc;

