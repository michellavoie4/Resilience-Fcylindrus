%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global sensitivity analysis of the metabolic model using the method of
% Morris or Elmentary Effect Test (EET)
% This method is efficient and is particulalry adapted for model with a
% large number of inputs.

% The method is described in Pianosi et al. (2015) and implemented
% in the Sensitivity Analysis For Everyone (SAFE) Toolbox developed 
% by F. Pianosi, F. Sarrazin and T. Wagener at Bristol University.
% See Pianosi et al (2015) for details to acquire for free the SAFE toolbox.

% Pianosi, F.; Sarrazin, F.; Wagener, T., 2015. A Matlab toolbox for global sensitivity analysis. Environ. Model. Software 70, 80-85.

% This script uses the function 'FBAlight.m' and the script 'Model_setup.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: set paths
clear all

my_dir = '/Users/mlavoie/Documents/MATLAB/safe_R1.1' ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/..../Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2: setup the model and define input ranges

% load the model
run Model_setup.m

% Define parameter values and probable ranges
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
param8 = [C_N, N_P, Ps];

% param9 : g carbon per cell (Ccell) and g dry weight per cell (DWcell)
param9 = [Ccell, DWcell];

% Bridge all parameter values in one 'param' vector
param = [param1, param2, param3, param4, param5, param6, param7, param8, param9];

% Define input distribution and ranges:
M  = length(param) ; % number of uncertain parameters [ Sm beta alfa Rs Rf ]
DistrFun  = 'unif'  ; % Parameter distribution

i = 0
DistrPar = cell(length(param),1)
factor = 1.4
for i = 1:length(DistrPar)
 DistrPar([i],1) = {[param(i) / factor, param(i) * factor]}
end

% Set ranges for parameter equal to 0
IndZero = find(param == 0) % Indices of param for which the initial parameters equal 0
inc = 0.4

i = 0
for i = 1:length(IndZero)
    DistrPar([IndZero(i)],1) = {[param(IndZero(i)), param(IndZero(i)) + inc]}  
end
DistrPar = DistrPar' % gives ranges for each parameter                                                                                                                                                                     


% Plot results:
param_names = {'chla', 'chlc', 'Fucoxan', 'BCarot', 'DD', 'GCprop',...
    'Genomesize', 'RNADNA', 'lip1', 'lip2', 'lip3', 'lip4', 'lip5',...
    'lip6', 'lip7', 'lip8', 'lip9', 'Lip10', 'lip11', 'Lip12',...
    'Lip13', 'Lip14', 'Lip15', 'lip16', 'Lip17', 'Lip18' 'Lip19',...
    'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'hist', 'ile',...
    'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val',...
    'totalProt', 'totalCarb', 'totalLip', 'totalDNA', 'totalRNA',...
    'totalPigm', 'PropGlu', 'PropTGA', 'glu', 'galac', 'mann', 'xylu',...
    'arab', 'fuco', 'rhamn', 'glucur', 'mannSul', 'CN', 'NP', 'Ps', 'Ccell', 'DWcell'};

X_Labels = param_names %{'chla', 'chlc', 'Fucoxan', 'B_carot', 'DD', 'GC_prop'} ;

% Number of uncertain parameters subject to SA:                                
M    = length(param) ; 
% Parameter ranges (from literature):
% Here I create two vectors of parameter values (xmin and xmax)
format long; % 'format long' makes it easier to see the very small numbers
i = 0;
min = zeros(1,length(param));
max = zeros(1,length(param));
for i = 1:length(DistrPar)
    val = DistrPar([1],i);
    valex = val{1,1};
    xmin(i) = valex(1);
    xmax(i) = valex(2);
end
                                                   
% Parameter distributions:                                                     
DistrFun  = 'unif'  ;                                                          
%DistrPar = cell(M,1);                                                          
%for i=1:M; DistrPar{i} = [ xmin(i) xmax(i) ] ; end                             
% Name of parameters (will be used to costumize plots):                        
X_labels = X_Labels ;                                                                                          
                                                                               
% Define output:                                                               
myfun = 'FBAlight' ;  

%% Step 3 (sample inputs space)                                                
                                                                               
r = 100 ; % Number of Elementary Effects                                       
% [notice that the final number of model evaluations will be equal to          
% r*(M+1)]                                                                     
                                                                               
% option 1: use the sampling method originally proposed by Morris (1991):      
% L = 6  ; % number of levels in the uniform grid                              
% design_type  = 'trajectory'; % (note used here but required later)           
% X = Morris_sampling(r,xmin,xmax,L); % (r*(M+1),M)                            
                                                                               
% option 2: Latin Hypercube sampling strategy                                  
SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';                                                        
% other options for design type:                                               
%design_type  = 'trajectory';                                                  
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);              
                                                            
%% Step 4 (run the model)                                                      
tic
Y = model_evaluation(myfun,X,model1) ; % size (r*(M+1),1)              
toc
% This takes around 45 minutes (r = 100, M = 67)
% If each model_evaluation takes 0.3 sec

%% Step 5 (Computation of the Elementary effects)                              
                                                                               
% Compute Elementary Effects:                                                  
[ mi, sigma ] = EET_indices(r,xmin,xmax,X,Y,design_type);                      
                                                                               
% Plot results in the plane (mean(EE),std(EE)):   
% mi : tells the relative importance of each parameter on the output
% sigma : s used to detect factors involved in interaction with other factors or whose effect is non-linear
% EET_plot(mi, sigma,X_labels ) 

% Try on log-log axes
figure (1)
loglog(mi, sigma, 'o', 'MarkerSize', 0.01)
xlabel('Mean of EEs','FontSize',14) %,'FontName',fn)
ylabel('Standard deviation of EEs','FontSize',14) %,'FontName',fn)
grid on
a = linspace(1,length(param),length(param))'
b = num2str(a)
c = cellstr(b)   % Why not simply use c = num2cell(a) ?
%dx = 0.1
%dy = 0.1
%text(mi + dx, sigma + dy, c, 'Fontsize', 12);
text(mi, sigma, c, 'Fontsize', 10);
saveas(gcf,'Morris1.eps')

% A closer look
xlim([8E-06 1])
saveas(gcf,'Morrissmall.eps')

% Legend as a table
Legend = [c, X_Labels']
% create a table
Legend = cell2table(Legend,'variableNames', {'Number','Components'})
% write to file
writetable(Legend,'LegendMorris.xlsx')

% Use bootstrapping to derive confidence bounds:                               
Nboot=100;                                                                     
[mi,sigma,EE,mi_sd,sigma_sd,mi_lb,sigma_lb,mi_ub,sigma_ub] = ...               
EET_indices(r,xmin,xmax,X,Y,design_type,Nboot);                                
                                                                               
% Plot bootstrapping results in the plane (mean(EE),std(EE)):                  
EET_plot(mi,sigma,X_labels,mi_lb,mi_ub,sigma_lb,sigma_ub)                      

% Plot with error bars
figure (2)
errorbar(mi, sigma, sigma_lb,sigma_ub, mi_lb,mi_ub, 'o', 'MarkerSize', 2)
xlabel('Mean of EEs','FontSize',14) %,'FontName',fn)
ylabel('Standard deviation of EEs','FontSize',14) %,'FontName',fn)
set(gca,'YScale','log', 'XScale', 'log');
saveas(gcf,'Morris2.eps')

% Repeat computations using a decreasing number of samples so as to assess     
% if convergence was reached within the available dataset: 
% If we see a straight line, then convergence was reached when increasing
% the number of model evaluations...
rr = [ r/5:r/5:r ] ;                                                           
m_r = EET_convergence(EE,rr);                                                  
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_r,rr*(M+1),[],[],[],...                             
'no of model evaluations','mean of EEs',X_labels)                              
                                                                               
% Repeat convergence analysis using bootstrapping:                             
Nboot = 100;                                                                   
rr = [ r/5:r/5:r ] ;                                                           
[m_r,s_r,m_lb_r,m_ub_r] = EET_convergence(EE,rr,Nboot);                        
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],...                     
'no of model evaluations','mean of EEs',X_labels)                              
     
% time the code
toc;

%% Step 6 (Adding up new samples)                                              
%{                                                                               
r2 = 200 ; % increase of base sample size                                      
[X2,Xnew] = OAT_sampling_extend(X,r2,DistrFun,DistrPar,design_type);% extended 
% sample (it includes the already evaluated sample 'X' and the new one)        
                                                                               
% Evaluate model against the new sample                                        
Ynew = model_evaluation(myfun,Xnew,rain,evap,flow) ; % size((r2-r)*(M+1),1)    
                                                                               
% Put new and old results together                                             
Y2=[Y;Ynew]; % size (r2*(M+1),1)                                               
                                                                               
% Recompute indices                                                            
Nboot=100;                                                                     
[mi_n,sigma_n,EEn,mi_sdn,sigma_sdn,mi_lbn,sigma_lbn,mi_ubn,sigma_ubn] = ...      
EET_indices(r2,xmin,xmax,X2,Y2,design_type,Nboot);                             
EET_plot(mi_n,sigma_n,X_labels,mi_lbn,mi_ubn,sigma_lbn,sigma_ubn)                
                                                                               
% Repeat convergence analysis                                                  
Nboot = 100;                                                                   
rr2 = [ r2/5:r2/5:r2 ] ;                                                       
[m_rn,s_rn,m_lb_rn,m_ub_rn] = EET_convergence(EEn,rr2,Nboot);                  
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_rn,rr2*(M+1),m_lb_rn,m_ub_rn,[],...                 
'no of model evaluations','mean of EEs',X_labels)                              
%}
