This code requires the following softwares and toolboxes:
 - MATLAB
 - The open-source MATLAB toolbox called ‘ Constraint-Based Reconstruction and Analysis (COBRA)Toolbox’
 - A quadratic programming solver such as GUROBI, MOSEK, TOMLAB or IBM ILOG CPLEX
 - The Sensitivity Analysis For Everyone MATLAB Toolbox (SAFE)
 - Python with packages 'networkx' and 'bumpy'.

1. Instruction for downloading the COBRA Toolbox and quadratic solver
All instructions can be found at : https://opencobra.github.io/cobratoolbox/stable/installation.html

2. The SAFE toolbox is freely available for non-commercial purposes at : https://www.safetoolbox.info/register-for-download/

3. The open-source Python software can be downloaded at : https://www.python.org/

4. Paste all .m and .py files in your MATLAB and Python working directory

5. Quality control analysis can be produced by running the file 'Quality_control.m'.

5. Figure 1 and quantitative sensitivity analysis of parameter combinations can be performed with the script 'localsens_script_vf.m', which load the scripts 'model_setup.m' as well as the function 'FBAlocalsens.m' and 'FBAlight.m'.

6. Figure 2 is produced by running the script 'globalSA_Morris.m'

7. Figure 3 and 5 are obtained after running the script 'ReactionDeletionSensitivity.m'

8. Network metrics (Figure 4) in Python (i.e., degree distribution and centrality indices) are obtained by running the script 'bipartite.m' in MATLAB and running the script 'main_network_analysis.py' in Python, which import all functions in the file 'network_metrics.py'.



 



