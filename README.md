# licq-path-following

This code is to accompany paper entitled "Sensitivity-Based Economic NMPC with a Path-Following Approach". 
The paper appears at Processes journal (http://www.mdpi.com/journal/processes/special_issues/real_time_optimization).

The code requires:
- MATLAB version 2014b or above
- CasADi version 3.1 or above
- TOMLAB optimization package (for QP solver)

Before running the code, you need to edit startup.m file at the root directory.
Make sure you put your CasADi installation folder in the startup.m. Additionally, you need to run startup file provided by TOMLAB too.
Failing to run those startups file, the code will not be able to recognize CasADi and TOMLAB commands.

The code is an implementation of economic MPC scheme. workflow to run the code include:
1. steady-state optimization (yielding set-point)
2. dyanamic optimization (MPC controller)

The steady-state code is located at \models\disColA\ folder, named distACstrSS.m. Run this file to get set-point for MPC controller.

There are two type of MPC controllers (located at folder \nmpc, namely:
1. ideal NMPC (iNMPC) --> file iNMPC.m
2. path-following (pf-NMPC) --> file pfNMPC.m

The main file for code is NMPCDistCstr.m (located at \models\disColA). Here, one can choose to use either iNMPC of pfNMPC controller.


Enjoy!
