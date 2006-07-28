global codonFactor=0.2068662102090525;
global S_2=0.1155503670328316;
S_2:>1e-07;
global S_1=0.2503723792999368;
S_1:>1e-07;
global S_0=0.4304272712957556;
S_0:>1e-07;
global R_0=0.02278861620075845;
R_0:<1;
global P_2=1;
P_2:<1;
global R_2=2;
R_2:>1;
global P_1=0.7993122949640235;
P_1:<1;
global AC=0.8889337399658167;
global AT=0.6551486104035469;
global CT=1.351450480894584;
global GT=0.6638807095900224;
global CG=0.5717979133032365;
global NS_0:=R_0*S_0;
global NS_1:=R_1*S_1;
global NS_2:=R_2*S_2;
global c_scale:=S_0*P_1+S_1*(1-P_1)*P_2+S_2*(1-P_1)*(1-P_2);
global R_1:=1;
tree_0.Dsim.synRate:=0.0180166/codonFactor;
tree_0.Dsec.synRate:=0.0535959/codonFactor;
tree_0.Node5.synRate:=1e-30/codonFactor;
tree_0.Dmel.synRate:=0.0386392/codonFactor;
tree_0.Node4.synRate:=0.00675294/codonFactor;
tree_0.Dyak.synRate:=0.00462499/codonFactor;
tree_0.Dere.synRate:=0.0258819/codonFactor;
tree_0.Node9.synRate:=0.00346849/codonFactor;
tree_0.Node3.synRate:=0.111341/codonFactor;
tree_0.Dana.synRate:=0.0627428/codonFactor;
tree_0.Node2.synRate:=0.0696995/codonFactor;
tree_0.Dpse.synRate:=0.0662834/codonFactor;
tree_0.Node1.synRate:=0.148251/codonFactor;
tree_0.Dvir.synRate:=0.102542/codonFactor;
tree_0.Dmoj.synRate:=0.0480886/codonFactor;
tree_0.Node14.synRate:=0.0165478/codonFactor;
tree_0.Dgri.synRate:=0.0778238/codonFactor;
tree_1.Dsim.synRate:=0.0180166/codonFactor;
tree_1.Dsec.synRate:=0.0535959/codonFactor;
tree_1.Node5.synRate:=1e-30/codonFactor;
tree_1.Dmel.synRate:=0.0386392/codonFactor;
tree_1.Node4.synRate:=0.00675294/codonFactor;
tree_1.Dyak.synRate:=0.00462499/codonFactor;
tree_1.Dere.synRate:=0.0258819/codonFactor;
tree_1.Node9.synRate:=0.00346849/codonFactor;
tree_1.Node3.synRate:=0.111341/codonFactor;
tree_1.Dana.synRate:=0.0627428/codonFactor;
tree_1.Node2.synRate:=0.0696995/codonFactor;
tree_1.Dpse.synRate:=0.0662834/codonFactor;
tree_1.Node1.synRate:=0.148251/codonFactor;
tree_1.Dvir.synRate:=0.102542/codonFactor;
tree_1.Dmoj.synRate:=0.0480886/codonFactor;
tree_1.Node14.synRate:=0.0165478/codonFactor;
tree_1.Dgri.synRate:=0.0778238/codonFactor;
tree_2.Dsim.synRate:=0.0180166/codonFactor;
tree_2.Dsec.synRate:=0.0535959/codonFactor;
tree_2.Node5.synRate:=1e-30/codonFactor;
tree_2.Dmel.synRate:=0.0386392/codonFactor;
tree_2.Node4.synRate:=0.00675294/codonFactor;
tree_2.Dyak.synRate:=0.00462499/codonFactor;
tree_2.Dere.synRate:=0.0258819/codonFactor;
tree_2.Node9.synRate:=0.00346849/codonFactor;
tree_2.Node3.synRate:=0.111341/codonFactor;
tree_2.Dana.synRate:=0.0627428/codonFactor;
tree_2.Node2.synRate:=0.0696995/codonFactor;
tree_2.Dpse.synRate:=0.0662834/codonFactor;
tree_2.Node1.synRate:=0.148251/codonFactor;
tree_2.Dvir.synRate:=0.102542/codonFactor;
tree_2.Dmoj.synRate:=0.0480886/codonFactor;
tree_2.Node14.synRate:=0.0165478/codonFactor;
tree_2.Dgri.synRate:=0.0778238/codonFactor;
Tree tree_0=((((((Dsim:0.0148487,Dsec:0.0441719)Node5:8.24166e-31,Dmel:0.0318451)Node4:0.00556555,(Dyak:0.00381176,Dere:0.021331)Node9:0.00285861)Node3:0.0917634,Dana:0.0517105)Node2:0.057444,Dpse:0.0546285)Node1:0.122184,(Dvir:0.0845117,Dmoj:0.039633)Node14:0.0136382,Dgri:0.0641397);
Tree tree_1=((((((Dsim:0.0327464,Dsec:0.097414)Node5:1.81757e-30,Dmel:0.0702293)Node4:0.0122739,(Dyak:0.00840622,Dere:0.0470421)Node9:0.00630421)Node3:0.202369,Dana:0.114039)Node2:0.126683,Dpse:0.120474)Node1:0.269456,(Dvir:0.186377,Dmoj:0.0874041)Node14:0.0300768,Dgri:0.14145);
Tree tree_2=((((((Dsim:0.0264991,Dsec:0.0788296)Node5:1.47081e-30,Dmel:0.0568312)Node4:0.00993233,(Dyak:0.0068025,Dere:0.0380675)Node9:0.0051015)Node3:0.163762,Dana:0.0922831)Node2:0.102515,Dpse:0.0974906)Node1:0.21805,(Dvir:0.15082,Dmoj:0.0707294)Node14:0.0243388,Dgri:0.114464);


# Model fit summary

# Log likelihood:    -1209.77722
# Parameters    :             37
# AIC           :     2493.55444
# c-AIC         :     2524.45554
# Class 1
# 	dS    =      1.092
# 	dN    =      0.025
# 	dN/dS =      0.023
# 	Prob  =      0.799
# Class 2
# 	dS    =      0.635
# 	dN    =      0.635
# 	dN/dS =      1.000
# 	Prob  =      0.201
# Class 3
# 	dS    =      0.293
# 	dN    =      0.586
# 	dN/dS =      2.000
# 	Prob  =      0.000


