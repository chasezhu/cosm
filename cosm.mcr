# Example for point mutation analysis.
# This is a qualitative approximation and starting point for...
# work, a fully automatic free energy predictor is under...
# see also the user manual at Recipes > Calculate free energies.
Clear
# Load the structure
LoadSce '192_W_R_mined.sce'
#NiceOri
#DelRes Hetgroup
# Set the wild type residue to mutate, here Arg 98 in Molecule A
wtresnamelist='Tyr','Trp','Leu','Phe','Ala'
wtresnumberlist='20','58','57','86','230'
wtmolist='Mol A','Mol A','Mol A','Mol B','Mol A'
wtresnummol=wtresnumberlist+' '+wtmolist
wtres=wtresnamelist+' '+wtresnumberlist+' '+wtmolist
ligname='PLT Mol A'
# List the mutants, wt must appear again as the first one
#alascan mutlist='(wtresname)','Ala'
#all 20  mutlist='(wtresname)','Arg','Lys','His','Asp','Glu','Asn','Gln','Ser','Thr','Cys','Gly','Ala','Val','Leu','Ile','Pro','Met','Tyr','Phe','Trp'
#11 mix  mutlist='(wtresname)','Lys','His','Asp','Gln','Thr','Ala','Leu','Pro','Met','Phe'
# actual mutation list
mutlist1='(wtresnamelist(1))','Leu'
mutlist2='(wtresnamelist(2))','Gly','Ala','Phe','Tyr','Met'
mutlist3='(wtresnamelist(3))','Thr','Tyr'
mutlist4='(wtresnamelist(4))','Leu'
mutlist5='(wtresnamelist(5))','Pro','Thr'
# Count the mutants per position
Console On
for i=1 to 5
  mutants(i)=count mutlist(i)
  wtres(i)=wtresnamelist(i)+' '+wtresnumberlist(i)+' '+wtmolist(i)
  wtresnummol(i)=wtresnumberlist(i)+' '+wtmolist(i)
# Calculate total number of runs
runs=mutants1*mutants2*mutants3*mutants4*mutants5
print 'Number of runs (runs)'
# 0.65 kJ/Mol is a guesstimate of the entropic cost of exposing...
# We divide by 6.02214e20 to obtain J and multiply with JToUnit
# to get the currently selected energy unit
surfcost=(0.65e0/6.02214e20)*JToUnit
# Prepare for simulation, add missing hydrogens etc.
Clean
Style Ribbon,Stick
# Add a simulation cell for energy minimization
Cell Auto,Extension=10
# Set the force field parameters
ForceField NOVA2,SetPar=Yes
# Create and minimize water shell, if force field is Amber/Yamber
#Experiment Neutralization
#  NaCl 0
#  Speed fast
#Experiment On
#Wait ExpEnd
#DelRes Water with distance>6 from !Water
# Energy minimize the protein
#Experiment Minimization
#Experiment On
#Wait ExpEnd
# Save the wildtype for reuse later
#SaveSce wildtype
# Ensure that the surface resolution is high enough for...
SurfPar Resolution=3,Molecular=Numeric
# Try all mutations

#for i=1 to count mutlist
# Loop over the runs
for i=0 to runs-1
  # Get the list of mutants for this run
  runmutlist=mutlist1(1+i%mutants1),mutlist2(1+(i//mutants1)%mutants2),mutlist3(1+i//(mutants1*mutants2)%mutants3),mutlist4(1+i//(mutants1*mutants2*mutants3)%mutants4),mutlist5(1+i//(mutants1*mutants2*mutants3*mutants4)%mutants5)
  print 'Run (i): (runmutlist)'
  Console Off
  epotiso(i) = 0
  esoliso(i)= 0
  for j=1 to 5
    mutresname=runmutlist(j)
    ShowMessage 'Analyzing mutation (wtres(j)) -> (mutresname)'
    # First calculate the energy of the isolated/unfolded amino...
    # (include acetyl and N-methyl caps) so that we can detect
    # the energetic cost of burying a polar/charged residue
    Clear
    BuildRes (mutresname)
    AddCap ACE+NME
    Clean
    Cell Auto,Extension=10
    # Energy minimize the isolated/unfolded amino acid. Since we...
    # have a solvent shell, only do a quick steepest descent...
    # Test WITHOUT MINIMIZATION
    Experiment Minimization
      Convergence immediate
    Experiment On
    Wait ExpEnd
    # Ignore the net charge when calculating the energies
    ChargeObj all,0
    energyref=Energy
    epotiso(i)=epotiso(i)+energyref
    esolcoulomb,esolvdw = SolvEnergy
    molsurf = Surf molecular
    esoliso(i)=epotiso(i)+esolcoulomb+esolvdw+molsurf*surfcost
  # Load the protein again
  LoadSce '192_W_R_mined.sce'
  Fixres all
  #ZoomAtom CA Res (wtres)
  # Remove the water shell temporarily
  #RemoveObj Water
  # Fix everything not close by
  for j=1 to 5
    FreeAtom all with distance <10 from Res (wtres(j))
    AddSpring Atom NZ res Lys 287 Mol A,Atom C4A res PLT 1001 Mol A,Len=3.2,SFC=10
    FreeRes PLT  Mol A
    mutresname=runmutlist(j) 
    # Make the mutation
    SwapRes (wtres(j)),(mutresname)
    if Structure and mutresname!='Gly' and mutresname!='Ala'
      # Additionally optimize the side-chain rotamer in YASARA...
      OptimizeRes (mutresname) (wtresnummol(j)),Method=SCWALL
    # Optimization changed the cell boundaries
    Boundary periodic
  # Add the water shell again
  #AddObj Water
  # Delete the waters that bump into our mutated residue
  DelRes Water with distance<3 from (mutresname) (wtresnummol(j))
  # Energy minimize the environment of the residue
  # Test WITHOUT MINIMIZATION
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  RemoveObj Water
  # Calculate potential and solvation energy
  ChargeObj all,0
  epot(i) = Energy
  esolcoulomb,esolvdw = SolvEnergy
  molsurf = Surf molecular
  esol(i)=esolcoulomb+esolvdw+molsurf*surfcost
  e_lig(i)= Energyres (ligname)
  # Save the scene for later manual inspection, if ligand energy is stabilizing
  if (e_lig(i)-e_lig(0))<=0 
    SaveYob obj 1, '192_W_R_mined_(i)'
  # Print the results
Console Off 
Print 'Results of the simple mutation analysis, negative "Difference to wt" are better than the wildtype:'
Print '[This is a qualitative estimate, quantitative free energies will be supported in the future]'
for i=0 to runs-1
  # Calculate folding energy, folded minus isolated/unfolded...
  runmutlist=mutlist1(1+i%mutants1),mutlist2(1+(i//mutants1)%mutants2),mutlist3(1+i//(mutants1*mutants2)%mutants3),mutlist4(1+i//(mutants1*mutants2*mutants3)%mutants4),mutlist5(1+i//(mutants1*mutants2*mutants3*mutants4)%mutants5)
  efold(i)=epot(i)+esol(i)-epotiso(i)-esoliso(i)
  result='Mutation (i) -> (runmutlist) '+ ' Folded protein: (0.00+epot(i)+esol(i)) (EnergyUnit) '+ '[Potential energy = (0.00+epot(i)) , Solvation energy = (0.00+esol(i)) ]'+ ' Unfolded residue: (0.00+epotiso(i)+esoliso(i)) (EnergyUnit) '+ '[Potential energy = (0.00+epotiso(i)) , Solvation energy = (0.00+esoliso(i)) ]'+ ' Folding energy: (0.00+efold(i)) (EnergyUnit)'+ ' Difference to wt: (0.00+efold(i)-efold(0)) (EnergyUnit) + Ligand Stabilization: (0.00+e_lig(i)-e_lig(0)) (EnergyUnit) '
  Print '(result)'
  Tabulate (i)
  #mutation(i)='(wtresnummol)'+' (mutlist(i))'
  Tabulate (0.00+efold(i)-efold(0))
  Tabulate (0.00+e_lig(i)-e_lig(0))
  Tabulate '(runmutlist)'
header=" run number  Protein Stab E  Ligand Stab E     WT Res Mut"   
SaveTab default,'combi_SSM_lig_analysis',Format=Text,Columns=4, NumFormat=14.3f,Header=(header)
Console On