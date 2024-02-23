#Using parflow

# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*

set proj "/home/rappl/PROJETS/Wank1D_ex"
set forc "${proj}/forcings/forc_MT"
set simus "${proj}/simus"
set tcl_pfidb "${proj}/tcl_pfidb"
cd $simus

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
pfset FileVersion 4

#-----------------------------------------------------------------------------
# Process Topology
#-----------------------------------------------------------------------------
pfset Process.Topology.P        1
pfset Process.Topology.Q        1
pfset Process.Topology.R        1



#-----------------------------------------------------------------------------
# Results output
#-----------------------------------------------------------------------------

set simu "test_LRH"
file delete -force $simu
file mkdir $simu
cd $simu

#-----------------------------------------------------------------------------
# Prepare input files
#-----------------------------------------------------------------------------

file copy ${tcl_pfidb}/Wank1D_ex.tcl ./Wank1D_ex.tcl

###############################################################################

file copy ${forc}/forc.txt ./forcagePF.txt.0
file copy ${forc}/lai.dat ./lai.dat
file copy ${forc}/drv_vegp.dat ./drv_vegp.dat
file copy ${forc}/sai.dat ./sai.dat
file copy ${forc}/z0m.dat ./z0m.dat
file copy ${forc}/displa.dat  ./displa.dat
file copy ${forc}/veg_map.pfb ./veg_map.pfb
file copy ${forc}/drv_clmin.dat ./drv_clmin.dat
file copy ${forc}/drv_vegm.dat ./drv_vegm.dat
pfdist veg_map.pfb   
#-----------------------------------------------------------------------------
# Computational Grid
#-----------------------------------------------------------------------------

pfset ComputationalGrid.Lower.X                 0.0
pfset ComputationalGrid.Lower.Y                	0.0
pfset ComputationalGrid.Lower.Z                	0.0

pfset ComputationalGrid.DX	               	10.0
pfset ComputationalGrid.DY	                10.0
pfset ComputationalGrid.DZ	               	1.0

pfset ComputationalGrid.NX	     	       	1
pfset ComputationalGrid.NY	     		1
pfset ComputationalGrid.NZ                  	1

pfset ComputationalGrid.NZ                  	30
#-----------------------------------------------------------------------------
# The Names of the GeomInputs
#-----------------------------------------------------------------------------
pfset GeomInput.Names "domain_input H1_input H2_input H3_input H4_input H5_input"

#-----------------------------------------------------------------------------
# Domain Geometry Input
#-----------------------------------------------------------------------------
# 3 méthodes possibles (pfsol, indicator files et les confifurations normales comme ici) 
pfset GeomInput.domain_input.InputType           Box
pfset GeomInput.domain_input.GeomName            domain

pfset Geom.domain.Lower.X                        0
pfset Geom.domain.Lower.Y                        0
pfset Geom.domain.Lower.Z                        0

pfset Geom.domain.Upper.X                        10
pfset Geom.domain.Upper.Y                        10
pfset Geom.domain.Upper.Z                        30

# les patches sont les bords de la boite
pfset Geom.domain.Patches                        "x-lower x-upper y-lower y-upper z-lower z-upper"
#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
# c'est la geometrie de l'ensemble du domaine
pfset Domain.GeomName                            domain

#-----------------------------------------------------------------------------
# H1 Geometry Input
#-----------------------------------------------------------------------------
# Ici pour chaque horizon on definit une boite ( on a 3 horizon)

# the vertical coordinate origin is at the bottom
# horizon H1 surface 4cm permet de representer la croute des sol du Niger
pfset GeomInput.H1_input.InputType               Box
pfset GeomInput.H1_input.GeomName                H1

pfset Geom.H1.Lower.X                            0
pfset Geom.H1.Lower.Y                            0
pfset Geom.H1.Lower.Z                            25
#corresponds to
#pfset Geom.H1.Lower.Z                            23.96

pfset Geom.H1.Upper.X                            10
pfset Geom.H1.Upper.Y                            10
pfset Geom.H1.Upper.Z                            30

#-----------------------------------------------------------------------------
# H2 Geometry Input
#-----------------------------------------------------------------------------
# horizon H2 surface ==> 50cm zone racinaire des herbaceees
pfset GeomInput.H2_input.InputType               Box
pfset GeomInput.H2_input.GeomName                H2

pfset Geom.H2.Lower.X                            0
pfset Geom.H2.Lower.Y                            0
#corresponds to
#pfset Geom.H2.Lower.Z                            23.5
pfset Geom.H2.Lower.Z                            19

pfset Geom.H2.Upper.X                            10
pfset Geom.H2.Upper.Y                            10
#corresponds to
#pfset Geom.H2.Upper.Z                            23.96 
pfset Geom.H2.Upper.Z                            25

#-----------------------------------------------------------------------------
# H3 Geometry Input
#-----------------------------------------------------------------------------
# horizon H3  ==> 130 permet de representer le gradient de proprietes dans la ZNS 
pfset GeomInput.H3_input.InputType               Box
pfset GeomInput.H3_input.GeomName                H3

pfset Geom.H3.Lower.X                            0
pfset Geom.H3.Lower.Y                            0
#corresponds to
#pfset Geom.H3.Lower.Z                            22.7
pfset Geom.H3.Lower.Z                            13

pfset Geom.H3.Upper.X                            10
pfset Geom.H3.Upper.Y                            10
#corresponds to
#pfset Geom.H3.Upper.Z                            23.5
pfset Geom.H3.Upper.Z                            19

#-----------------------------------------------------------------------------
# H4 Geometry Input
#-----------------------------------------------------------------------------
# horizon H4  ==> 250 Zone argileuse des sols au benin
pfset GeomInput.H4_input.InputType               Box
pfset GeomInput.H4_input.GeomName                H4

pfset Geom.H4.Lower.X                            0
pfset Geom.H4.Lower.Y                            0
#corresponds to
#pfset Geom.H4.Lower.Z                            21.5
pfset Geom.H4.Lower.Z                            7

pfset Geom.H4.Upper.X                            10
pfset Geom.H4.Upper.Y                            10
#corresponds to
#pfset Geom.H4.Upper.Z                            22.7
pfset Geom.H4.Upper.Z                            13

#-----------------------------------------------------------------------------
# H5 Geometry Input
#-----------------------------------------------------------------------------
# horizon H5  ==> 24 eppaisseur du Regolhite
pfset GeomInput.H5_input.InputType               Box
pfset GeomInput.H5_input.GeomName                H5

pfset Geom.H5.Lower.X                            0
pfset Geom.H5.Lower.Y                            0
pfset Geom.H5.Lower.Z                            0

pfset Geom.H5.Upper.X                            10
pfset Geom.H5.Upper.Y                            10
#corresponds to
#pfset Geom.H5.Upper.Z                            21.5
pfset Geom.H5.Upper.Z                            7

#-----------------------------------------------------------------------------
# Perm = Ksat (m/h)
# On associe une valeur à chaque horizon
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                           "H1 H2 H3 H4 H5"

##### Millet
pfset Geom.H1.Perm.Type                         Constant
pfset Geom.H1.Perm.Value                        0.0012

pfset Geom.H2.Perm.Type                         Constant
pfset Geom.H2.Perm.Value                        0.18

pfset Geom.H3.Perm.Type                         Constant
pfset Geom.H3.Perm.Value                        0.18

pfset Geom.H4.Perm.Type                         Constant
pfset Geom.H4.Perm.Value                        0.18

pfset Geom.H5.Perm.Type                         Constant
pfset Geom.H5.Perm.Value                        0.18




pfset Perm.TensorType                           TensorByGeom
pfset Geom.Perm.TensorByGeom.Names              "domain"

# pas d'anisotropie ( tensorVal =1)
pfset Geom.domain.Perm.TensorValX               1.0
pfset Geom.domain.Perm.TensorValY               1.0
pfset Geom.domain.Perm.TensorValZ               1.0


#-----------------------------------------------------------------------------
# Specific Storage (/m)
#-----------------------------------------------------------------------------
pfset SpecificStorage.Type                      Constant
pfset SpecificStorage.GeomNames                 "domain"

pfset Geom.domain.SpecificStorage.Value         1.0e-3

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------
pfset Phase.Names                               "water"

pfset Phase.water.Density.Type	                Constant
pfset Phase.water.Density.Value	                1.0

pfset Phase.water.Viscosity.Type	        Constant
pfset Phase.water.Viscosity.Value	        1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
pfset Contaminants.Names		        ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------
pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
# Le pas de temps (ici 1/2h)
#-----------------------------------------------------------------------------
pfset TimingInfo.BaseUnit        0.5
pfset TimingInfo.StartCount      0.0
pfset TimingInfo.StartTime       0.0

pfset TimingInfo.DumpInterval    24
pfset TimeStep.Type              Constant
pfset TimeStep.Value             0.5
#2days
#pfset TimingInfo.StopTime	48
#1 year  
pfset TimingInfo.StopTime	 8760
#2 year
#pfset TimingInfo.StopTime	17520 
#3 years
#pfset TimingInfo.StopTime	26280
#4 year
#pfset TimingInfo.StopTime	35064
#4 year  + 2006 at beginning
#pfset TimingInfo.StopTime       43824
#pfset TimingInfo.StopTime	480
#6 year
#pfset TimingInfo.StopTime       52584
#11yrs
#ipfset TimingInfo.StopTime	96360
# 10 ans
#pfset TimingInfo.StopTime       87672

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

#pfset Geom.Porosity.GeomNames      "domain H1"
pfset Geom.Porosity.GeomNames      "H1 H2 H3 H4 H5"

pfset Geom.H1.Porosity.Type             Constant
pfset Geom.H1.Porosity.Value            0.358

pfset Geom.H2.Porosity.Type             Constant
pfset Geom.H2.Porosity.Value   	        0.358

pfset Geom.H3.Porosity.Type             Constant
pfset Geom.H3.Porosity.Value   	        0.321

pfset Geom.H4.Porosity.Type             Constant
pfset Geom.H4.Porosity.Value   	        0.358

pfset Geom.H5.Porosity.Type             Constant
pfset Geom.H5.Porosity.Value   	        0.34

#-----------------------------------------------------------------------------
# Mobility
#-----------------------------------------------------------------------------
pfset Phase.water.Mobility.Type        Constant
pfset Phase.water.Mobility.Value       1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                        ""

#-----------------------------------------------------------------------------2005
# Time Cycles
#-----------------------------------------------------------------------------

#pfset Cycle.Names                       "constant rainrec"
pfset Cycle.Names                       "constant"

pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat              -1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type                   DirEquilRefPatch
pfset Patch.z-lower.BCPressure.RefPatch               "z-lower"
pfset Patch.z-lower.BCPressure.RefGeom               "domain"
pfset Patch.z-lower.BCPressure.Cycle                  "constant"
pfset Patch.z-lower.BCPressure.alltime.Value          7.0

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

pfset Patch.z-upper.BCPressure.Type                   OverlandFlow
pfset Patch.z-upper.BCPressure.Cycle                  "constant"
pfset Patch.z-upper.BCPressure.alltime.Value          0.0


#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
pfset TopoSlopesX.Type 				"Constant"
pfset TopoSlopesX.GeomNames                     "domain"
pfset TopoSlopesX.Geom.domain.Value 		0.00

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
pfset TopoSlopesY.Type 				"Constant"
pfset TopoSlopesY.GeomNames                     "domain"
#Bele:
pfset TopoSlopesY.Geom.domain.Value 		0.003

#---------------------------------------------------------
# Mannings coefficient (min^1/3/m)
#---------------------------------------------------------
pfset Mannings.Type 				"Constant"
pfset Mannings.GeomNames 			"domain"
#hr
pfset Mannings.Geom.domain.Value 		0.0000056

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "H1 H2 H3 H4 H5"


## Millet

pfset Geom.H1.RelPerm.N               1.77
pfset Geom.H1.RelPerm.Alpha           1.15

pfset Geom.H2.RelPerm.Alpha           3.33
pfset Geom.H2.RelPerm.N               2

pfset Geom.H3.RelPerm.Alpha            3.33
pfset Geom.H3.RelPerm.N                2

pfset Geom.H4.RelPerm.N                5
pfset Geom.H4.RelPerm.Alpha            2

pfset Geom.H5.RelPerm.Alpha            5
pfset Geom.H5.RelPerm.N                2


#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type            VanGenuchten
#pfset Phase.Saturation.GeomNames       "domain H1"
pfset Phase.Saturation.GeomNames       "H1 H2 H3 H4 H5"

## Millet 

pfset Geom.H1.Saturation.Alpha         1.15
pfset Geom.H1.Saturation.N             1.77
pfset Geom.H1.Saturation.SRes          0.028
pfset Geom.H1.Saturation.SSat          0.9

pfset Geom.H2.Saturation.Alpha        3.33
pfset Geom.H2.Saturation.N            2
pfset Geom.H2.Saturation.SRes          0.064
pfset Geom.H2.Saturation.SSat          0.9

pfset Geom.H3.Saturation.Alpha        3.33
pfset Geom.H3.Saturation.N            2
pfset Geom.H3.Saturation.SRes          0.13
pfset Geom.H3.Saturation.SSat          0.9

pfset Geom.H4.Saturation.Alpha         5
pfset Geom.H4.Saturation.N             2
pfset Geom.H4.Saturation.SRes          0.13
pfset Geom.H4.Saturation.SSat          0.9

pfset Geom.H5.Saturation.Alpha         5
pfset Geom.H5.Saturation.N             2
pfset Geom.H5.Saturation.SRes          0.16
pfset Geom.H5.Saturation.SSat          0.9

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------
pfset PhaseSources.water.Type                         "Constant"
#pfset PhaseSources.water.GeomNames                    "domain"
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value            0.0


#----------------------------------------------------------------
# CLM Settings:
# ---------------------------------------------------------------

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
pfset ICPressure.Type                                 HydroStaticPatch
pfset ICPressure.GeomNames                            "domain"

#pfset Geom.domain.ICPressure.RefGeom                  "domain"
pfset Geom.domain.ICPressure.RefGeom                  domain

pfset Geom.domain.ICPressure.Value                    -10
pfset Geom.domain.ICPressure.RefPatch                 z-upper

#-----------------------------------------------------------------------------
# Set Outputs
#-----------------------------------------------------------------------------
pfset Solver.PrintDZMultiplier				True
pfset Solver.PrintOverlandSum				True
pfset Solver.PrintSlopes				True
pfset Solver.PrintEvapTrans				True
pfset Solver.PrintEvapTransSum				True
pfset Solver.PrintSubsurf                               True
pfset Solver.PrintMannings                              True
pfset Solver.PrintPressure                              True
pfset Solver.PrintSaturation                            True
pfset Solver.PrintMask                                  True
pfset Solver.PrintSpecificStorage 		        True
pfset Solver.PrintSlopes 		        	True


#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------
pfset KnownSolution                                      NoKnownSolution

#-----------------------------------------------------------------------------
# Variable DZ
#-----------------------------------------------------------------------------
pfset Solver.TerrainFollowingGrid                        True
#pfset Solver.Nonlinear.VariableDz                     False

pfset Solver.Nonlinear.VariableDz                        True

pfset dzScale.GeomNames				        domain

pfset dzScale.GeomNames                                  domain
pfset dzScale.Type                                       nzList

#0 is bottom layer
pfset dzScale.nzListNumber                               30

# H5 horizon : aquifer 1.2 - 30m
pfset Cell.0.dzScale.Value                               20.0
pfset Cell.1.dzScale.Value                               3.0
pfset Cell.2.dzScale.Value                               3.0
pfset Cell.3.dzScale.Value                               1.0
pfset Cell.4.dzScale.Value                               0.7
pfset Cell.5.dzScale.Value                               0.3
pfset Cell.6.dzScale.Value                               0.7
#H4 horizon : potential clay layer 0.7 - 1.2m
pfset Cell.7.dzScale.Value                               0.083
pfset Cell.8.dzScale.Value                               0.083
pfset Cell.9.dzScale.Value                               0.083
pfset Cell.10.dzScale.Value                              0.083
pfset Cell.11.dzScale.Value                              0.083
pfset Cell.12.dzScale.Value                              0.083
#H3 horizon 0.2 - 0.7
pfset Cell.13.dzScale.Value                              0.083
pfset Cell.14.dzScale.Value                              0.083
pfset Cell.15.dzScale.Value                              0.083
pfset Cell.16.dzScale.Value                              0.083
pfset Cell.17.dzScale.Value                              0.083
pfset Cell.18.dzScale.Value                              0.083
#H2 horizon (herbaceous root zone) 0.01 - 0.2m
pfset Cell.19.dzScale.Value                              0.031
pfset Cell.20.dzScale.Value                              0.031
pfset Cell.21.dzScale.Value                              0.031
pfset Cell.22.dzScale.Value                              0.031
pfset Cell.23.dzScale.Value                              0.031
pfset Cell.24.dzScale.Value                              0.031
#H1 horizon (crust) 0 - 0.01m
pfset Cell.25.dzScale.Value                              0.002
pfset Cell.26.dzScale.Value                              0.002
pfset Cell.27.dzScale.Value                              0.002
pfset Cell.28.dzScale.Value                              0.002
pfset Cell.29.dzScale.Value                              0.002

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
 
pfset Solver                                             Richards

pfset Solver.MaxConvergenceFailures			 6
pfset Solver.MaxIter                                     1000000
pfset Solver.AbsTol                                      1E-8
pfset Solver.Drop                                        1E-20
pfset Solver.Nonlinear.MaxIter                           200 
pfset Solver.Nonlinear.ResidualTol                       1e-9
pfset Solver.Nonlinear.StepTol                           1e-30
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-8
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      200
pfset Solver.Linear.MaxRestart                           2
pfset Solver.Linear.Preconditioner                       PFMG

#dis line below new
pfset Solver.Linear.Preconditioner.PCMatrixType     FullJacobian


pfset OverlandFlowSpinUp				 0



#-----------------------------------------------------------------------------
# Set CLM parameters
#-----------------------------------------------------------------------------
 
pfset Solver.LSM                                         CLM



pfset Solver.CLM.MetForcing                              1D
pfset Solver.CLM.MetFileName                             "forcagePF.txt.0"
pfset Solver.CLM.MetFilePath                             ./

pfset Solver.CLM.CLMDumpInterval                         1

pfset Solver.CLM.ForceVegetation			 True

pfset Solver.CLM.RootZoneNZ				 25
pfset Solver.CLM.BinaryOutDir				 False

pfset Solver.CLM.SingleFile				 True

pfset Solver.PrintCLM                                    True
pfset Solver.WriteCLMBinary				 False

pfset Solver.CLM.WriteLogs				 True
pfset Solver.CLM.WriteLastRST				 True

pfset Solver.CLM.EvapBeta				 "none"
pfset Solver.CLM.ResSat					 0.028
pfset Solver.CLM.VegWaterStress                          "Pressure"
pfset Solver.CLM.WiltingPoint                           -150
pfset Solver.CLM.FieldCapacity                          0.0
pfset Solver.CLM.SoiLayer                                25  

/home/rappl/PFTree/Lucie/parflow-dev_MT/

#-----------------------------------------------------------------------------
# Run and Unload the Parflow output files
#-----------------------------------------------------------------------------
 
pfwritedb /home/rappl/PROJETS/Wank1D_ex/tcl_pfidb/Wank1D_ex
puts {before run OK}
set time_deb [clock seconds]
puts [time {
pfrun Wank1D_ex
}]
set time_fin [clock seconds]
puts [expr $time_fin - $time_deb] 
puts {second}

cd ..

