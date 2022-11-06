"
3D Eikonal equation solver
"
module eikonal
export acoustic_eikonal_forward,acoustic_eikonal_adjoint
##
include("..//utilities/utilities.jl");
## using packages
using Random,MAT,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,
FileIO,JLD2
## acoustic traveltime forward
include("./acoustic_eikonal_forward.jl");
## acoustic traveltime adjoint
include("./acoustic_eikonal_adjoint.jl");
end
