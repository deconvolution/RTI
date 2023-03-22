## Compute model mosfit
using RTI
## gc
#=
tt=RTI.readmat("./inversion_process_gd/final/final_model.mat","data");
nx=tt["nx"];
ny=tt["ny"];
nz=tt["nz"];
v0=ones(nx,ny,nz)*1000.0;
v_gc=tt["v"];
=#
## lbfgs
tt=RTI.readmat("./inversion_process_lbfgs/final/final_model.mat","data");
v_lbfgs=tt["v"];
## vt
tt=RTI.readmat("./toy_model/true_model.mat","data");
nx=tt["nx"];
ny=tt["ny"];
nz=tt["nz"];
vt=tt["v"];
## v0
#tt=RTI.readmat("../tutorial_inversion_toy_model_strong_contrast/inversion_process_lbfgs/final/final_model.mat","data");
#v0=tt["v"];
v0=copy(vt);
v0[:] .=1000;
##
# misfit initial model
mis0=.5*sum((v0-vt).^2)

# misfit gc
#mis_gc=.5*sum((v_gc-vt).^2)

# misfit lbfgs
mis_lbfgs=.5*sum((v_lbfgs-vt).^2)
