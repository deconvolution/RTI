## using packages
using RTI,MATLAB,Optim,LineSearches,Plots
## define empty arrays for loading data
# observed data
R_true=Vector{Vector{Float64}}();
# source grid location
s1=Vector{Vector{Int64}}();
s2=Vector{Vector{Int64}}();
s3=Vector{Vector{Int64}}();
# receiver grid location
r1=Vector{Vector{Int64}}();
r2=Vector{Vector{Int64}}();
r3=Vector{Vector{Int64}}();
# source true location
s1t=Vector{Vector{Float64}}();
s2t=Vector{Vector{Float64}}();
s3t=Vector{Vector{Float64}}();
# receiver true location
r1t=Vector{Vector{Float64}}();
r2t=Vector{Vector{Float64}}();
r3t=Vector{Vector{Float64}}();
T0=0;
## read data
tt=readdir("./obs/");
file_name=tt;
for I=1:size(tt,1)
    global R_true,s1,s2,s3,r1,r2,r3;
    tt2=RTI.readmat(string("./obs/",tt[I]),"data");
    R_true=push!(R_true,tt2["Rp"][:,4]);
    s1=push!(s1,round.(Int64,tt2["S"][:,1]));
    s2=push!(s2,round.(Int64,tt2["S"][:,2]));
    s3=push!(s3,round.(Int64,tt2["S"][:,3]));
    r1=push!(r1,round.(Int64,tt2["Rp"][:,1]));
    r2=push!(r2,round.(Int64,tt2["Rp"][:,2]));
    r3=push!(r3,round.(Int64,tt2["Rp"][:,3]));
end
## define parameters for the forward problem
# load initial velocity model
tt=RTI.readmat("./m.mat","data");
# number of grids
nx=round(Int64,tt["nx"]);
ny=round(Int64,tt["ny"]);
nz=round(Int64,tt["nz"]);

# Coordinate
X=tt["X"];
Y=tt["Y"];
Z=tt["Z"];

# grid spacing
h=tt["dx"];

# initial velocity model
v=zeros(nx,ny,nz);
v[:] .=5333;

# reshape the velocity to a 1D array
vc=reshape(v,nx*ny*nz,);
## create path for inversions
p2=@__FILE__;
p3=replace(p2,".jl"=>"");
if isdir(p3)==0
    mkdir(p3);
end
if isdir(string(p3,"/temp"))==0
    mkdir(string(p3,"/temp"))
end
if isdir(string(p3,"/final"))==0
    mkdir(string(p3,"/final"))
end
## define cost function
mutable struct data2
    nx
    ny
    nz
    h
    v
    T
end
data=data2(0,0,0,0,0,0);
function data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3)
    """
    Target: to implement the forward problem.
    Comments: shots are parallelized.
    Output:
      T: travel time for all events.
      chi: L2-norm of the data misfit.
    """
    foreach(rm,filter(endswith(".jld2"),readdir(string(p3,"/temp"),join=true)));
    v=reshape(vc,nx,ny,nz);
    M=[0:40:size(s1,1);size(s1,1)];
    E=zeros(length(s1,));
    for m=1:size(M,1)-1
        Threads.@threads for I=(M[m]+1):M[m+1]
            input_s1=s1[I][1];
            input_s2=s2[I][1];
            input_s3=s3[I][1];

            # Compute forward travel time
            T,R_cal=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            v=v,
            s1=input_s1,
            s2=input_s2,
            s3=input_s3,
            T0=T0,
            r1=r1[I],
            r2=r2[I],
            r3=r3[I]);
            E[I]=.5*sum((R_cal-R_true[I]).^2);
            RTI.JLD2.save(string(p3,"/temp/source_",I,".jld2"), "data",T);
        end
    end
    chi=sum(E);
    return chi
end
## define gradient function, the gradient is computed using adjoint state method
function g!(storage,vc)
    """
    Target: to supply the routine for the update direction of the velocity.
    Comments: shots are parallelized.
    storage: the update direction.
    """
    DV=zeros(nx,ny,nz);
    M=[0:40:size(s1,1);size(s1,1)];
    for m=1:size(M,1)-1
        DV2=zeros(nx,ny,nz,length((M[m]+1):M[m+1]));
        @Threads.threads for I=(M[m]+1):M[m+1]
            temp_t=RTI.JLD2.load(string(p3,"/temp/source_",I,".jld2"))["data"];
            input_s1=s1[I][1];
            input_s2=s2[I][1];
            input_s3=s3[I][1];

            # Compute adjoint travel time
            lambda=RTI.eikonal.acoustic_eikonal_adjoint(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            T=temp_t,
            r1=r1[I],
            r2=r2[I],
            r3=r3[I],
            s1=s1[I][1],
            s2=s2[I][1],
            s3=s3[I][1],
            R_cal=temp_t[CartesianIndex.(r1[I],r2[I],r3[I])],
            R_true=R_true[I],
            N=ones(size(temp_t[CartesianIndex.(r1[I],r2[I],r3[I])]))*(N2));
            if m==1
                tt=0;
            else
                tt=(m-1)*length((M[m-1]+1):M[m]);
            end
            DV2[:,:,:,I-tt]=DV2[:,:,:,I-tt]+lambda ./v .^3;
        end
        DV[:,:,:]=DV[:,:,:]+sum(DV2,dims=4);
    end
    mat"""
    $DV(:,:,:)=imgaussfilt3($DV,$fu);
    """
    storage[:]=DV[:];
end
## Optimization
# define the inciden direction for rays to receivers, 1 for +z and -1 for -z.
N2=-1;

# Test inversion.
fu=6;
data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3);
sca=1;
test_storage=zeros(size(vc));
g!(test_storage,vc);
sca=1/maximum(abs.(test_storage));

# Perform inversion
fu=6;
opt1=optimize(vc->data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3)[1],
g!,vc,LBFGS(m=5,alphaguess=LineSearches.InitialQuadratic(α0=sca*50.0,αmin=sca*10.0),
linesearch=LineSearches.Static()),
Optim.Options(iterations=6,store_trace=true,show_trace=true,
x_tol=0,g_tol=0));
## write final model to vtk
vtkfile=RTI.vtk_grid(string(p3,"/final/final_model"),X,Y,Z);
vtkfile["v"]=reshape(opt1.minimizer,nx,ny,nz);
RTI.vtk_save(vtkfile);
## create CSV for source and receivers
#=
CSV.write(string(p3,"/final/receiver location.csv"),DataFrame([r1t' r2t' r3t'],:auto));
CSV.write(string(p3,"/final/source location.csv"),DataFrame([s1t' s2t' s3t'],:auto));
=#
## write velocity to mat file
mutable struct data3
    X
    Y
    Z
    nx
    ny
    nz
    h
    v
end
data=data3(0,0,0,0,0,0,0,0);

Y,X,Z=RTI.meshgrid((1:ny)*h,(1:nx)*h,(1:nz)*h);
data.X=X;
data.Y=Y;
data.Z=Z;
data.v=reshape(opt1.minimizer,nx,ny,nz);
data.h=h;
data.nx=nx;
data.ny=ny;
data.nz=nz;

file=RTI.matopen(string(p3,"/final/final_model.mat"), "w");
write(file,"data",data);
close(file);
## plot data misfit
trace=opt1.trace;
trace_err=zeros(size(trace));
trace_time=zeros(size(trace));
for i in 1:length(trace)
    trace_err[i]=parse(Float64, split(string(trace[i]))[2]);
    trace_time[i]=parse(Float64, split(string(trace[i]))[end]);
end

ax=plot(trace_time,trace_err,seriestype=:scatter,
xlabel="time [s]",
ylabel="data misfit [s^2]",yscale=:log10);
display(ax)
