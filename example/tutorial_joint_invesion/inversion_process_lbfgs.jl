## using packages
using RTI,MATLAB,Optim,LineSearches,Plots,Random
## create path for inversions
p2=@__FILE__;
p3=replace(p2,".jl"=>"");
if isdir(p3)==0
    mkdir(p3);
end
if isdir(string(p3,"/temp"))==0
    mkdir(string(p3,"/temp"))
end
if isdir(string(p3,"/temp2"))==0
    mkdir(string(p3,"/temp2"))
end
if isdir(string(p3,"/final"))==0
    mkdir(string(p3,"/final"))
end
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
## give perturbation to s
Random.seed!(1);
ds1=rand(-5:5,length(s1),);
ds2=rand(-5:5,length(s2),);
ds3=rand(-5:5,length(s3),);
for I=1:length(s1)
    s1[I]=s1[I]+[ds1[I]];
    s2[I]=s2[I]+[ds2[I]];
    s3[I]=s3[I]+[ds3[I]];
end
## define parameters for the forward problem
# load initial velocity model
tt=RTI.readmat("./toy_model/true_model.mat","data");
# number of grids
nx=round(Int64,tt["nx"]);
ny=round(Int64,tt["ny"]);
nz=round(Int64,tt["nz"]);

# Coordinate
X=tt["X"];
Y=tt["Y"];
Z=tt["Z"];

# grid spacing
h=tt["h"];


# compute average velocity
va=0;
ne=0;
for i=1:length(s1)
    for j=1:length(r1[i])
        global va,ne;
        va=va+((s1[i][1]*h-r1[i][j]*h)^2+
        (s2[i][1]*h-r2[i][j]*h)^2+
        (s3[i][1]*h-r3[i][j]*h)^2)^.5/R_true[i][j];
        ne=ne+1;
    end
end
va=va/ne;
# initial velocity model
v=zeros(nx,ny,nz);
v[:] .=va;

# reshape the velocity to a 1D array
vc=reshape(v,nx*ny*nz,);

# initial move of sources
l1=zeros(length(s1),);
l2=zeros(length(s1),);
l3=zeros(length(s1),);
for I=1:length(l1)
    RTI.JLD2.save(string(p3,"/temp2/source_",I,".jld2"), "data",[l1[I],l2[I],l3[I]]);
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
function data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3,R_true,update_source)
    """
    Target: to implement the forward problem.
    Comments: shots are parallelized.
    Output:
    T: travel time for all events.
    chi: L2-norm of the data misfit.
    """
    foreach(rm,filter(endswith(".jld2"),readdir(string(p3,"/temp"),join=true)));
    v=reshape(vc,nx,ny,nz);
    M=[0:8:size(s1,1);size(s1,1)];
    E=zeros(length(s1,));
    for m=1:size(M,1)-1
        Threads.@threads for I=(M[m]+1):M[m+1]
            input_l1,input_l2,input_l3=RTI.JLD2.load(string(p3,"/temp2/source_",I,".jld2"))["data"];

            input_s1=s1[I][1];
            input_s2=s2[I][1];
            input_s3=s3[I][1];

            # Compute forward travel time
            T,R_cal,~=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
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
            tx,ty,tz=RTI.G(T,h);
            E[I]=.5*sum((R_cal+input_l1*tx[CartesianIndex.(r1[I],r2[I],r3[I])]
            +input_l2*ty[CartesianIndex.(r1[I],r2[I],r3[I])]
            +input_l3*tz[CartesianIndex.(r1[I],r2[I],r3[I])]-R_true[I]).^2);
            RTI.JLD2.save(string(p3,"/temp/source_",I,".jld2"), "data",T);
            if update_source==1
                # update souce move
                input_l1,input_l2,input_l3=RTI.eikonal.acoustic_eikonal_compute_source_move(nx=nx,
                ny=ny,
                nz=nz,
                h=h,
                v=v,
                T=T,
                r1=r1[I],
                r2=r2[I],
                r3=r3[I],
                s1=s1[I][1],
                s2=s2[I][1],
                s3=s3[I][1],
                R_cal=T[CartesianIndex.(r1[I],r2[I],r3[I])],
                R_true=R_true[I]);
                RTI.JLD2.save(string(p3,"/temp2/source_",I,".jld2"), "data",[input_l1,input_l2,input_l3]);
            end
        end
    end

    chi=sum(E);
    return chi,l1,l2,l3
end
## define gradient function, the gradient is computed using adjoint state method
function g!(storage,vc)
    """
    Target: to supply the routine for the update direction of the velocity.
    Comments: shots are parallelized.
    storage: the update direction.
    """
    DV=zeros(nx,ny,nz);
    M=[0:8:size(s1,1);size(s1,1)];
    for m=1:size(M,1)-1
        DV2=zeros(nx,ny,nz,length((M[m]+1):M[m+1]));
        @Threads.threads for I=(M[m]+1):M[m+1]
            temp_t=RTI.JLD2.load(string(p3,"/temp/source_",I,".jld2"))["data"];
            input_s1=s1[I][1];
            input_s2=s2[I][1];
            input_s3=s3[I][1];
            tx,ty,tz=RTI.G(temp_t,h);

            input_l1,input_l2,input_l3=RTI.JLD2.load(string(p3,"/temp2/source_",I,".jld2"))["data"];
            #=
            input_l1=0;
            input_l2=0;
            input_l3=0;
            =#
            # Compute adjoint travel time
            lambda,~=RTI.eikonal.acoustic_eikonal_adjoint(nx=nx,
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
            R_cal=temp_t[CartesianIndex.(r1[I],r2[I],r3[I])]+input_l1*tx[CartesianIndex.(r1[I],r2[I],r3[I])]+input_l2*ty[CartesianIndex.(r1[I],r2[I],r3[I])]+input_l3*tz[CartesianIndex.(r1[I],r2[I],r3[I])],
            R_true=R_true[I],
            N=ones(size(temp_t[CartesianIndex.(r1[I],r2[I],r3[I])]))*(N2));
            if m==1
                tt=0;
            else
                tt=(m-1)*length((M[m-1]+1):M[m]);
            end
            DV2[:,:,:,I-tt]=DV2[:,:,:,I-tt]-2*lambda ./v .^3;
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
N2=1;

# Test inversion.
fu=3;
E0=data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3,R_true,0);
sca=1;
test_storage=zeros(size(vc));
g!(test_storage,vc);
sca=1/maximum(abs.(test_storage));

# Perform inversion
fu=3;
opt1=optimize(vc->data_cost_L2_norm(vc,nx,ny,nz,h,s1,s2,s3,T0,r1,r2,r3,p3,R_true,1)[1],
g!,vc,LBFGS(m=5,alphaguess=LineSearches.InitialQuadratic(α0=sca*10.0,αmin=sca*.5),
linesearch=LineSearches.BackTracking(c_1=10.0^-8)),
Optim.Options(iterations=20,store_trace=true,show_trace=true,
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
for i=1:length(trace)
    trace_err[i]=parse(Float64, split(string(trace[i]))[2]);
    trace_time[i]=parse(Float64, split(string(trace[i]))[end]);
end

ax=plot(trace_time,trace_err,seriestype = :scatter,
xlabel="time [s]",
ylabel="data misfit [s^2]",yscale=:log10);
display(ax)
##
for i=1:length(l1)
    l1[i],l2[i],l3[i]=input_l1,input_l2,input_l3=RTI.JLD2.load(string(p3,"/temp2/source_",i,".jld2"))["data"];
end
.5*sum((l1-ds1*h).^2+(l2-ds2*h).^2+(l3-ds3*h).^2)
.5*sum((ds1*h).^2+(ds2*h).^2+(ds3*h).^2)
