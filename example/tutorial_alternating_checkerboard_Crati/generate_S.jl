using RTI,Random,CSV,DataFrames,MAT
## create path for inversions
p2=@__FILE__;
p3=replace(p2,".jl"=>"");
if isdir(p3)==0
    mkdir(p3);
end

## define empty arrays for loading data
# observed data
R_true=Vector{Vector{Float64}}();
# source grid location data
s1d=Vector{Vector{Int64}}();
s2d=Vector{Vector{Int64}}();
s3d=Vector{Vector{Int64}}();
# source grid location initial guess
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
    global R_true,s1d,s2d,s3d,r1,r2,r3;
    tt2=RTI.readmat(string("./obs/",tt[I]),"data");
    if length(tt2["Rs"])!=0
        R_true=push!(R_true,tt2["Rs"][:,4]);
        s1d=push!(s1d,round.(Int64,tt2["S"][:,1]));
        s2d=push!(s2d,round.(Int64,tt2["S"][:,2]));
        s3d=push!(s3d,round.(Int64,tt2["S"][:,3]));
        r1=push!(r1,round.(Int64,tt2["Rs"][:,1]));
        r2=push!(r2,round.(Int64,tt2["Rs"][:,2]));
        r3=push!(r3,round.(Int64,tt2["Rs"][:,3]));
    end
end
s1=copy(s1d);
s2=copy(s2d);
s3=copy(s3d);
##
R_cal=copy(R_true);
## define parameters for the forward problem
tt=RTI.readmat("./m.mat","data");
# number of grids
nx=round(Int64,tt["nx"]);
ny=round(Int64,tt["ny"]);
nz=round(Int64,tt["nz"]);
# create a checkerboard
w=zeros(nx,ny,nz);
# create a Coordinate
X=tt["X"];
Y=tt["Y"];
Z=tt["Z"];
h=tt["dx"];

# Change the edge length of each block
edge_length=6000;
w=RTI.cb(h,nx,ny,nz,edge_length);

v=zeros(nx,ny,nz);
v[:] .=3000;

v=v+300*w;

data=copy(tt);
data["vs"]=v;

file=RTI.matopen("./c_vs.mat", "w");
write(file,"data",data);
close(file);
## write final model to vtk
vtkfile=RTI.vtk_grid(string("./c_vs"),X,Y,Z);
vtkfile["v"]=v;
RTI.vtk_save(vtkfile);
##
data=nothing;
for I=1:length(file_name)
    tt=readdir("./obs/");
    tt2=RTI.readmat(string("./obs/",tt[I]),"data");
    if length(tt2["Rs"])!=0
        T,R_cal,~=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
        ny=ny,
        nz=nz,
        h=h,
        v=v,
        s1=round.(Int64,tt2["S"][:,1])[1],
        s2=round.(Int64,tt2["S"][:,2])[1],
        s3=round.(Int64,tt2["S"][:,3])[1],
        T0=T0,
        r1=round.(Int64,tt2["Rs"][:,1]),
        r2=round.(Int64,tt2["Rs"][:,2]),
        r3=round.(Int64,tt2["Rs"][:,3]));

        data=copy(tt2);
        data["Rs"][:,4]=R_cal;
        file=RTI.matopen(string(p3,"/",file_name[I]), "w");
        write(file,"data",data);
        close(file);
        println("\n  l=",I);
    end
end
