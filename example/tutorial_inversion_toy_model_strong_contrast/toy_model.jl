using RTI,Random,CSV,DataFrames
## generate toy model
# dimensions
nx=60;
ny=50;
nz=55;
# grid spacing
h=10;

# create a checkerboard
w=zeros(nx,ny,nz);
edge_length=150;
w=RTI.cb(h,nx,ny,nz,edge_length);

# create a Coordinate
Y,X,Z=RTI.meshgrid(h*(1:ny),h*(1:nx),h*(1:nz));

# createa a model
v=ones(nx,ny,nz);
v=1000 .+200*w;
## define observations
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
## Create empty vectors to store source and receivers
Random.seed!(1);
tt=rand(2:nx-1,90,);
tt2=rand(2:ny-1,90,);
tt3=ones(size(tt))*2;
for I=1:100
    global R_true,s1,s2,s3,r1,r2,r3;
    s1=push!(s1,round.(Int64,rand(2:nx-1,1,)));
    s2=push!(s2,round.(Int64,rand(2:ny-1,1,)));
    s3=push!(s3,round.(Int64,rand(20:nz-1,1,)));
    r1=push!(r1,round.(Int64,tt));
    r2=push!(r2,round.(Int64,tt2));
    r3=push!(r3,round.(Int64,tt3));
end
## create folders
p2=@__FILE__;
p3=replace(p2,".jl"=>"");
if isdir(p3)==0
    mkdir(p3);
end
if isdir(string(p3,"/../obs"))==0
    mkdir(string(p3,"/../obs"))
end
## save toy model
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
data=data3(0,0,0,0,0,0,0,0)
data.X=X;
data.Y=Y;
data.Z=Z;
data.nx=nx;
data.ny=ny;
data.nz=nz;
data.v=v;
data.h=h;
file=RTI.matopen(string(p3,"/true_model.mat"), "w");
write(file,"data",data);
close(file);
## write sources and receivers
r1t=r1*h;
r2t=r2*h;
r3t=r3*h;
s1t=zeros(length(s1),);
s2t=zeros(length(s2),);
s3t=zeros(length(s3),);
for i=1:length(s1t)
    s1t[i]=h*s1[i][1];
    s2t[i]=h*s2[i][1];
    s3t[i]=h*s3[i][1];
end

CSV.write(string(p3,"/receiver location.csv"),DataFrame([r1t[1] r2t[1] r3t[1]],:auto));
CSV.write(string(p3,"/source location.csv"),DataFrame([s1t s2t s3t],:auto));
## write true model to vtk
vtkfile=RTI.vtk_grid(string(p3,"/true_model"),X,Y,Z);
vtkfile["v"]=reshape(v,nx,ny,nz);
RTI.vtk_save(vtkfile);
## generate true observations
mutable struct data4
    S
    Rp
    Rs
end
data=data4(0,0,0);

M=[0:50:size(s1,1);size(s1,1)];
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
        T0=0,
        r1=r1[I],
        r2=r2[I],
        r3=r3[I]);

        data.S=[input_s1 input_s2 input_s3];
        data.Rp=[reshape(r1[I],length(r1[I]),1) reshape(r2[I],length(r2[I]),1) reshape(r3[I],length(r3[I]),1) reshape(R_cal,length(R_cal),1)];
        file=RTI.matopen(string(p3,"/../obs/source_",I,".mat"), "w");
        write(file,"data",data);
        close(file);
    end
end
