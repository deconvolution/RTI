## using packages
using RTI,CSV,DataFrames
## define parameters for the forward problem
# grid dimensions
nx=100;
ny=90;
nz=120;
# grid spacing
h=1;
# velocity
v=zeros(nx,ny,nz);
for i=1:nz
    v[:,:,i] .=1000+1*i;
end

# receiver grid location
r1=[10,50,80,90,80];
r2=[50,50,50,20,80];
r3=[5,5,5,5,5];
# source grid location
s1d=90;
s2d=90;
s3d=100;

# Compute forward travel time, this is the true model
T,R_true=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=s1d,
s2=s2d,
s3=s3d,
T0=0,
r1=r1,
r2=r2,
r3=r3,
n_iteration=2,tol=.1);
## Iterative source locating
s1=s1d-3;
s2=s2d-5;
s3=s3d-4;

s1_range=[2,nx-1];
s2_range=[2,ny-1];
s3_range=[2,nz-1];
s1i,s2i,s3i,E_source_location=RTI.eikonal.acoustic_eikonal_compute_source_location(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=s1,
s2=s2,
s3=s3,
r1=r1,
r2=r2,
r3=r3,
s1_range=s1_range,
s2_range=s2_range,
s3_range=s3_range,
n_iteration=10,
R_true=R_true);
## 1-step source locating
T,R_cal,nt=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=s1,
s2=s2,
s3=s3,
T0=0,
r1=r1,
r2=r2,
r3=r3);

l1,l2,l3=RTI.eikonal.acoustic_eikonal_compute_source_move(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
T=T,
s1=s1,
s2=s2,
s3=s3,
r1=r1,
r2=r2,
r3=r3,
R_cal=R_cal,
R_true=R_true);
