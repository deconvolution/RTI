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
v[:] .=1;
v[:,:,60:end] .=2;
# receiver grid location
r1=[10,50,80,90,80];
r2=[50,50,50,20,80];
r3=[5,5,5,5,5];
# source grid location
s1=90;
s2=90;
s3=100;

# Compute forward travel time
T,R_true=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
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
##
s1=s1-20;
s2=s2-10;
s3=s3-10;

s1_range=[2,nx-1];
s2_range=[2,ny-1];
s3_range=[2,nz-1];
s1,s2,s3,=RTI.eikonal.acoustic_eikonal_update_source_location(nx=nx,
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
##
