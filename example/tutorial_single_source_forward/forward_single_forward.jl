## using packages
using RTI,CSV
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
r1=[10,50,80];
r2=[50,50,50];
r3=[5,5,5];
# source grid location
s1=90;
s2=90;
s3=100;
# boundary condition for the source
T0=0;

# Compute forward travel time
T,R_cal=RTI.eikonal.acoustic_eikonal_forward(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=s1,
s2=s2,
s3=s3,
T0=T0,
r1=r1,
r2=r2,
r3=r3);

## Output
# travel time at receivers
display(R_cal)
# save travel time to vts
Y3D,X3D,Z3D=RTI.meshgrid(1:ny,1:nx,1:nz);
vtkfile=RTI.vtk_grid(string("./forward"),X3D,Y3D,Z3D);
vtkfile["v"]=reshape(v,nx,ny,nz);
vtkfile["T"]=reshape(T,nx,ny,nz);
RTI.vtk_save(vtkfile);
## create CSV for source and receivers
r1t=r1*h;
r2t=r2*h;
r3t=r3*h;
s1t=s1*h;
s2t=s2*h;
s3t=s3*h;
CSV.write(string("./receiver location.csv"),DataFrame([r1t r2t r3t],:auto));
CSV.write(string("./source location.csv"),DataFrame([s1t s2t s3t],:auto));
