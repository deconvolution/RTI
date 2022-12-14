## mini utilities
## using packages
export meshgrid,write2mat,readmat,rickerWave,LinearAlgebra,JLD2
## import packages
using Random,MAT,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,
LinearAlgebra,JLD2
##
"
Same with meshgrid in Matlab

X,Y=meshgrid(x,y) returns 2-D grid coordinates based on the
coordinates contained in vectors x and y. X is a matrix where each row
is a copy of x, and Y is a matrix where each column is a copy of y. The
grid represented by the coordinates X and Y has length(y) rows and
length(x) columns.

X,Y,Z=meshgrid(x,y,z) returns 3-D grid coordinates defined by the
vectors x, y, and z. The grid represented by X, Y, and Z has size
length(y)-by-length(x)-by-length(z).
"
function meshgrid(x,y,z=1)
    if z==1
        x2=zeros(length(x),length(y));
        y2=copy(x2);
        x2=repeat(reshape(x,1,length(x)),length(y),1);
        y2=repeat(reshape(y,length(y),1),1,length(x));
        z2=nothing;
    end
    if z!=1
        x2=zeros(length(x),length(y),length(z));
        y2=copy(x2);
        x2=repeat(reshape(x,1,length(x),1),length(y),1,length(z));
        y2=repeat(reshape(y,length(y),1,1),1,length(x),length(z));
        z2=repeat(reshape(z,1,1,length(z)),length(y),length(x),1);
    end
    return x2,y2,z2
end
##
"
write a mat file.

Example: write2mat(path,var).
"
function write2mat(path,var)
    file=matopen(path,"w");
    write(file,"data",data);
    close(file);
    return nothing
end
##
"
read a mat file.

Example: readmat(path,var).
"
function readmat(path,var)
    file=matopen(path);
    tt=read(file,var);
    close(file)
    return tt
end
##
"
ricker wavelet.

Example:rickerWave(frequency,dt,nt,magnitude).
"
function rickerWave(freq,dt,ns,M)
    ## calculate scale
    E=10 .^(5.24+1.44 .*M);
    s=sqrt(E.*freq/.299);

    t=dt:dt:dt*ns;
    t0=1 ./freq;
    t=t .-t0;
    ricker=s .*(1 .-2*pi^2*freq .^2*t .^2).*exp.(-pi^2*freq^2 .*t .^2);
    ricker=ricker;
    ricker=Float32.(ricker);
    return ricker
end
##
"
Create checkerboard.

Example: cb(h,nx,ny,nz,edge_length)
h: grid spacing.
nx,ny,nz: Dimensions.
edge_length: edge length of each block.
"
function cb(h,nx,ny,nz,edge_length)
    Y,X,Z=RTI.meshgrid(h*(1:ny),h*(1:nx),h*(1:nz));
    wx=sin.(pi/edge_length*(X .-minimum(X)));
    wy=sin.(pi/edge_length*(Y .-minimum(Y)));
    wz=sin.(pi/edge_length*(Z .-minimum(Z)));
    w=(wx.*wy.*wz);
    w[w .>=0] .=1;
    w[w .<0] .=-1;
    return w
end
##
"
Compute 3D gradient
"
function G(u,h)
    nx,ny,nz=size(u);
    tx=zeros(nx,ny,nz);
    tx[1:end-1,:,:]=(u[2:end,:,:]-u[1:end-1,:,:])/h;
    ty=zeros(nx,ny,nz);
    ty[:,1:end-1,:]=(u[:,2:end,:]-u[:,1:end-1,:])/h;
    tz=zeros(nx,ny,nz);
    tz[:,:,1:end-1]=(u[:,:,2:end]-u[:,:,1:end-1])/h;
    return tx,ty,tz
end
