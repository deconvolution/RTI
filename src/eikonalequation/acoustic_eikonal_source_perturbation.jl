##
"
Eikonal equation source perturbation
"

function acoustic_eikonal_source_perturbation(;tau,h)
    tx=zeros(size(tau));
    ty=zeros(size(tau));
    tz=zeros(size(tau));
    tx[1:end-1,:,:]=(tau[2:end,:,:]-tau[1:end-1,:,:])/h;
    ty[:,1:end-1,:]=(tau[:,2:end,:]-tau[:,1:end-1,:])/h;
    tz[:,:,1:end-1]=(tau[:,:,2:end]-tau[:,:,1:end-1])/h;
    return tx,ty,tz
end
