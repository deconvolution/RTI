##
"
Compute source move.
"

function acoustic_eikonal_compute_source_move2(;nx,ny,nz,h,v,T,s1,s2,s3,r1,r2,r3,
    R_cal,R_true)

    tx,ty,tz=acoustic_eikonal_source_perturbation(tau=T,h=h);

    dl1=0;
    dl2=0;
    dl3=0;
    for i=1:length(r1)
        dl1=dl1+tx[r1[i],r2[i],r3[i]]*(R_cal[i]-R_true[i]);
        dl2=dl2+ty[r1[i],r2[i],r3[i]]*(R_cal[i]-R_true[i]);
        dl3=dl3+tz[r1[i],r2[i],r3[i]]*(R_cal[i]-R_true[i]);
    end
    return dl1,dl2,dl3
end
