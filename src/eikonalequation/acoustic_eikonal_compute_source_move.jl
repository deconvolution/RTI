##
"
Compute source move.
"

function acoustic_eikonal_compute_source_move(;nx,ny,nz,h,v,T,s1,s2,s3,r1,r2,r3,
    R_cal,R_true)

    tau=T;
    tx,ty,tz=acoustic_eikonal_source_perturbation(tau=tau,h=h);

    Ds1=zeros(length(r1),);
    Ds2=zeros(length(r2),);
    Ds3=zeros(length(r3),);

    for i=1:length(r1)
        Ds1[i]=tx[r1[i],r2[i],r3[i]];
        Ds2[i]=ty[r1[i],r2[i],r3[i]];
        Ds3[i]=tz[r1[i],r2[i],r3[i]];
    end

    l1,l2,l3=[reshape(Ds1,length(Ds1),1) reshape(Ds2,length(Ds2),1) reshape(Ds3,length(Ds3),1)]\(R_true-R_cal);
    return l1,l2,l3
end
