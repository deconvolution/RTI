##
"
Update source location.
"

function acoustic_eikonal_compute_source_location(;nx,ny,nz,h,v,s1,s2,s3,r1,r2,r3,
    s1_range=nothing,s2_range=nothing,s3_range=nothing,
    n_iteration,R_true)

    E=NaN*ones(n_iteration,);
    l=1;
    for l=1:n_iteration
        s1_last=s1;
        s2_last=s2;
        s3_last=s3;

        t,R_cal=acoustic_eikonal_forward(nx=nx,
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

        tx,ty,tz=acoustic_eikonal_source_perturbation(tau=t,h=h);

        Ds1=zeros(length(r1),);
        Ds2=zeros(length(r2),);
        Ds3=zeros(length(r3),);

        for i=1:length(r1)
            Ds1[i]=tx[r1[i],r2[i],r3[i]];
            Ds2[i]=ty[r1[i],r2[i],r3[i]];
            Ds3[i]=tz[r1[i],r2[i],r3[i]];
        end

        E[l]=.5*sum((R_cal-R_true).^2);

        l1,l2,l3=[reshape(Ds1,length(Ds1),1) reshape(Ds2,length(Ds2),1) reshape(Ds3,length(Ds3),1)]\(R_true-R_cal);

        s1=s1-round(Int64,l1/h);
        s2=s2-round(Int64,l2/h);
        s3=s3-round(Int64,l3/h);

        if s1_range!==nothing
            if s1<s1_range[1]
                s1=s1_range[1];
            else
                if s1>s1_range[2]
                    s1=s1_range[2];
                end
            end
        end

        if s2_range!==nothing
            if s2<s2_range[1]
                s2=s2_range[1];
            else
                if s2>s2_range[2]
                    s2=s2_range[2];
                end
            end
        end

        if s3_range!==nothing
            if s3<s3_range[1]
                s3=s3_range[1];
            else
                if s3>s3_range[2]
                    s3=s3_range[2];
                end
            end
        end
        if s1==s1_last && s2==s2_last && s3==s3_last
            break;
        end
    end

    return s1,s2,s3,E
end
