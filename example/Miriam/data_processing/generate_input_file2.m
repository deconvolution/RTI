close all;
clear all;

origin=[-2.85,35.85,-19000];
h=800;
x_min=3;
y_min=3;
z_min=3;
x_max=3;
y_max=3;
z_max=3;
%%
p2='./temp/';
listing=dir(p2);
p3='./traveltime_input/';
mkdir(p3);
%% for input
for i=1:size(listing,1)-2
    %% generate empty file
    data.Rp=[];
    data.Rs=[];
    data.S=[];
    tt=load([p2,listing(i+2).name]);
    if 1
        %% S
        [x2,y2]=latlon2local(tt.data.S(:,2),tt.data.S(:,1),tt.data.S(:,3),origin);
        DX=x2;
        DY=y2;
        data.S=[round(DX/h),round(DY/h),round((tt.data.S(:,3)-origin(3))/h)];
        %%
        x_min=min(x_min,data.S(1));
        y_min=min(y_min,data.S(2));
        z_min=min(z_min,data.S(3));

        x_max=max(x_max,data.S(1));
        y_max=max(y_max,data.S(2));
        z_max=max(z_max,data.S(3));
        %% Rp
        [x2,y2]=latlon2local(tt.data.Rp(:,2),tt.data.Rp(:,1),tt.data.Rp(:,3),origin);
        DX=x2;
        DY=y2;
        data.Rp=[round(DX/h),round(DY/h),round((tt.data.Rp(:,3)-origin(3))/h),tt.data.Rp(:,4)];
        %%
        x_min=min(x_min,min(data.Rp(:,1)));
        y_min=min(y_min,min(data.Rp(:,2)));
        z_min=min(z_min,min(data.Rp(:,3)));

        x_max=max(x_max,max(data.Rp(:,1)));
        y_max=max(y_max,max(data.Rp(:,2)));
        z_max=max(z_max,max(data.Rp(:,3)));
         %% Rs
        [x2,y2]=latlon2local(tt.data.Rs(:,2),tt.data.Rs(:,1),tt.data.Rs(:,3),origin);
        DX=x2;
        DY=y2;
        data.Rs=[round(DX/h),round(DY/h),round((tt.data.Rs(:,3)-origin(3))/h),tt.data.Rs(:,4)];
        %%
        x_min=min(x_min,min(data.Rs(:,1)));
        y_min=min(y_min,min(data.Rs(:,2)));
        z_min=min(z_min,min(data.Rs(:,3)));

        x_max=max(x_max,max(data.Rs(:,1)));
        y_max=max(y_max,max(data.Rs(:,2)));
        z_max=max(z_max,max(data.Rs(:,3)));
        %%
        save([p3 listing(i+2).name],'data');
    end
end
%% verify data
listing=dir(p3);
data3=[];
data4=[];
for i=1:size(listing,1)-2
    tt=load([p3,listing(i+2).name]);
    data4=[data4;
        ((tt.data.Rp(:,1)-tt.data.S(1)).^2+(tt.data.Rp(:,2)-tt.data.S(2)).^2 ...
        +(tt.data.Rp(:,3)-tt.data.S(3)).^2).^.5*h,tt.data.Rp(:,4)];
    data3=[data3;
        ((tt.data.Rs(:,1)-tt.data.S(1)).^2+(tt.data.Rs(:,2)-tt.data.S(2)).^2 ...
        +(tt.data.Rs(:,3)-tt.data.S(3)).^2).^.5*h,tt.data.Rs(:,4)];
end
figure;
ax=plot(data4(:,1),data4(:,2),'.','color','red');
hold on;
ax2=plot(data3(:,1),data3(:,2),'.','color','blue');
legend([ax,ax2],'P','S','Location','northwest','orientation','horizontal');
xlabel('distance [m]');
ylabel('t [s]');