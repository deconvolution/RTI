%%
%{
This script generates information for each event.
%}
%%
clear all;
close all;
%% read file
tt=load('./picks_focmecs.mat');
data2=tt.picks;
%% error data
data2(find(data2.stat_lon==0 & data2.stat_lat==0),:)=[];
%% create a folder to store temp
p2='./temp/';
mkdir(p2);
%% assign data
data.S=[];
data.Rp=[];
data.Rs=[];
ns=1;
i=1;
while i<=size(data2,1)
    data.S=[data2(i,3).long,data2(i,2).lat,data2(i,4).depth];
    tt=max(find((data.S(1)==data2.long) & (data.S(2)==data2.lat) & (data.S(3)==data2.depth) ...
        & (data2(i,1).source_time==data2.source_time) ));
    ind_p=find(data2(i:tt,:).phase=="IP");
    ind_s=find(data2(i:tt,:).phase=="IS");
    tt2=data2(i:tt,:);
    data.Rp=[tt2(ind_p,:).stat_lon,tt2(ind_p,:).stat_lat,tt2(ind_p,:).stat_h,seconds(tt2(ind_p,:).relative_time)];
    data.Rs=[tt2(ind_s,:).stat_lon,tt2(ind_s,:).stat_lat,tt2(ind_s,:).stat_h,seconds(tt2(ind_s,:).relative_time)];
    data.S(3)=data.S(3)*1000;
    data.Rp(:,3)=data.Rp(:,3)*1000;
    data.Rs(:,3)=data.Rs(:,3)*1000;
    if data.S(3)>=-17000
        save([p2 num2str(ns) '.mat']);
    end
    data.S=[];
    data.Rp=[];
    data.Rs=[];
    i=tt+1;
    ns=ns+1;
end
%% plot source distribution
figure;
plot(data2.long,data2.lat,'v','color','red');
hold on;
plot(data2.stat_lon,data2.stat_lat,'^','color','green');

figure;
plot(data2.depth,'.');
%% verify data
listing=dir(p2);
data3=[];
data4=[];
for i=1:size(listing,1)-2
    tt=load([p2,listing(i+2).name]);
    for j=1:size(tt.data.Rp,1)
        origin=[tt.data.S(2),tt.data.S(1),tt.data.S(3)];
        [x,y,z]=latlon2local(tt.data.Rp(j,2),tt.data.Rp(j,1),tt.data.Rp(j,3),origin);
        data4=[data4;
            (x^2+y^2+z^2)^.5,tt.data.Rp(j,4)];
    end
    for j=1:size(tt.data.Rs,1)
        origin=[tt.data.S(2),tt.data.S(1),tt.data.S(3)];
        [x,y,z]=latlon2local(tt.data.Rs(j,2),tt.data.Rs(j,1),tt.data.Rs(j,3),origin);
        data3=[data3;
            (x^2+y^2+z^2)^.5,tt.data.Rs(j,4)];
    end
end
figure;
ax=plot(data4(:,1),data4(:,2),'.','color','red');
hold on;
ax2=plot(data3(:,1),data3(:,2),'.','color','blue');
xlabel('distance [m]');
ylabel('t [s]');

legend([ax,ax2],'P','S','Location','northwest','orientation','horizontal');