clear all;
close all;
tt=load('./inversion_process_lbfgs_P/final/final_model.mat');
vp=tt.data.v;

X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;
h=X(2)-X(1);

tt=load('./inversion_process_lbfgs_S/final/final_model.mat');
vs=tt.data.v;

vp_vs=vp./vs;
tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_P/final/final_model.mat']);
cip=tt.data.v;
tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_S/final/final_model.mat']);
cis=tt.data.v;

tt=load(['../tutorial_alternating_checkerboard_Crati/c_vp.mat']);
cipt=tt.data.vp;

tt=load(['../tutorial_alternating_checkerboard_Crati/c_vs.mat']);
cist=tt.data.vs;

tt=load('./m.mat');
vp0=tt.data.vp;
vs0=tt.data.vs;
turbo2=flip(turbo,1);
redblue2=flip(redblue,1);
%%
% type the edge length
edge_length=6000;
%%
receiver=load('./receiver_overview.csv');
%% test conversione coord receiver
receiver(:,2)=(receiver(:,2)/1000)/111+39.08;
receiver(:,1)=(receiver(:,1)/1000)/(111*cos(39.08*pi/180))+15.93;
%% create checkerboard
[nx,ny,nz]=size(X);

tt=load('../tutorial_alternating_checkerboard_Crati/c_vp.mat');
checkerboard_P=tt.data.vp;
%% please remove this section if you already have an inverted vp and vs
%{
vp=randn(nx,ny,nz);
vs=randn(nx,ny,nz);
%}
%% slection with checkerboard P
tt=cipt(:,10,10);
tt2=diff(tt);
x_line=find(tt2~=0);
x_center=floor((x_line(1:end-1)+x_line(2:end))*.5);

tt=cipt(10,:,10);
tt2=diff(tt);
y_line=find(tt2~=0);
y_center=floor((y_line(1:end-1)+y_line(2:end))*.5);

tt=cipt(10,10,:);
tt2=diff(tt);
z_line=find(tt2~=0);
z_center=floor((z_line(1:end-1)+z_line(2:end))*.5);

s=zeros(length(x_line),...
    length(y_line), ...
    length(z_line));

s(4,4,1)=1;


% Change k from 1 to max and assign s
k=3;

figure(60)
subplot(2,2,1)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],checkerboard_P(:,:,z_center(k))');
set(gca,'YDir','normal');
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
subplot(2,2,2)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],cip(:,:,z_center(k))');
set(gca,'YDir','normal');
hold on;
for i=1:length(x_line)
    plot([X(x_line(i),1,1),X(x_line(i),1,1)],[min(Y(:)),max(Y(:))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

for i=1:length(y_line)
    plot([min(X(:)),max(X(:))],[Y(1,y_line(i),1),Y(1,y_line(i),1)],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

subplot(2,2,3)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],s(:,:,k)')
set(gca,'YDir','normal');
%% slection with checkerboard S
s2=zeros(length(x_line),...
    length(y_line), ...
    length(z_line));

s2(4,4,1)=1;

% Change k from 1 to max and assign s2
k=1;

figure(60)
subplot(2,2,1)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],cist(:,:,z_center(k))');
set(gca,'YDir','normal');
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
subplot(2,2,2)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],cis(:,:,z_center(k))');
set(gca,'YDir','normal');
hold on;
for i=1:length(x_line)
    plot([X(x_line(i),1,1),X(x_line(i),1,1)],[min(Y(:)),max(Y(:))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

for i=1:length(y_line)
    plot([min(X(:)),max(X(:))],[Y(1,y_line(i),1),Y(1,y_line(i),1)],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

subplot(2,2,3)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],s2(:,:,k)')
set(gca,'YDir','normal');
%% evaluate inversion
%{
tt2=ones(size(vp))*nan;
for i=1:nx
    for j=1:ny
        for k=1:nz
            tt4=floor(X(i,1,1)/edge_length);
            tt5=floor(Y(1,j,1)/edge_length);
            tt6=floor((Z(1,1,k)-min(Z,[],"all"))/edge_length);
            if tt4~=0 && tt5~=0 && tt6~=0 && ...
                    tt4<=size(s,1) && tt5<=size(s,2) && ...
                    tt6<=size(s,3) && s(tt4,tt5,tt6)==1
                tt2(i,j,k)=vp(i,j,k);
            end
        end
    end
end
vp=tt2;

tt2=ones(size(vs))*nan;
for i=1:nx
    for j=1:ny
        for k=1:nz
            tt4=floor(X(i,1,1)/edge_length);
            tt5=floor(Y(1,j,1)/edge_length);
            tt6=floor((Z(1,1,k)-min(Z,[],"all"))/edge_length);
            if tt4~=0 && tt5~=0 && tt6~=0 && ...
                    tt4<=size(s2,1) && tt5<=size(s2,2) && ...
                    tt6<=size(s2,3) && s2(tt4,tt5,tt6)==1
                tt2(i,j,k)=vs(i,j,k);
            end
        end
    end
end
vs=tt2;
%}
%% save evaluate inversion
tt=ones(size(vp));
tt(isnan(vp))=nan;
data=tt;
save('P_selction.mat','data');
tt=ones(size(vs));
tt(isnan(vs))=nan;
save('s_selction.mat','data');
%% compute relative v
clear data;
data.X=X;
data.Y=Y;
data.Z=Z;
data.nx=nx;
data.ny=ny;
data.nz=nz;
data.h=h;
data.r_vp=vp-vp0;
save('r_vp.mat','data');

clear data;
data.X=X;
data.Y=Y;
data.Z=Z;
data.nx=nx;
data.ny=ny;
data.nz=nz;
data.h=h;
data.r_vs=vs-vs0;
save('r_vs.mat','data');
%% test conversione coordinate
Y=(Y/1000)/111+39.08;
X=(X/1000)/(111*cos(39.08*pi/180))+15.93;
%% vp horizontal xy plane
%z=-2000:-2000:-20000;
z=-2000:-2000:-16000;
%a=ceil(sqrt(length(z)));

figure;
set(gcf,'position',[80,80,1000,1500]);
% compute mean
vp0=zeros(size(vp));
for k=1:size(vp0,3)
    vp0(:,:,k)=mean(vp(:,:,k),'all');
end

%limps=[1.5,2.1];
for i=1:length(z)
    tt=find(abs(Z(1,1,:)-z(i))==min(abs(Z(1,1,:)-z(i))));
    tt=tt(1);
    %subplot(a,a,i)
    subplot(5,3,i)
    imAlpha=ones(size(vp(:,:,tt)'));
    imAlpha(isnan(vp(:,:,tt)'))=0;
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        vp(:,:,tt)'-vp0(:,:,tt)','AlphaData',imAlpha);
    %imagesc(vp(:,:,tt)'-vp0(:,:,tt)','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp z =' num2str(z(i)) 'm']});
    
    %xlabel('X [m]');
    %ylabel('Y [m]');
    %caxis(limp);
    caxis([-600 600])
    colormap(redblue2);
    if i>6 %8
        xlabel('lon [deg]');
        cb=colorbar('southoutside');
        cb.Position=cb.Position+[-0.007 -0.09 0.009 0.009];
    end
    %if i==1 || i==5 || i==9
    if i==1 || i==4 || i==7 || i==10 || i==13
        ylabel('lat [deg]');
    end
%     xlim([16.1 16.5])
%     ylim([39.2 39.6])
    %colorbar;
    hold on;
    %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');    
end
print(gcf,['./vp_xy'],'-djpeg','-r400');
%% vp vertical yz plane
% x=3000:10000:44800;
% a=ceil(sqrt(length(x)));
% 
% figure;
% set(gcf,'position',[80,80,1500,1500]);
% 
% %limps=[1.5,2.1];
% for i=1:length(x)
%     tt=find(abs(X(:,1,1)-x(i))==min(abs(X(:,1,1)-x(i))));
%     tt=tt(1);
%     subplot(a,a,i)
%     imAlpha=reshape(ones(size(vp(tt,:,:))),[ny,nz])';
%     imAlpha(isnan(reshape(vp(tt,:,:),[ny,nz])'))=0;
%     imagesc([min(Y(:)),max(Y(:))], ...
%         [min(Z(:)),max(Z(:))], ...
%         reshape(vp(tt,:,:)-vp0(tt,:,:),[ny,nz])','AlphaData',imAlpha);
%     set(gca,'color',1*[1 1 1]);
%     set(gca,'ydir','normal');
%     title({['vp x =' num2str(x(i)) 'm']});
%     
%     xlabel('Y [m]');
%     ylabel('Z [m]');
%     %caxis(limp);
%     colormap(redblue2);
%     colorbar;
%     %hold on;
%     %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
% end
% print(gcf,['./vp_yz'],'-djpeg','-r400');
% %% vp vertical xz plane
% y=3000:10000:60400;
% a=ceil(sqrt(length(x)));
% 
% figure;
% set(gcf,'position',[80,80,1500,1500]);
% 
% %limps=[1.5,2.1];
% for i=1:length(y)
%     tt=find(abs(Y(1,:,1)-y(i))==min(abs(Y(1,:,1)-y(i))));
%     tt=tt(1);
%     subplot(a,a,i)
%     imAlpha=reshape(ones(size(vp(:,tt,:))),[nx,nz])';
%     imAlpha(isnan(reshape(vp(:,tt,:),[nx,nz])'))=0;
%     imagesc([min(X(:)),max(X(:))], ...
%         [min(Z(:)),max(Z(:))], ...
%         reshape(vp(:,tt,:)-vp0(:,tt,:),[nx,nz])','AlphaData',imAlpha);
%     set(gca,'color',1*[1 1 1]);
%     set(gca,'ydir','normal');
%     title({['vp y =' num2str(y(i)) 'm']});
%     
%     xlabel('X [m]');
%     ylabel('Z [m]');
%     %caxis(limp);
%     colormap(redblue2);
%     colorbar;
%     %hold on;
%     %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
% end
% print(gcf,['./vp_xz'],'-djpeg','-r400');
%% vs horizontal xy plane
z=-2000:-2000:-20000;
a=ceil(sqrt(length(z)));

figure;
set(gcf,'position',[80,80,1500,1500]);

vs0=zeros(size(vs));
for k=1:size(vp0,3)
    vs0(:,:,k)=mean(vs(:,:,k),'all');
end

%limps=[1.5,2.1];
for i=1:length(z)
    tt=find(abs(Z(1,1,:)-z(i))==min(abs(Z(1,1,:)-z(i))));
    tt=tt(1);
    j=1;
    subplot(a,a,i)
    imAlpha=ones(size(vs(:,:,tt)'));
    imAlpha(isnan(vs(:,:,tt)'))=0;
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        vs(:,:,tt)'-vs0(:,:,tt)','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vs z =' num2str(z(i)) 'm']});
    
    %xlabel('X [m]');
    %ylabel('Y [m]');
    %caxis(limp);
    caxis([-400 400]);
    colormap(redblue2);
    if i>8
        xlabel('lon [deg]');
        cb=colorbar('southoutside');
        cb.Position=cb.Position+[-0.007 -0.09 0.009 0.009];
    end
    if i==1 || i==5 || i==9
        ylabel('lat [deg]');
    end
    %colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vs_xy'],'-djpeg','-r400');
%% vs vertical yz plane
% x=3000:10000:44800;
% a=ceil(sqrt(length(x)));
% 
% figure;
% set(gcf,'position',[80,80,1500,1500]);
% 
% %limps=[1.5,2.1];
% for i=1:length(x)
%     tt=find(abs(X(:,1,1)-x(i))==min(abs(X(:,1,1)-x(i))));
%     tt=tt(1);
%     subplot(a,a,i)
%     imAlpha=reshape(ones(size(vs(tt,:,:))),[ny,nz])';
%     imAlpha(isnan(reshape(vs(tt,:,:),[ny,nz])'))=0;
%     imagesc([min(Y(:)),max(Y(:))], ...
%         [min(Z(:)),max(Z(:))], ...
%         reshape(vs(tt,:,:)-vs0(tt,:,:),[ny,nz])','AlphaData',imAlpha);
%     set(gca,'color',1*[1 1 1]);
%     set(gca,'ydir','normal');
%     title({['vs x =' num2str(x(i)) 'm']});
%     
%     xlabel('Y [m]');
%     ylabel('Z [m]');
%     %caxis(limp);
%     colormap(redblue2);
%     colorbar;
%     %hold on;
%     %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
% end
% print(gcf,['./vs_yz'],'-djpeg','-r400');
% %% vs vertical xz plane
% y=3000:10000:60400;
% a=ceil(sqrt(length(x)));
% 
% figure;
% set(gcf,'position',[80,80,1500,1500]);
% 
% %limps=[1.5,2.1];
% for i=1:length(y)
%     tt=find(abs(Y(1,:,1)-y(i))==min(abs(Y(1,:,1)-y(i))));
%     tt=tt(1);
%     subplot(a,a,i)
%     imAlpha=reshape(ones(size(vs(:,tt,:))),[nx,nz])';
%     imAlpha(isnan(reshape(vs(:,tt,:),[nx,nz])'))=0;
%     imagesc([min(X(:)),max(X(:))], ...
%         [min(Z(:)),max(Z(:))], ...
%         reshape(vs(:,tt,:)-vs0(:,tt,:),[nx,nz])','AlphaData',imAlpha);
%     set(gca,'color',1*[1 1 1]);
%     set(gca,'ydir','normal');
%     title({['vs y =' num2str(y(i)) 'm']});
%     
%     xlabel('X [m]');
%     ylabel('Z [m]');
%     %caxis(limp);
%     colormap(redblue2);
%     colorbar;
%     %hold on;
%     %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
% end
% print(gcf,['./vs_xz'],'-djpeg','-r400');

%% vp_vs horizontal xy plane
z=-2000:-2000:-20000;
a=ceil(sqrt(length(z)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(z)
    tt=find(abs(Z(1,1,:)-z(i))==min(abs(Z(1,1,:)-z(i))));
    tt=tt(1);
    j=1;
    subplot(a,a,i)
    imAlpha=ones(size(vp_vs(:,:,tt)'));
    imAlpha(isnan(vp_vs(:,:,tt)'))=0;
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        vp_vs(:,:,tt)','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp / vs z =' num2str(z(i)) 'm']});
    
%     xlabel('X [m]');
%     ylabel('Y [m]');
    %caxis(limp);
    caxis([1.6 1.9])
    colormap(redblue2);
    if i>8
        xlabel('lon [deg]');
        cb=colorbar('southoutside');
        cb.Position=cb.Position+[-0.007 -0.09 0.009 0.009];
    end
    if i==1 || i==5 || i==9
        ylabel('lat [deg]');
    end

    %colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vp_vs_xy'],'-djpeg','-r400');