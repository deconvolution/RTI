clear all;
close all;
tt=load('./inversion_process_lbfgs_P/final/final_model.mat');
vp=tt.data.v;

X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;
h=X(2)-X(1);
tt=load('./inversion_process_lbfgs_P/final/final_model.mat');
vs=tt.data.v;

vp_vs=vp./vs;
tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_P/final/final_model.mat']);
cip=tt.data.v;
tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_S/final/final_model.mat']);
cis=tt.data.v;

tt=load('./m.mat');
vp0=tt.data.vp;
vs0=tt.data.vs;
turbo2=flip(turbo,1);
redblue2=flip(redblue,1);
%%
% type the edge length
edge_length=5000;
%%
receiver=load('./receiver_overview.csv');
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
n=round(edge_length/h);
s=zeros(length(1:floor(nx/n)+1) ...
    ,length(1:floor(ny/n)+1), ...
    length(1:floor(nz/n)));

s(4,4,1)=1;
s(4,5,1)=1;
s(4,6,1)=1;
s(5,2,1)=1;
s(5,3,1)=1;
s(5,4,1)=1;
s(5,5,1)=1;

s(3,3,2)=1;
s(3,4,2)=1;
s(3,5,2)=1;
s(4,4,2)=1;
s(4,5,2)=1;
s(4,6,2)=1;
s(4,7,2)=1;
s(5,2,2)=1;
s(5,3,2)=1;
s(5,4,2)=1;
s(5,5,2)=1;
s(5,6,2)=1;
s(5,7,2)=1;
s(6,3:6,2)=1;

s(3,4,3)=1;
s(3,3,4)=1;
s(5,5,5)=1;
s(6,6,6)=1;


% Change k from 1 to max and assign s
k=6;

tt=k*n-floor(n/2);

figure(60)
subplot(2,2,1)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],checkerboard_P(:,:,end-tt)');
set(gca,'YDir','normal');
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
subplot(2,2,2)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],cip(:,:,end-tt)');
set(gca,'YDir','normal');
hold on;
for i=1:n:nx
    plot([X(i,1,1),X(i,1,1)],[min(Y(:)),max(Y(:))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

for i=1:n:ny
    plot([min(X(:)),max(X(:))],[(Y(1,i,1)),(Y(1,i,1))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

subplot(2,2,3)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],s(:,:,k)')
set(gca,'YDir','normal');
%% slection with checkerboard S
tt=load('../tutorial_alternating_checkerboard_Crati/c_vs.mat');
checkerboard_S=tt.data.vs;

n=round(edge_length/h);
s2=zeros(length(1:floor(nx/n)+1) ...
    ,length(1:floor(ny/n)+1), ...
    length(1:floor(nz/n)));

s2(5,3,1)=1;

s2(4,4:5,2)=1;
s2(5,2:5,2)=1;
s2(6,2:6,2)=1;

s2(3,2,3)=1;
s2(4,6,3)=1;
s2(5,6:7,3)=1;


% Change k from 1 to max and assign s
k=2;

tt=k*n-floor(n/2);

figure(1)
subplot(2,2,1)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],checkerboard_S(:,:,end-tt)');
set(gca,'YDir','normal');
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
subplot(2,2,2)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],cis(:,:,end-tt)');
set(gca,'YDir','normal');
hold on;
for i=1:n:nx
    plot([X(i,1,1),X(i,1,1)],[min(Y(:)),max(Y(:))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

for i=1:n:ny
    plot([min(X(:)),max(X(:))],[(Y(1,i,1)),(Y(1,i,1))],'color','black');
    hold on;
end
hold on;
plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');

subplot(2,2,3)
imagesc([min(X(:)),max(X(:))],[min(Y(:)),max(Y(:))],s2(:,:,k)')
set(gca,'YDir','normal');
%% evaluate inversion
%{
s(:,:,end-1:-1:1)=s(:,:,1:end-1);
s2(:,:,end-1:-1:1)=s2(:,:,1:end-1);
tt2=vp;
for i=1:floor(nx/n)+1
    for j=1:floor(ny/n)+1
        for k=1:floor(nz/n)
            if s(i,j,k)==0
                tt=nan;
            else
                tt=1;
            end
            if i==floor(nx/n)+1 || j==floor(ny/n)+1 || k==floor(nz/n)
                tt2(((1+(i-1)*n):end), ((1+(j-1)*n):end), ((1+(k-1)*n):end))=tt*tt2(((1+(i-1)*n):end), ((1+(j-1)*n):end), ((1+(k-1)*n):end));
            else
                tt2(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n))=tt*tt2(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n));
            end
        end
    end
end
vp=tt2;

tt2=vs;
for i=1:floor(nx/n)+1
    for j=1:floor(ny/n)+1
        for k=1:floor(nz/n)
            tt=1;
            if s2(i,j,k)==0
                tt=nan;
            end
            if i==floor(nx/n)+1 || j==floor(ny/n)+1 || k==floor(nz/n)
                tt2(((1+(i-1)*n):end), ((1+(j-1)*n):end), ((1+(k-1)*n):end))=tt*tt2(((1+(i-1)*n):end), ((1+(j-1)*n):end), ((1+(k-1)*n):end));
            else
                tt2(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n))=tt*tt2(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n));
            end
        end
    end
end
vs=tt2;
%}
%% evaluate inversion
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

%% vp horizontal xy plane
z=2500:-2000:-26000;
a=ceil(sqrt(length(z)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(z)
    tt=find(abs(Z(1,1,:)-z(i))==min(abs(Z(1,1,:)-z(i))));
    tt=tt(1);
    subplot(a,a,i)
    imAlpha=ones(size(vp(:,:,tt)'));
    imAlpha(isnan(vp(:,:,tt)'))=0;
    imagesc(vp(:,:,tt)'-vp0(:,:,tt)','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp z =' num2str(z(i)) 'm']});
    
    xlabel('X [m]');
    ylabel('Y [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vp_xy'],'-djpeg','-r400');
%% vp vertical yz plane
x=3000:10000:44800;
a=ceil(sqrt(length(x)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(x)
    tt=find(abs(X(:,1,1)-x(i))==min(abs(X(:,1,1)-x(i))));
    tt=tt(1);
    subplot(a,a,i)
    imAlpha=reshape(ones(size(vp(tt,:,:))),[ny,nz])';
    imAlpha(isnan(reshape(vp(tt,:,:),[ny,nz])'))=0;
    imagesc([min(Y(:)),max(Y(:))], ...
        [min(Z(:)),max(Z(:))], ...
        reshape(vp(tt,:,:)-vp0(tt,:,:),[ny,nz])','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp x =' num2str(x(i)) 'm']});
    
    xlabel('Y [m]');
    ylabel('Z [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    %hold on;
    %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vp_yz'],'-djpeg','-r400');
%% vp vertical xz plane
y=3000:10000:60400;
a=ceil(sqrt(length(x)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(y)
    tt=find(abs(Y(1,:,1)-y(i))==min(abs(Y(1,:,1)-y(i))));
    tt=tt(1);
    subplot(a,a,i)
    imAlpha=reshape(ones(size(vp(:,tt,:))),[nx,nz])';
    imAlpha(isnan(reshape(vp(:,tt,:),[nx,nz])'))=0;
    imagesc([min(X(:)),max(X(:))], ...
        [min(Z(:)),max(Z(:))], ...
        reshape(vp(:,tt,:)-vp0(:,tt,:),[nx,nz])','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp y =' num2str(y(i)) 'm']});
    
    xlabel('X [m]');
    ylabel('Z [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    %hold on;
    %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vp_xz'],'-djpeg','-r400');
%% vs horizontal xy plane
z=5000:-2000:-18000;
a=ceil(sqrt(length(z)));

figure;
set(gcf,'position',[80,80,1500,1500]);

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
    
    xlabel('X [m]');
    ylabel('Y [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vs_xy'],'-djpeg','-r400');
%% vs vertical yz plane
x=3000:10000:44800;
a=ceil(sqrt(length(x)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(x)
    tt=find(abs(X(:,1,1)-x(i))==min(abs(X(:,1,1)-x(i))));
    tt=tt(1);
    subplot(a,a,i)
    imAlpha=reshape(ones(size(vs(tt,:,:))),[ny,nz])';
    imAlpha(isnan(reshape(vs(tt,:,:),[ny,nz])'))=0;
    imagesc([min(Y(:)),max(Y(:))], ...
        [min(Z(:)),max(Z(:))], ...
        reshape(vs(tt,:,:)-vs0(tt,:,:),[ny,nz])','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vs x =' num2str(x(i)) 'm']});
    
    xlabel('Y [m]');
    ylabel('Z [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    %hold on;
    %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vs_yz'],'-djpeg','-r400');
%% vs vertical xz plane
y=3000:10000:60400;
a=ceil(sqrt(length(x)));

figure;
set(gcf,'position',[80,80,1500,1500]);

%limps=[1.5,2.1];
for i=1:length(y)
    tt=find(abs(Y(1,:,1)-y(i))==min(abs(Y(1,:,1)-y(i))));
    tt=tt(1);
    subplot(a,a,i)
    imAlpha=reshape(ones(size(vs(:,tt,:))),[nx,nz])';
    imAlpha(isnan(reshape(vs(:,tt,:),[nx,nz])'))=0;
    imagesc([min(X(:)),max(X(:))], ...
        [min(Z(:)),max(Z(:))], ...
        reshape(vs(:,tt,:)-vs0(:,tt,:),[nx,nz])','AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vs y =' num2str(y(i)) 'm']});
    
    xlabel('X [m]');
    ylabel('Z [m]');
    %caxis(limp);
    colormap(redblue2);
    colorbar;
    %hold on;
    %plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
end
print(gcf,['./vs_xz'],'-djpeg','-r400');