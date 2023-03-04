clear all;
close all;
tt=load('./m.mat');
vp=tt.data.vp;
X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;
h=X(2)-X(1);

tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_P/final/final_model.mat']);
cip=tt.data.v;
tt=load(['../tutorial_alternating_checkerboard_Crati/' ...
    'inversion_process_lbfgs_S/final/final_model.mat']);
cis=tt.data.v;

tt=load('../tutorial_alternating_checkerboard_Crati/c_vp.mat');
checkerboard_P=tt.data.vp;
tt=load('../tutorial_alternating_checkerboard_Crati/c_vs.mat');
checkerboard_S=tt.data.vs;

turbo2=flip(turbo,1);
redblue2=flip(redblue,1);
%%
receiver=load('receiver_overview.csv');
source=load('source_overview.csv');
%% create checkerboard
[nx,ny,nz]=size(X);

% type the edge length of the block
edge_length=5000;
%
mx=[];
n=0;
while 1
    tt=floor(((.5+n)*edge_length)/h);
    if tt>=nx
        break;
    end
    mx=[mx,tt];
    n=n+1;
end

my=[];
n=0;
while 1
    tt=floor(((.5+n)*edge_length)/h);
    if tt>=ny
        break;
    end
    my=[my,tt];
    n=n+1;
end

mz=[];
n=0;
while 1
    tt=floor(((.5+n)*edge_length)/h);
    if tt>=nz
        break;
    end
    mz=[mz,tt];
    n=n+1;
end
%% vp checkerboard horizontal xy plane
D=mz*h+min(Z(:));
figure;
set(gcf,'position',[80,80,800,1500]);
%%
figure;
imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        checkerboard_P(:,:,2)');
%%
for i=1:length(mz)
    subplot(length(mz),2,2*i-1)
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        checkerboard_P(:,:,mz(i))');
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp checkerboard true model z =' num2str(D(i)) 'm']});
    xlabel('X [m]');
    ylabel('Y [m]');
    caxis([min(checkerboard_P(:)),max(checkerboard_P(:))]);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
    
    subplot(length(mz),2,2*i)
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        cip(:,:,mz(i))');
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vp z =' num2str(mz(i)) 'm']});
    xlabel('X [m]');
    ylabel('Y [m]');
    caxis([min(checkerboard_P(:)),max(checkerboard_P(:))]);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
    title({['vp checkerboard inversion z =' num2str(D(i)) 'm']});
end
print(gcf,['./checkerboard_P_xy'],'-djpeg','-r400');
%% vs checkerboard horizontal xy plane
figure;
set(gcf,'position',[80,80,800,1500]);

for i=1:length(mz)
    subplot(length(mz),2,2*i-1)
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        checkerboard_S(:,:,mz(i))');
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vs checkerboard true model z =' num2str(D(i)) 'm']});
    xlabel('X [m]');
    ylabel('Y [m]');
    caxis([min(checkerboard_S(:)),max(checkerboard_S(:))]);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
    
    subplot(length(mz),2,2*i)
    imagesc([min(X(:)),max(X(:))], ...
        [min(Y(:)),max(Y(:))], ...
        cis(:,:,mz(i))');
    set(gca,'color',1*[1 1 1]);
    set(gca,'ydir','normal');
    title({['vs checkerboard true model z =' num2str(D(i)) 'm']});
    xlabel('X [m]');
    ylabel('Y [m]');
    caxis([min(checkerboard_S(:)),max(checkerboard_S(:))]);
    colormap(redblue2);
    colorbar;
    hold on;
    plot(receiver(:,1),receiver(:,2),'^','MarkerFaceColor','black');
    
end
print(gcf,['./checkerboard_S_xy'],'-djpeg','-r400');