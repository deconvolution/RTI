clear all;
close all;
tt=load('./final_model.mat');
data=tt.data;
reference=[35.85,-2.85];
a = [km2deg(double(data.X(:))/1000)+reference(1),...
    km2deg(double(data.Y(:))/1000)+reference(2),...
    double((data.Z(:) + int64(-19000)))/1000,double(data.v(:))/1000];
b = sortrows(a,[1, 2 -3]);
writematrix(b,'Oldoinyo_modVp.txt');