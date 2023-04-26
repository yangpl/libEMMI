graphics_toolkit('gnuplot')
set(0, "defaulttextfontname", "Helvetica")
set(0, "defaultaxesfontname", "Helvetica")

clear all

fid = fopen('x1nu');
x1 = fread(fid, 'float');
fclose(fid);

fid = fopen('x2nu');
x2 = fread(fid, 'float');
fclose(fid);

fid = fopen('x3nu');
x3 = fread(fid, 'float');
fclose(fid);

fid = fopen('rho11');
rho = fread(fid, 'float');
fclose(fid);

nx = length(x1);
ny = length(x2);
nz = length(x3);
V=reshape(rho, [nx ny nz]);

iy=46;
figure(1),clf
[xx, zz] = meshgrid(x1, x3);
surf(xx, zz, squeeze(V(:,iy,:)).'); shading interp
set(gca,'Ydir','reverse')
view(0,90)
colormap(jet(2000))
axis tight
xlabel('X[m]')
ylabel('Z[m]')
%ylim([0 3000])
set(gca,'Ydir','reverse')
%set(gca, 'Position', [0.1 0.1 0.8 0.8])
colorbar("NorthOutside")
caxis([0.3 10])

%print -color -deps resistivity.eps
%print -color  resistivity.pdf
%print -djpg resistivity.jpg
print -dpng  resistivity.png

topo = load('topo.txt');
acqui = load('acquisition.txt');

figure(2),clf
plot(xx, zz, 'k')
hold on
plot(xx.', zz.', 'k')
hold on
plot(topo(:,1), topo(:,2), '-')
hold on
plot(acqui(:,1), acqui(:,3),'r*')


set(gca,'Ydir','reverse')
axis tight
xlabel('X[m]')
ylabel('Z[m]')

#print -deps mesh.eps
print -dpng mesh.png



