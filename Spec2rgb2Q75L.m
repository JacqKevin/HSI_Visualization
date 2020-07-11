function [RGBim,R700_450, L]=Spec2rgb2Q75L(M,d,wl)
% Function that convert hyperspectral image into RGB, and to L*a*b* to compute a equivalent to the Q7/4 diagram.
% INPUT:
%           M: Hyperspectral datacube (n*m*p)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
% OUTPUT:
%           RGBim: RGB image
%           R700_450: ratio of 700 nm and 450 nm
%           L: Luminance
% Debret, M., Sebag, D., Desmet, M., Balsam, W., Copard, Y., Mourier, B., Susperrigui, A.-S., Arnaud, F., Bentaleb, I., Chapron, E., Lallier-Vergès, E., Winiarski, T., 2011. Spectrocolorimetric interpretation of sedimentary dynamics: The new “Q7/4 diagram.�? Earth-Science Rev. 109, 1–19. https://doi.org/10.1016/j.earscirev.2011.07.002

if size(M,3)>1
    S=reshape(M,[],size(M,3));
end
if median(S(:))>1000
    S=S/10000;
end

[x, y, z, ~, ~] = Get_xyz();
close all

nmi=x(1,1):0.01:x(end,1);
xi=interp1(x(:,1)',x(:,2)',nmi);
yi=interp1(y(:,1)',y(:,2)',nmi);
zi=interp1(z(:,1)',z(:,2)',nmi);

wlr=round(wl,2);

xs=zeros(1,length(wlr));
ys=zeros(1,length(wlr));
zs=zeros(1,length(wlr));
for i=1:length(wlr)
    [~,b]=find(abs(wlr(i)-nmi)==min(abs(wlr(i)-nmi)));
    if abs(wlr(i)-nmi(b))<0.1
        xs(i)=xi(b);
        ys(i)=yi(b);
        zs(i)=zi(b);
    else
        xs(i)=0;
        ys(i)=0;
        zs(i)=0;
    end
end

X=S*xs';
Y=S*ys';
Z=S*zs';

XYZ=[X/sum(xs) Y/sum(ys) Z/sum(zs)];

LAB=xyz2lab(XYZ);
RGB=xyz2rgb(XYZ);

wl450=find(abs(wlr-450)==min(abs(wlr-450)));
wl700=find(abs(wlr-700)==min(abs(wlr-700)));

L=reshape(LAB(:,1),size(M,1),size(M,2));
R700_450=reshape(S(:,wl700)./S(:,wl450),size(M,1),size(M,2));
RGBim=reshape(RGB,size(M,1),size(M,2),3);

figure;
plot(LAB(:,1),S(:,wl700)./S(:,wl450),'.')
grid on
xlabel('L* (%)')
ylabel('Q700/450')
title('Q700/450 vs L* (global)')
set(gca,'fontsize',14)

figure;
plot(L(round(size(L,1)/2),:),R700_450(round(size(L,1)/2),:),'.')
grid on
xlabel('L* (%)')
ylabel('Q700/450')
title('Q700/450 vs L* (central)')
set(gca,'fontsize',14)

figure
ha(1)=subplot(311);
imagesc(d,d(1:size(RGBim,1)),RGBim)
colorbar
xlabel('Depth (cm)')
ylabel('Width (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(312);
imagesc(d,d(1:size(RGBim,1)),R700_450)
colorbar
% caxis([nanmedian(R700_450(:))-3*nanstd(R700_450(:)) ...
%     nanmedian(R700_450(:))+3*nanstd(R700_450(:))])
title('700/450 ratio')
xlabel('Depth (cm)')
ylabel('Width (cm)')
set(gca,'fontsize',14)
ha(3)=subplot(313);
imagesc(d,d(1:size(RGBim,1)),L)
caxis([max([0 nanmedian(R700_450(:))-3*nanstd(R700_450(:))])...
    min([nanmedian(R700_450(:))+3*nanstd(R700_450(:)) 100])])
colormap(jet)
colorbar
title('L* (%)')
xlabel('Depth (cm)')
ylabel('Width (cm)')
set(gca,'fontsize',14)
linkaxes(ha,'xy')
end