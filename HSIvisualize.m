function HSIvisualize(M,RGB,d,wl)
% Function to visualize the HSI cube based on :
%      * Q7/4 vs L* diagram
%      * Abberant pixels
%      * Median and standard deviation spectra
%      * Raw, Continuum removed, Continuum and First derivative  
%            spectral profiles along the sample.
%      * Minimum Noise Fraction compression to create a composite image
%      * Wavelength correlations.
% INPUT:
%           M: Hyperspectral datacube (n*m*p)
%           RGB: Associated RGB image (n*m*3)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
        
m=median(M(:));
if m>1000
    M=M/10000;
end

% RGB enhancement
RGB=RGB*(0.5/mean(RGB(:)));

% Check unique depth
[~,ia,~] = unique(d);

%% Q700/400 vs L*

if median(wl)<1000
    [RGB,~, ~]=Spec2rgb2Q75L(M,d,wl);
end

%% Aberrant pixels

AbberantPixels(M,RGB,d,wl,1);

%% Median and standard deviation spectra

figure
yyaxis left
plot(wl,median(squeeze(M(round(size(M,1)),:,:))),'linewidth',2)
ylabel('Median spectrum')
yyaxis right
plot(wl,std(squeeze(M(round(size(M,1)),:,:))),'linewidth',2)
ylabel('Standard deviation spectrum')
title('Median and standard deviation spectra')
grid on
set(gca,'fontsize',14)

%% Raw spectra

figure;
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(1,3,2:3);
imagesc(wl,d,squeeze(median(M,1)))
title('Raw spectra')
colormap(flipud(jet))
colorbar
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(ha,'y')
set(gca,'fontsize',14)

figure;
hb(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d(ia),imrotate(RGB(:,ia,:),-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
hb(2)=subplot(1,3,2:3);
surf(wl,d(ia),squeeze(median(M(:,ia,:),1)),'EdgeColor','none')
title('Raw spectra')
colormap(flipud(jet))
colorbar
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(hb,'y')
set(gca,'fontsize',14,'Ydir','reverse')

%% Continuum removed

figure;
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(1,3,2:3);
[Scr,cr,~]=continuum_removal(wl,squeeze(median(M,1)));
imagesc(wl,d,Scr),
title('Continuum removed spectra')
colormap(flipud(jet))
colorbar
caxis([nanmedian(Scr(:))-3*nanstd(Scr(:)) 1])
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(ha,'y')
set(gca,'fontsize',14)

figure;
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(1,3,2:3);
imagesc(wl,d,Scr),
title('Continuum removed spectra')
colormap(flipud(jet))
colorbar
caxis([0.8 1])
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(ha,'y')
set(gca,'fontsize',14)

figure;
hb(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d(ia),imrotate(RGB(:,ia,:),-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
hb(2)=subplot(1,3,2:3);
surf(wl,d(ia),Scr(ia,:),'EdgeColor','none')
title('Continuum removed spectra')
colormap(flipud(jet))
colorbar
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(hb,'y')
set(gca,'fontsize',14,'Ydir','reverse')

%% Continuum 

figure;
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(1,3,2:3);
imagesc(wl,d,cr),
title('Continuum spectra')
colormap(jet)
colorbar
caxis([nanmedian(cr(:))-3*nanstd(cr(:)) nanmedian(cr(:))+3*nanstd(cr(:))])
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(ha,'y')
set(gca,'fontsize',14)

figure;
hb(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d(ia),imrotate(RGB(:,ia,:),-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
hb(2)=subplot(1,3,2:3);
surf(wl,d(ia),cr(ia,:),'EdgeColor','none')
title('Continuum spectra')
colormap(flipud(jet))
colorbar
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(hb,'y')
set(gca,'fontsize',14,'Ydir','reverse')

%% FDS

figure;
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(1,3,2:3);
Sfds=savgol(squeeze(median(M,1)),7,2,1);
imagesc(wl,d,Sfds)
title('First derivative spectra')
colormap(flipud(jet))
colorbar
caxis([nanmedian(Sfds(:))-2*nanstd(Sfds(:)) nanmedian(Sfds(:))+2*nanstd(Sfds(:))])
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(ha,'y')
set(gca,'fontsize',14)

figure;
hb(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d(ia),imrotate(RGB(:,ia,:),-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
hb(2)=subplot(1,3,2:3);
surf(wl,d(ia),Sfds(ia,:),'EdgeColor','none')
title('First derivative spectra')
colormap(flipud(jet))
colorbar
ylabel('Depth (cm)')
xlabel('Wavelength (nm)')
linkaxes(hb,'y')
set(gca,'fontsize',14,'Ydir','reverse')

%% MNF 

[mnf, A, noiseFractions] = hyperMnf(reshape(M,[],size(M,3))', size(M,1), size(M,2));
mnf=reshape(mnf',size(M,1),size(M,2),size(M,3));
mnf=mnf(:,:,end-2:end);
mnf_min=squeeze(min(min(mnf)));
mnf_max=squeeze(max(max(mnf)));
mnf_min_m=ones(size(M,1),size(M,2),3);
mnf_max_m=ones(size(M,1),size(M,2),3);
for i=1:3
    mnf_min_m(:,:,i)=squeeze(mnf_min_m(:,:,i))*mnf_min(i);
    mnf_max_m(:,:,i)=squeeze(mnf_max_m(:,:,i))*mnf_max(i);
end
MNF=(mnf-mnf_min_m)./(mnf_max_m-mnf_min_m);

figure;
ha(1)=subplot(3,3,[1 4 7]);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(3,3,[2 5 8]);
imagesc(d(1:size(RGB,1)),d,imrotate(MNF,-90))
title('MNF composite image')
colormap(jet)
ylabel('Depth (cm)')
xlabel('Width (cm)')
set(gca,'fontsize',14)
subplot(3,3,3);
plot(wl,A(:,end),'r','linewidth',2)
grid on
xlabel('Wavelength (nm)')
ylabel('Reflectance')
xlim([wl(1) wl(end)])
title(['MNF1 ',num2str(noiseFractions(end))])
set(gca,'fontsize',14)
subplot(3,3,6);
plot(wl,A(:,end-1),'g','linewidth',2)
grid on
xlabel('Wavelength (nm)')
ylabel('Reflectance')
xlim([wl(1) wl(end)])
title(['MNF2 ',num2str(noiseFractions(end-1))])
set(gca,'fontsize',14)
subplot(3,3,9);
plot(wl,A(:,end-2),'b','linewidth',2)
grid on
xlabel('Wavelength (nm)')
ylabel('Reflectance')
xlim([wl(1) wl(end)])
title(['MNF3 ',num2str(noiseFractions(end-2))])
set(gca,'fontsize',14)
linkaxes(ha,'y')

%% Wavelength correlation

figure
imagesc(wl,wl,corr(reshape(M,[],size(M,3))))
colormap(jet)
colorbar
ylabel('Wavelength (nm)')
xlabel('Wavelength (nm)')
set(gca,'fontsize',14)

%% Clustering

% S
evaS = evalclusters(squeeze(median(M,1)),'kmeans','CalinskiHarabasz','KList',[1:10]);
% Scr
evaScr = evalclusters(Scr,'kmeans','CalinskiHarabasz','KList',[1:10]);
% Cr
evaCr = evalclusters(cr,'kmeans','CalinskiHarabasz','KList',[1:10]);
% Sfds
evaSfds = evalclusters(Sfds,'kmeans','CalinskiHarabasz','KList',[1:10]);

Km=[evaS.OptimalY evaScr.OptimalY evaCr.OptimalY evaSfds.OptimalY];
Kmxlab={strcat('S (',num2str(evaS.OptimalK),')'),...
    strcat('Scr (',num2str(evaScr.OptimalK),')'),...
    strcat('Cr (',num2str(evaCr.OptimalK),')'),...
    strcat('Sfds (',num2str(evaSfds.OptimalK),')')};

% S
evaS = evalclusters(squeeze(median(M,1)),'linkage','CalinskiHarabasz','KList',[1:10]);
% Scr
evaScr = evalclusters(Scr,'linkage','CalinskiHarabasz','KList',[1:10]);
% Cr
evaCr = evalclusters(cr,'linkage','CalinskiHarabasz','KList',[1:10]);
% Sfds
evaSfds = evalclusters(Sfds,'linkage','CalinskiHarabasz','KList',[1:10]);

HACm=[evaS.OptimalY evaScr.OptimalY evaCr.OptimalY evaSfds.OptimalY];
HACmxlab={strcat('S (',num2str(evaS.OptimalK),')'),...
    strcat('Scr (',num2str(evaScr.OptimalK),')'),...
    strcat('Cr (',num2str(evaCr.OptimalK),')'),...
    strcat('Sfds (',num2str(evaSfds.OptimalK),')')};

figure
ha(1)=subplot(131);
imagesc(d(1:size(RGB,1)),d,imrotate(RGB,-90))
xlabel('Width (cm)')
ylabel('Depth (cm)')
set(gca,'fontsize',14)
ha(2)=subplot(132);
imagesc(1:4,d,Km)
colormap(jet)
colorbar
xlabel('Width (cm)')
ylabel('Data cluster')
title('Kmeans')
set(gca,'fontsize',14,'xtick',1:4,'xticklabel',Kmxlab)
ha(3)=subplot(133);
imagesc(1:4,d,HACm)
colormap(jet)
colorbar
xlabel('Width (cm)')
ylabel('Data cluster')
title('Agglomerative Hierarchical Tree')
set(gca,'fontsize',14,'xtick',1:4,'xticklabel',HACmxlab)
linkaxes(ha,'y')
end