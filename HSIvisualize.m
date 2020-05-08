function HSIvisualize(M,RGB,d,wl)
% Function to visualize the HSI cube based on :
        % * Raw, Continuum removed, Continuum and First derivative  
        %       spectral profiles along the sample.
        % * Minimum Noise Fraction compression to create a composite image
        % * Wavelength correlations.
        
m=median(M(:));
if m>1000
    M=M/10000;
end

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

end