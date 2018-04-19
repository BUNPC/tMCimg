%%
% Parameters
filenm1 = 'test2b.2pt';
filenm2 = 'test3.2pt';

iFx = 51;
iFy = 51;
iFz = 4;
iF2z = 53;

nPhotons = 1e7;
mua = 1; % 1/mm (100uM = 0.02 for 800nm, 1/mm for 570nm)
h = 0.01;
hVol = h^3;

nx = 100;
ny = 100;
nz = 100;
nA3 = 4;
nA1 = 8;

%%
fid = fopen(filenm1);
Io = fread(fid,'float32');
IoOut = Io((end-1):end);
Io = reshape(Io(1:end-2),[nx ny nz nA3 nA1]);
fclose(fid);

% normalize the 2pt
lstIn = find(Io>0);
lstOut = find(Io<0);
Io(lstOut) = Io(lstOut) / nPhotons;
IoOut(2) = IoOut(2) / nPhotons;
nPhotonDet = -(sum(Io(lstOut))+IoOut(2))
k = (1-nPhotonDet) / (mua * hVol * (sum(Io(lstIn))+IoOut(1)))
Io(lstIn) = Io(lstIn) * k;
IoOut(1) = IoOut(1) * k;

clear lstIn lstOut

% plot
I = sum(sum(Io,5),4);

figure(1);
imagesc(log10(squeeze(max(I(:,iFy,:),1e-10))),[0 4])
colormap(jet)

%%
% load second 2pt
fid = fopen(filenm2);
I2o = fread(fid,'float32');
I2oOut = I2o((end-1):end);
I2o = reshape(I2o(1:end-2),[nx ny nz nA3 nA1]);
fclose(fid);

% normalize the 2pt
lstIn = find(I2o>0);
lstOut = find(I2o<0);
I2o(lstOut) = I2o(lstOut) / nPhotons;
I2oOut(2) = I2oOut(2) / nPhotons;
nPhotonDet2 = -(sum(I2o(lstOut))+I2oOut(2))
k2 = (1-nPhotonDet2) / (mua * hVol * (sum(I2o(lstIn))+I2oOut(1)))
I2o(lstIn) = I2o(lstIn) * k2;
IoOut(1) = IoOut(1) * k2;

clear lstIn lstOut


% plot
I2 = sum(sum(I2o,5),4);


figure(2);
imagesc(log10(squeeze(max(I2(:,iFy,:),1e-10))),[0 4])
colormap(jet)



%%
% calculate plane wave 2pt from Io 
%absorption
Ip = zeros(size(Io,1),size(Io,3),size(Io,4),size(Io,5));
for iy = 1:size(Io,2)
    Isub = squeeze(Io(:,iy,:,:,:));
    for ii=1:(size(Io,1)-1)
        lst = [ii:size(Isub,1)];
        Ip(1:length(lst),:,:,:) = Ip(1:length(lst),:,:,:) + Isub(lst,:,:,:);
    end
    for ii=2:(size(I,1)-1)
        lst = [1:(size(Isub,1)-ii+1)];
        Ip(ii:end,:,:,:) = Ip(ii:end,:,:,:) + Isub(lst,:,:,:);
    end
end
Ip = Ip / (size(Io,1)*size(Io,2));

%%
% calculate 3pts
Isub = squeeze(Io(:,iFy,:,:,:));
I2sub = squeeze(I2o(:,iFy,:,:,:));

Isub = Isub(:,:,end:-1:1,:); % reverse direction along elevation and then azimuth
I2sub = I2sub(:,:,end:-1:1,:); % reverse direction along elevation and then azimuth
nAzi = size(Io,5); % assumed to be even
lst = [nAzi/2+[1:nAzi/2] 1:nAzi/2];
Isub = Isub(:,:,:,lst);
I2sub = I2sub(:,:,:,lst);

% absorption 3pt
Iap3pt = sum(sum(Ip .* Isub,4),3) * hVol / sum(sum(Ip(iFx,iFz,:,:),4),3);
I2ap3pt = sum(sum(Ip .* I2sub,4),3) * hVol / sum(sum(Ip(iFx,iF2z,:,:),4),3);

%fluorescence 3pt
Ifp = sum(sum(Ip,4),3);
Ifsub = sum(sum(Isub,4),3);
I2fsub = sum(sum(I2sub,4),3);
Ifp3pt = Ifp .* Ifsub * hVol / Ifp(iFx,iFz);
I2fp3pt = Ifp .* I2fsub * hVol / Ifp(iFx,iF2z);

%%
% plot the 3pts
figure(3);
ha(1)=subplot(2,4,1);
foo = log10(squeeze(max(Ifp,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('Source')
axis image

ha(2)=subplot(2,4,2);
foo = log10(squeeze(max(Ifsub,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('Detector Pixel')
axis image

ha(3)=subplot(2,4,3);
foo = log10(squeeze(max(I2fsub,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('Detector Pixel')
axis image

ha(4)=subplot(2,4,5);
foo = log10(squeeze(max(Ifp3pt,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('3pt Fluorescence')
axis image

ha(5)=subplot(2,4,6);
foo = log10(squeeze(max(I2fp3pt,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('3pt Fluorescence')
axis image

ha(6)=subplot(2,4,7);
foo = log10(squeeze(max(Iap3pt,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('3pt Absorption')
axis image

ha(7)=subplot(2,4,8);
foo = log10(squeeze(max(I2ap3pt,1e-10)));
imagesc(foo,[-3 0]+max(foo(:)))
title('3pt Absorption')
axis image

linkaxes(ha,'xy')

colormap(jet)


%%
% plot the lateral profiles of the 3pts at different depths

iZ = [20:20:100];

ha = [];

figure(4)
ha(1)=subplot(1,2,1);
h1=semilogy( Ifp3pt(:,iZ) );
hold on
h2=semilogy( I2fp3pt(:,iZ), '--' );
hold off
title('3pt Fluorescence')
for ii=1:length(h1)
    set(h2(ii),'color',get(h1(ii),'color'))
end

ha(2)=subplot(1,2,2);
h1=semilogy( Iap3pt(:,iZ) );
hold on
h2=semilogy( I2ap3pt(:,iZ), '--' );
hold off
for ii=1:length(h1)
    set(h2(ii),'color',get(h1(ii),'color'))
end
title('3pt Absorption')
legend('20','40','60','80','100');

linkaxes(ha,'x')


%%
% compare abs and flr lateral profiles at different depths
iZ = 20;

figure(5)
semilogy( [Ifp3pt(:,iZ)/max(Ifp3pt(:,iZ)) Iap3pt(:,iZ)/max(Iap3pt(:,iZ))] )
legend( 'fluorescence', 'absorption')
title( sprintf('Depth = %d',iZ) )

%%
% calculate depth sensitivity for a uniform slab change in abs
% need to recalculate 3pts for all y indices

IA3pt = zeros(size(I));
I2A3pt = zeros(size(I));

IF3pt = zeros(size(I));
I2F3pt = zeros(size(I));

for iy=1:ny
    
    % calculate 3pts
    Isub = squeeze(Io(:,iy,:,:,:));
    I2sub = squeeze(I2o(:,iy,:,:,:));
    
    Isub = Isub(:,:,end:-1:1,:); % reverse direction along elevation and then azimuth
    I2sub = I2sub(:,:,end:-1:1,:); % reverse direction along elevation and then azimuth
    nAzi = size(Io,5); % assumed to be even
    lst = [nAzi/2+[1:nAzi/2] 1:nAzi/2];
    Isub = Isub(:,:,:,lst);
    I2sub = I2sub(:,:,:,lst);
    
    % absorption 3pt
    IA3pt(:,iy,:) = permute( sum(sum(Ip .* Isub,4),3) * hVol / sum(sum(Ip(iFx,iFz,:,:),4),3), [1 3 2]);
    I2A3pt(:,iy,:) = permute( sum(sum(Ip .* I2sub,4),3) * hVol / sum(sum(Ip(iFx,iF2z,:,:),4),3), [1 3 2]);
    
    %fluorescence 3pt
    Ifp = sum(sum(Ip,4),3);
    Ifsub = sum(sum(Isub,4),3);
    I2fsub = sum(sum(I2sub,4),3);
    IF3pt(:,iy,:) = permute( Ifp .* Ifsub * hVol / Ifp(iFx,iFz), [1 3 2]);
    I2F3pt(:,iy,:) = permute( Ifp .* I2fsub * hVol / Ifp(iFx,iF2z), [1 3 2]);

end

%%
% plot depth sensitivity for a uniform slab change in abs

dIA = squeeze(sum(sum(IA3pt,2),1));
dI2A = squeeze(sum(sum(I2A3pt,2),1));

dIF = squeeze(sum(sum(IF3pt,2),1));
dI2F = squeeze(sum(sum(I2F3pt,2),1));

figure(5)
subplot(1,3,1)
h1=plot([dIF dI2F]);
title('fluorescence')

subplot(1,3,2)
h2=plot([dIA dI2A],'--');
title('absorption')

subplot(1,3,3)
h1=plot([dIF dI2F]/max(dI2F(:)));
hold on
h2=plot([dIA dI2A]/max(dI2A(:)),'--');
hold off
for ii=1:length(h1)
    set(h2(ii),'color',get(h1(ii),'color'))
end



