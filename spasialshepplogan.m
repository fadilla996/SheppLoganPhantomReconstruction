%% Program Rekonstruksi Shepp-Logan Phantom 512 Pixel%%
%%Fadilla Sofa Amatullah%%
%%10217012%%

clc
clear all

%%Inisiasi phantom uji

%Menentukan parameter phantom : Jumlah pixel, ukuran pixel, posisi titik 
%pusat, radius dan koefisien atenuasi (gray level)

nx = 512;                       
ny = 512;                       
ps = 1;                         
cx = [0; 0; 56; -56; 0; 0; 0; -20; 0; 15];                  
cy = [0; 3; 0; 0; -90; -26; 26; 154; 154; 154];
osxo = [92; 87; 31; 41; 25; 4; 4; 4; 2; 4];
osx = 2.56*osxo;  
osyo = [69; 66; 11; 16; 21; 4; 4; 2; 2; 2];
osy = 2.56*osyo;
at = [2; -1.5; -0.45; -0.45; 0.45; 0.45; 0.45; 0.45; 0.45; 0.45];                  
deg = [0; 0; -72*pi/180; 72*pi/180; 0; 0; 0; 0; 0; 0];

x = zeros(1,nx);
y = zeros(ny,1);
xx = zeros(nx,ny);
yy = zeros(nx,ny);
phantom = zeros(nx,ny);

x(1,1) = -nx/2;
y(1,1) = ny/2;

for i = 2 : nx
    x(1,i)= x(1,1)+(i-1)*ps;
    y(i,1)= y(1,1)-(i-1)*ps;
end

for i = 1 : nx
    for j = 1 : ny
        xx(i,j) = x(1,j);
        yy(i,j) = y(i,1);
    end
end

% Proses penentuan nilai gray level/ nilai koefisien atenuasi phantom
for i = 1:2
    sx = cx(i,1);
    sy = cy(i,1);
    oscx = osx(i,1);
    oscy = osy(i,1);
    a = at(i,1);
    s = tan(deg(i,1));
    for j=1:nx
        for k=1:ny
            if ((((xx(j,k)-sx)/oscy).^2 + ((yy(j,k)-sy)/oscx).^2) <= 1);
                phantom(j,k)=a+phantom(j,k);
            end
        end
    end
end

for i = 3:4
    sx = cx(i,1);
    sy = cy(i,1);
    oscx = osx(i,1);
    oscy = osy(i,1);
    a = at(i,1);
    s = tan(deg(i,1));
    for j=1:nx
        for k=1:ny
            if ((((((yy(j,k)-sy)-s*(xx(j,k)-sx))^2)/(oscy^2*(1+(s^2))))+((((s*(yy(j,k)-sy))+(xx(j,k)-sx))^2)/(oscx^2*(1+(s^2)))))<= 1);
                phantom(j,k)=a+phantom(j,k);
            end
        end
    end
end

for i = 5:7
    sx = cx(i,1);
    sy = cy(i,1);
    oscx = osx(i,1);
    oscy = osy(i,1);
    a = at(i,1);
    s = tan(deg(i,1));
    for j=1:nx
        for k=1:ny
            if ((((xx(j,k)-sx)/oscy).^2 + ((yy(j,k)-sy)/oscx).^2) <= 1);
                phantom(j,k)=a+phantom(j,k);
            end
        end
    end
end

for i = 8
    sx = cx(i,1);
    sy = cy(i,1);
    oscx = osx(i,1);
    oscy = osy(i,1);
    a = at(i,1);
    s = tan(deg(i,1));
    for j=1:nx
        for k=1:ny
            if ((((xx(j,k)-sx)/oscx).^2 + ((yy(j,k)-sy)/oscy).^2) <= 1);
                phantom(j,k)=a+phantom(j,k);
            end
        end
    end
end

for i = 9:10
    sx = cx(i,1);
    sy = cy(i,1);
    oscx = osx(i,1);
    oscy = osy(i,1);
    a = at(i,1);
    s = tan(deg(i,1));
    for j=1:nx
        for k=1:ny
            if ((((xx(j,k)-sx)/oscy).^2 + ((yy(j,k)-sy)/oscx).^2) <= 1)
                phantom(j,k)=a+phantom(j,k);
            end
        end
    end
end

% Menampilkan gambar phantom
figure(1)
   imagesc(x, y, phantom)               
    colormap('gray')
    title('Shepp-Logan Phantom')
    xlabel('Position')
    ylabel('Position')
    colorbar
%% Akuisisi data proyeksi (sinogram)
% Penentuan jumlah berkas sinar dan jumlah sudut proyeksi

disp(sprintf('Akuisisi data proyeksi...'));

[M N] = size(phantom);
rhomax = round(sqrt(M^2 + N^2));
if mod(rhomax,2)==0;
    nrho = rhomax + 1;
    rho = zeros(1,nrho);
    rho(1,1)=-rhomax/2*ps;
else
    nrho = rhomax + 2;
    rho = zeros(1,nrho);
    rho(1,1)=-(rhomax+1)/2*ps;
end

for i=2:nrho
    rho(1,i)=rho(1,1)+((i-1)*ps);
end

ntheta = 180;
maxtheta = 180;
theta = zeros(ntheta,1);
theta(1,1)=maxtheta/ntheta;
for i = 2:ntheta
    theta(i,1)=theta(1,1)+((i-1)*theta(1,1));
end

res = zeros(nrho, ntheta);
post = zeros(nx,ny);

% Proses penentuan data nilai tiap pixel pada sinogram
for i=1:44
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+phantom(j,k)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end

for i=45
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=((sqrt(ps^2+ps^2))/2);
                    res(l,i)=res(l,i)+(phantom(j,k)*2*abs(post(j,k)-rho(1,l)));
                end
            end
        end
    end
end
for i=46:134
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+phantom(j,k)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end

for i=135
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=((sqrt(ps^2+ps^2))/2);
                    res(l,i)=res(l,i)+(phantom(j,k)*2*abs(post(j,k)-rho(1,l)));
                end
            end
        end
    end
end

for i=136:180
    disp(sprintf('angle %g of %g', i, ntheta));
    for j = 1:nx
        for k = 1:ny
            post(j,k)=xx(j,k)*cos(theta(i,1)*pi/180)+yy(j,k)*sin(theta(i,1)*pi/180);
            for l = 1:nrho
                if abs(post(j,k)-rho(1,l))<=(ps/2);
                    res(l,i)=res(l,i)+phantom(j,k)* ps;
                else
                    res(l,i)=res(l,i);
                end
            end
        end
    end
end

% Menampilkan gambar hasil sinogram
figure(2)
imagesc(theta,rho,res);
colormap(gray);
title('Sinogram Shepp-Logan')
xlabel('Angle')
ylabel('Ray Position')
colorbar

figure(3)
    plot(rho, res(:,45))
    title('Projection data at 45 degrees');
    xlabel('Ray Position');
    ylabel('P_\theta (t)');
%% Rekonstruksi Citra dengan Proyeksi Balik Tanpa Filter

disp(sprintf('\nProyeksi Balik Citra...'));
% Proyeksi balik data sinogram
sinogram = res(108:619,:);
lamin = zeros(nx,ny);  
  for ia = 1:180
    projection_ia=sinogram(:,ia); 
    projection_smear=repmat (projection_ia,1,512);
    rot= imrotate(projection_smear', ia, 'bicubic','crop');
    lamin=lamin+(rot'/180);   
  end

%Kalibrasi nilai gray level
m = max(max(lamin))+abs(min(min(lamin)));
lamin = ((lamin + abs(min(min(lamin))))/m)*2;

% Menampilkan citra hasil proyeksi balik
figure(4)
imagesc(x, y, lamin'); colormap('gray'); 
axis('image');
title('Simple Backprojection Image');
xlabel('mm');  
ylabel('mm');
colorbar;

%Hitung Nilai PSNR
input = phantom;
output = lamin';
M = 512;
N = 512;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Ram-Lak

disp(sprintf('\nFilter Ram-Lak...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Ram-Lak
sinogramfiltered=sinogram;     
a = length(sinogram);
freq=linspace(-1, 1, a).';
Filter = abs(freq);
Filtersp = ifftshift(ifft(ifftshift(Filter)));
Filtersp = real(Filtersp);
sinogramfilt = zeros(512,180);
for i = 1:180
    sinogramfilt(:,i)=imfilter(sinogramfiltered(:,i),Filtersp,'same','conv');    
end  

%Proyeksi balik data sinogram
bpf_recon = zeros(nx,ny);
  for ia = 1:180
    bpf_ia=sinogramfilt(:,ia);
    bpf_smear=repmat(bpf_ia,1,512);
    rot1= imrotate(bpf_smear', ia, 'bicubic','crop');   
    bpf_recon=bpf_recon+(rot1'/(180));
  end
  
%Kalibrasi nilai gray level
m = max(max(bpf_recon))+abs(min(min(bpf_recon)));
bpf_recon = ((bpf_recon + abs(min(min(bpf_recon))))/m)*2;

figure(5)
imagesc(theta,rho, sinogramfilt)   
colormap('gray')                 
title('Sinogram Filtered Ram-Lak')
xlabel('Angle')
ylabel('Ray Position')
colorbar

% Plot Frekuensi Filter Ram-Lak
figure(6)
plot(freq, Filtersp);  
title('Filter Ram-Lak')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Ram-Lak
figure(7)
imagesc(x, y, bpf_recon'); colormap('gray'); 
axis('image')  
title('Filter Ram-Lak')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon';
M = 512;
N = 512;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Hamming  

disp(sprintf('\nFilter Hamming...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hamming
sinogramfiltered2=sinogram;
b = length(sinogram);
hamm = zeros(b,1);
freq=linspace(-1, 1, b).';
absfreq = abs(freq);
for i=1:b
    hamm(i,1)=0.54+0.46*cos(pi*freq(i,1));
end
hamm = hamm.*absfreq; 
hammsp = real(ifftshift(ifft(ifftshift(hamm))));
sinogramfit2 = zeros(512,180);
for i = 1:180
    sinogramfilt2(:,i)=imfilter(sinogramfiltered2(:,i),hammsp,'same','conv');    
end

% Proyeksi balik data sinogram
bpf_recon2 = zeros(nx,ny);
  for ia = 1:180
    bpf_ia2=sinogramfilt2(:,ia);
    bpf_smear2=repmat(bpf_ia2,1,512);
    rot2= imrotate(bpf_smear2', ia, 'bicubic','crop');   
    bpf_recon2=bpf_recon2+(rot2'/(180));
  end
  
%Kalibrasi nilai gray level
m = max(max(bpf_recon2))+abs(min(min(bpf_recon2)));
bpf_recon2 = ((bpf_recon2 + abs(min(min(bpf_recon2))))/m)*2;

% Respon Frekuensi Filter Hamming
figure(8)
plot(freq, hammsp);  
title('Filter Hamming')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Hamming
figure(9)
imagesc(x, y, bpf_recon2'); colormap('gray'); 
axis('image')  
title('Filter Hamming')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon2';
M = 512;
N = 512;
peakval = 2;
Fungsi(input,output,M,N,peakval);
%% Rekonstruksi Citra dengan Penerapan Filter Hanning

disp(sprintf('\nFilter Hanning...'));
% Perkalian sinogram dalam domain frekuensi dengan nilai frekuensi pada
% filter Hanning
sinogramfiltered3=sinogram;
c = length(sinogram);
hann = zeros(c,1);
freq=linspace(-1, 1, c).';
absfreq = abs(freq);
for i=1:c
    hann(i,1)=0.5*(1+cos(pi*freq(i,1)));
end
hann = hann.*absfreq;
hannsp = real(ifftshift(ifft(ifftshift(hann))));
sinogramfilt3 = zeros(512,180);
for i = 1:180
    sinogramfilt3=imfilter(sinogramfiltered3,hannsp,'same','conv');  
end   

%Proyeksi balik data sinogram
bpf_recon3 = zeros(nx,ny);
  for ia = 1:180
    bpf_ia3=sinogramfilt3(:,ia);
    bpf_smear3=repmat(bpf_ia3,1,512);
    rot3= imrotate(bpf_smear3', ia, 'bicubic','crop');   
    bpf_recon3=bpf_recon3+(rot3'/(180));
  end

%Kalibrasi nilai gray level
m = max(max(bpf_recon3))+abs(min(min(bpf_recon3)));
bpf_recon3 = ((bpf_recon3 + abs(min(min(bpf_recon3))))/m)*2;

% Respon Frekuensi Filter Hanning
figure(10)
plot(freq, hannsp);  
title('Filter Hanning')
xlabel('Freq(w)')
ylabel('H(w)')

% Menampilkan rekonstruksi citra dengan penerapan Filter Hanning
figure(11)
imagesc(x, y, bpf_recon3'); colormap('gray'); 
axis('image')  
title('Filter Hanning')
xlabel('Position')
ylabel('Position')
colorbar

%Hitung Nilai PSNR
input = phantom;
output = bpf_recon3';
M = 512;
N = 512;
peakval = 2;
Fungsi(input,output,M,N,peakval);

function [MSE,PSNR] = Fungsi(input,output,M,N,peakval)
    MSE = sum(sum((input-output).^2))/(M*N);
    PSNR = 10*log10((peakval^2)/MSE);
    fprintf('MSE : %7.2f ', MSE);
    fprintf('\nPSNR : %9.7f dB \n', PSNR);
end