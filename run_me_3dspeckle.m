n=2^10*1.5;
nz=n/1.5;
% NA=0.1;
% pupil_over_focallength=asin(tan(NA))*2;
% max_w0=(length_x/2).^2+(length_y/2).^2; % maximum pupil allow this grid
length_x=1; length_y=1; % pysical length in mm 
w0=1; % the standard deviation of the guassian beam size in mm
dx=length_x/n;
scatter_angle=1/180*pi; %% determined by the material
focal_length=4; % in mm > length_x,length_y
% NA=(w0*pupil)/(focal_length)
pupil=0.7; % unitless
NA=sin(atan(pupil*length_x/2/focal_length))
wav_length=0.0005;
assert(wav_length>=2*dx*NA);% wave length (should satisfy nyquist wave length/NA>2*dx)

zmin=0;zmax=2;
z_pos=linspace(zmin,zmax,nz);
z=z_pos.*focal_length;
x=linspace(0,length_x,n);
y=linspace(0,length_y,n);
max_angle=((wav_length/dx))/2; % maximum angle in Fourier grid (determined by x grid)
flag_diffuser=false;
slm_phase=false;
%
[XX,YY]=meshgrid(x,y); % creating the mesh
XX=XX-mean(XX(:)); 
YY=YY-mean(YY(:));
R=sqrt(XX.^2+YY.^2);

% Gaussian Beam: 
E0=exp(-(R/w0).^2);
% Diffuser: 
if flag_diffuser
T_diff=exp(1i*2*pi*rand(size(XX)));
T_K=fftshift(fft2(T_diff));
T_K=T_K.*(((R/max(XX(:)))*max_angle)<scatter_angle);  % determine scattering angle

T_diff=ifft2(fftshift(T_K));
T_diff=exp(1i*angle(T_diff)); % convert to thin phase plate (meaning takeing only the phase under fft and not the intensity)
if slm_phase
numerical_img=zeros(size(E0));
slm_gray_to_pi=200;
numerical_img(n/2-300:n/2+300-1,n/2-400:n/2+400-1)=exp(1i*2*pi*phase3/slm_gray_to_pi);
T_diff=numerical_img;
% zoomed_slm_sim_speckle3D=slm_sim_speckle3D(768-268:768+269,768-268:768+269,512-200:512+200);
end
E_in=E0.*T_diff; % un/comment for pointscan/spcekle
else    
E_in=E0; % un/comment for pointscan/spcekle
end
%figure;imagesc(abs(E_in))
% Lens Phase
phase_of_lens=focal_length.*ones(size(XX))-sqrt(focal_length^2.*ones(size(XX))-(XX.^2)-(YY.^2));
phase_of_lens=phase_of_lens.*(2*pi/wav_length);
A_lens=exp(-1i.*phase_of_lens);
A_lens=A_lens .* (R<(max(XX(:))*(pupil)));
E_until_lens=E_in.*A_lens;

figure;
subplot(2,2,1);imagesc(abs(E_until_lens))
subplot(2,2,2);imagesc(angle(E_until_lens))
subplot(2,2,3);imagesc(abs(fftshift(fft2(E_until_lens))));
subplot(2,2,4);imagesc(abs(volume_angular_spectrum_prop_absorbtion(E_until_lens,x,y,focal_length,wav_length,normalized)));

%%
normalized=false;
recursive=true;
absorber.put_before_focal=false; %change with/without absorbers, 
absorber.size.x=round(50); %change with/without absorbers, 
absorber.size.y=round(10); %change with/without absorbers, 
absorber.size.z=round(1); %change with/without absorbers, 
absorber.position.x=round(n/2); %change with/without absorbers, 
absorber.position.y=round(n/2); %change with/without absorbers, 
absorber.position.z=floor(nz/3);
% change absobers size in ***volume_angular_spectrum_prop.m***

PSF3D=volume_angular_spectrum_prop_absorbtion(E_until_lens,x,y,z,wav_length,normalized,recursive,absorber);
figure;imagesc(PSF3D(:,:,nz/2))
figure;imagesc(squeeze(PSF3D(:,n/2,:))')
f1=figure('Units','normalized','Position',[0.35 0.5642 0.5 0.3500]);
subplot(1,2,1);imagesc((PSF3D(:,:,nz/2)).^1); colormap(viridis);
set(gca,'XTick',[]);set(gca,'YTick',[]);title('XY','Fontsize',32)
subplot(1,2,2);imagesc(squeeze((PSF3D(:,n/2,:)).^1)'); colormap(viridis);
set(gca,'YTick',[]);set(gca,'XTick',[]);title('XZ','Fontsize',32);
% %% vertical prefernce
% A=(speckle3D(:,:,nz/2));
% B=circshift(A,[0,-50]);
% figure;imagesc(B+A);
% 
% f2=figure('Units','normalized','Position',[0.35 0.2 0.18 0.7]);
% hAxis(1)=subplot(2,1,1);imagesc(A.^0.7);
% set(gca,'XTick',[]);set(gca,'YTick',[]);title('XY','Fontsize',32)
% hAxis(2)=subplot(2,1,2);imagesc(B.^3);
% set(gca,'YTick',[]);set(gca,'XTick',[]);title('XZ','Fontsize',32);
% 
% 
if flag_diffuser
    speckle3D=PSF3D;
else
    point3D=PSF3D;
end
clear PSF3D

%%

%  figure;imagesc(angle(T_diff));
% colormap(viridis); axis image;
% set(gca,'XTick',[]);set(gca,'YTick',[]);
% %export_fig(sprintf('E0%i.png',1),'-transparent','-png',gcf)
% cbh = colorbar('v');
% %get(get(cbh, 'Children'))
% set(cbh, 'Limits', [-pi-eps,pi+eps])
% set(cbh, 'YTick', [-pi,0,pi])
% set(cbh,  'YTickLabel',{'-\pi','0','+\pi'})
% set(cbh,  'Fontsize',18)
