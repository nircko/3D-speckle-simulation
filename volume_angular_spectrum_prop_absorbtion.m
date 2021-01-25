function Intensity3d = volume_angular_spectrum_prop_absorbtion(E_in,x,y,z_prop,wl,normalized,recursive,change_in_propagation)

% DESCRIPTION:
% volume_angular_spectrum_prop  is an extension to 
% angular_spectrum_prop from Fourier Optic by Goddman
%
% USAGE:
%     speckles3d = volume_angular_spectrum_prop(E_in,x,y,3e2,632e-6,true)
%     speckles3d = volume_angular_spectrum_prop(E_in,x,y,[0:0.01:1].*3e2,632e-6,true)
%
% INPUTS:
%     E_in      - Feild to propagate, notice that it should be the same
%     size of x,y
%     x,y      - 1D arrays in mm units
%     z_prop      - scaler/1D array in mm units of the z,axial location to
%     propagate the feild
%     wl      - scaler, the wavelength in mm
%    normalized - bool. whether to normazlie each axial section of the 2D
%    feild
%    recursive - bool. whether to iterate recursively


%
% OUTPUTS:
%     Intensity3d            - 3D array of Intensity of the EM field.

if ~exist('normalized','var') || nargin<6
    normalized=false;
end
    

if ~exist('recursive','var') || nargin<7
    recursive=false;
end
   

if ~exist('change_in_propagation','var') || nargin<8
    change_in_propagation=false;
end


%% volume_angular_spectrum_prop return 3d array of propgating scaler EM feild
% using angular spectrum method
    s_z=length(z_prop);
    [kx,ky,s_x,s_y]=transverse_fourier_coordinates(x,y);
    alpha=(kx/(2*pi)).*wl;  % unitless normalized fourier coordinates
    beta=(ky/(2*pi)).*wl;  % unitless normalized fourier coordinates
    [Alpha, Beta]= meshgrid(alpha,beta);
    mu = ((2*pi)/wl).* sqrt((Alpha.^2+Beta.^2)-1);
    convolved_fourier=fftshift(fft2(E_in));
    Intensity3d=zeros(s_x,s_y,s_z);
    %angular spectrum of propogation distance z
    if recursive
    z_pos=z_prop(2)-z_prop(1); % step size
    angs_prop_z=exp(1.*mu.*z_pos);
    E_out=ifft2(ifftshift(convolved_fourier.*exp(1.*mu.*(z_prop(1)-z_pos))));
    for z_counter=1:s_z
      E_in=E_out;
        if change_in_propagation.put_before_focal
            E_before_absorber=E_in;
            E_in=put_absorber(E_before_absorber,z_counter,change_in_propagation);
        end
      convolved_fourier=fftshift(fft2(E_in));
      E_out=ifft2(ifftshift(convolved_fourier.*angs_prop_z));
      Iout=abs(E_out).^2;
      if normalized
      Iout=Iout./max(Iout(:)); % Normalized
      end
      Intensity3d(:,:,z_counter)=Iout;
    end
    
    else
    
    for z_counter=1:s_z
      z_pos=z_prop(z_counter);
      angs_prop_z=exp(1.*mu.*z_pos);
      E_out=ifft2(ifftshift(convolved_fourier.*angs_prop_z));
      Iout=abs(E_out).^2;
      if normalized
      Iout=Iout./max(Iout(:)); % Normalized
      end
      Intensity3d(:,:,z_counter)=Iout;
    end
    end
    
end
    
      function [kx,ky,s_x,s_y]=transverse_fourier_coordinates(x,y)
    s_x=length(x);
    s_y=length(y);
    temp1=[-s_x/2:s_x/2-1]/s_x; % [-0.5 0.5] unitless
    temp2=[-s_y/2:s_y/2-1]/s_y; % [-0.5 0.5] unitless
    kx=temp1.*(s_x*2*pi/(max(x)-min(x))); % fourier coordinates in 2pi/mm
    ky=temp2.*(s_y*2*pi/(max(y)-min(y)));
      end
      
      function E_out=put_absorber(E_in,z_counter,absorber)
      E_out=E_in;
        if z_counter==floor(absorber.position.z)
           size_x=([-absorber.size.x:absorber.size.x]); %absorber size    
           size_y=([-absorber.size.y:absorber.size.y]);   %absorber size
          %absorber in the middle
          absorber_loc_x=size_x+absorber.position.x; %absorber size in grid
          absorber_loc_y=size_y+absorber.position.y; %absorber size in grid
          E_out(absorber_loc_y,absorber_loc_x)=0;
        end
      end