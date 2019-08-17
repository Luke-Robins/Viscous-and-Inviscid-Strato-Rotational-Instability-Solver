function [ fig_vertical, fig_horizontal ] = Visualise( u_k, v_k, w_k, eta_, mu_, m_, k_, N_ )
%Dr Luke Robins 2019 luke.robins@cantab.net
%
%This function takes a result from Inv_SRI_solver or SRI_solver, with
%eigenfunctions already derived.
%The function will produce vertical (r,z) and horizontal (r,theta)
%cross-sections as MatLab plots. These will be saved in the current
%directory.
%   For the vertical cross-section, we hold theta=0.
%   For the horizontal cross-section, we hold z=0.
%
%Input Parameters:
% - [u_k, v_k, w_k]
%   These three terms are truncated Chebyshev expansions of r-dependent
%   complex amplitudes for velocity perturbations in the radial (u),
%   azimuthal (v) and vertical (z) directions (using cylindrical polar
%   co-ordinates).
%
% - eta_
%   This is the radius ratio of the Taylor-Couette system, defined such
%   that eta_ is equal to the inner radius divided by the outer radius.
%
% - mu_
%   This is the angular velocity ratio of the Taylor-Couette system,
%   defined such that mu_ is equal to the angular velocity of the outer
%   cylinder divided by the angular velocity of the inner cylinder.
%
% - m_
%   This is the horizontal (theta) wavenumber of the perturbation. This
%   term must be a real integer.
%
% - k_
%   This is the vertical (z) wavenumber of the system.
%
%Outputs
%   This function will produce two plots - the vertical and horizontal
%   cross-sections of the perturbation.
%       For the vertical cross-section, radial and vertical velocities are
%       shown as vectors, whereas azimuthal velocities are represented as
%       colour contours.
%       For the horizontal cross-section, radial and azimuthal velocities
%       are shown as vectors, whereas vertical velocities are represented
%       as colour contours.
%       The colour contours of the two plots use the same scale, which is
%       shown adjacent to the horizontal cross-section. This colour scale
%       is normalised to the maximum radial velocity perturbation.
% - fig_vertical
%   This is the Matlab-figure-index of the vertical cross-section.
% - fig_horizontal
%   This is the Matlab-figure-index of the horizontal cross-section.
%
%Note:
%   Throughout this file we typically work in the co-ordinates
%   [x_,theta_,z_]. The terms theta_ and z_ both directly correspond to
%   their typical polar co-ordinates.
%   The term x_ ranges from 0 to 1 and denotes co-ordinates across the
%   Taylor-Couette gap width, ranging from 0 at the inner boundary to 1 at
%   the outer boundary.
%       Hence x_=0 corresponds to r=eta_/(1-eta_) and x_=1 corresponds to
%       r=1/(1-eta_).

%Normalisation of perturbations:
x_s=0:0.0005:1;
for i_=length(x_s):-1:1
    u_s(i_)=ChebEval(u_k,1-2*x_s(i_));
end
Norm_=max(u_s); %We normalise by the maximum radial velocity perturbation.
u_k=u_k/Norm_; v_k=v_k/Norm_; w_k=w_k/Norm_;

%With the perturbations normalised, we can define functions that will yield
%the value of each perturbation at each combination of x_, theta_ and z_.
%Function definitions:
u_eval=@(x_,theta_,z_) real(exp(1i*(m_*theta_+k_*z_))*ChebEval(u_k,1-2*x_));
v_eval=@(x_,theta_,z_) real(exp(1i*(m_*theta_+k_*z_))*ChebEval(v_k,1-2*x_));
w_eval=@(x_,theta_,z_) real(exp(1i*(m_*theta_+k_*z_))*ChebEval(w_k,1-2*x_));

%We now need to develop grids for each cross-section. For the vertical
%cross-section we need grids in the [x,z]-plane, and for the horizontal
%cross-section we need grids in the [x,theta]-plane.
%   For conciseness, these will be denoted as xz and xt below.
%
%We will need a coarse grid (20 points in each dimension) for the quiver
%plot and a fine grid (100 points in each direction) for the contour plot.
%
%We will also need grid results for each velocity perturbation. This means
%a total of 12 arrays:
%   ["3 (perturbations)" times "2 (coarse and fine)" times "2 (xt and xz)"]

%Define minimum and maximum of each dimension:
r_min=eta_/(1-eta_);
r_max=   1/(1-eta_);
x_min=0; %r = eta_/(1-eta_)
x_max=1; %r =    1/(1-eta_)
if m_==0
    theta_min=-pi;
    theta_max=pi;
else
    theta_min=-pi/m_;
    theta_max=pi/m_;
    %This is scaled to show a single horizontal cell.
end
z_min=-pi/k_;
z_max=pi/k_;
%This is scaled to show a single vertical cell.

%Coarse grid dimensions:
x_step_C=(x_max-x_min)/20;
theta_step_C=(theta_max-theta_min)/20;
z_step_C=(z_max-z_min)/20;

x_vals_C=x_min:x_step_C:x_max;
theta_vals_C=theta_min:theta_step_C:theta_max;
z_vals_C=z_min:z_step_C:z_max;

x_LC=length(x_vals_C);
t_LC=length(theta_vals_C);
z_LC=length(z_vals_C);

%Fine grid dimensions:
x_step_F=(x_max-x_min)/100;
theta_step_F=(theta_max-theta_min)/100;
z_step_F=(z_max-z_min)/100;

x_vals_F=x_min:x_step_F:x_max;
theta_vals_F=theta_min:theta_step_F:theta_max;
z_vals_F=z_min:z_step_F:z_max;

x_LF=length(x_vals_F);
t_LF=length(theta_vals_F);
z_LF=length(z_vals_F);

%Coarse Grids:
u_gridC_xt=zeros(t_LC,x_LC);
v_gridC_xt=zeros(t_LC,x_LC);
w_gridC_xt=zeros(t_LC,x_LC);

u_gridC_xz=zeros(z_LC,x_LC);
v_gridC_xz=zeros(z_LC,x_LC);
w_gridC_xz=zeros(z_LC,x_LC);

%Fine Grids:
u_gridF_xt=zeros(t_LF,x_LF);
v_gridF_xt=zeros(t_LF,x_LF);
w_gridF_xt=zeros(t_LF,x_LF);

u_gridF_xz=zeros(z_LF,x_LF);
v_gridF_xz=zeros(z_LF,x_LF);
w_gridF_xz=zeros(z_LF,x_LF);

%Evaluate the horizontal grids:
%Coarse:
for i_=1:x_LC
    for j_=1:t_LC
        u_gridC_xt(j_,i_)=u_eval(x_vals_C(i_),theta_vals_C(j_),0);
        v_gridC_xt(j_,i_)=v_eval(x_vals_C(i_),theta_vals_C(j_),0);
        w_gridC_xt(j_,i_)=w_eval(x_vals_C(i_),theta_vals_C(j_),0);
    end
end
%Fine:
for i_=1:x_LF
    for j_=1:t_LF
        u_gridF_xt(j_,i_)=u_eval(x_vals_F(i_),theta_vals_F(j_),0);
        v_gridF_xt(j_,i_)=v_eval(x_vals_F(i_),theta_vals_F(j_),0);
        w_gridF_xt(j_,i_)=w_eval(x_vals_F(i_),theta_vals_F(j_),0);
    end
end

%Evaluate the vertical grids:
%Coarse:
for i_=1:x_LC
    for k_i=1:z_LC
        u_gridC_xz(k_i,i_)=u_eval(x_vals_C(i_),0,z_vals_C(k_i));
        v_gridC_xz(k_i,i_)=v_eval(x_vals_C(i_),0,z_vals_C(k_i));
        w_gridC_xz(k_i,i_)=w_eval(x_vals_C(i_),0,z_vals_C(k_i));
    end
end
%Fine:
for i_=1:x_LF
    for k_i=1:z_LF
        u_gridF_xz(k_i,i_)=u_eval(x_vals_F(i_),0,z_vals_F(k_i));
        v_gridF_xz(k_i,i_)=v_eval(x_vals_F(i_),0,z_vals_F(k_i));
        w_gridF_xz(k_i,i_)=w_eval(x_vals_F(i_),0,z_vals_F(k_i));
    end
end

r_vals_F=(eta_/(1-eta_))+x_vals_F;
r_vals_C=(eta_/(1-eta_))+x_vals_C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function Grids Calculated

%Calibrate min/max scale of colour contours:
max_=max(max([u_gridF_xz, u_gridF_xt, v_gridF_xz, v_gridF_xt, w_gridF_xz, w_gridF_xt]));
min_c=-max_; max_c=max_; %This needs to be THE maximum value across all datasets.

%Count of shades of colour in contour plots:
conts_=500;

%Load colour scheme:
load('parula_alt.mat')
%(This is a customised version of the parula colour scheme, intended to
%better display near-zero changes. The parula colour scheme itself should
%ensure good display even when printed in black and white.)

%Vertical Plot:
fig_vertical=figure();
axes_1 = axes('Parent',fig_vertical);
set(fig_vertical, 'Position', [400 500 280 280]);
set(axes_1, 'Position', [0.18 0.18 0.77 0.77]);
hold on;
box on;

[~,h] = contourf(r_vals_F,z_vals_F,v_gridF_xz,conts_);
set(h,'LineColor','none');
colormap(map);
quiver(r_vals_C,z_vals_C,u_gridC_xz,w_gridC_xz,'black');
x_label = xlabel('Radius'); xlim([r_min r_max]);
y_label = ylabel('Height'); ylim([z_min z_max]);
set(y_label,'position',[0,0,1]);
caxis([min_c max_c])
vh_='v';

%Save Figure:
saveas(fig_vertical,['eta_' strrep(num2str(eta_),'.','') '_mu_' strrep(num2str(mu_),'.','') '_N_' strrep(num2str(N_),'.','') '_' vh_ '.fig']);

%Horizontal:
fig_horizontal=figure();
axes_2 = axes('Parent',fig_horizontal);
set(fig_horizontal, 'Position', [800 500 372 280]);
set(axes_2, 'Position', [0.18 0.18 0.77 0.77]);
hold on;
box on;

[~,h] = contourf(r_vals_F,theta_vals_F,w_gridF_xt,conts_);
set(h,'LineColor','none');
colormap(map);
quiver(r_vals_C,theta_vals_C,u_gridC_xt,v_gridC_xt,'black');
xlabel('Radius'); xlim([r_min r_max]);
ylabel('\theta','Rotation',0); ylim([theta_min theta_max]);
caxis([min_c max_c])
colorbar;
vh_='h';

%Save Figure:
saveas(fig_horizontal,['eta_' strrep(num2str(eta_),'.','') '_mu_' strrep(num2str(mu_),'.','') '_N_' strrep(num2str(N_),'.','') '_' vh_ '.fig']);

end


function [ f_val ] = ChebEval( Coeff_vals, xk_ )
%CHEBEVAL Takes a list of coefficients and a location within the Chebyshev
%range [-1 to +1]. Return the value of the Chebyshev expansion of the
%function at that point.

L_ = length(Coeff_vals);

Cheb_vals=zeros(1,L_);

Cheb_vals(1)=1;
Cheb_vals(2)=xk_;

for i_=3:L_
    Cheb_vals(i_)=2*xk_*Cheb_vals(i_-1)-Cheb_vals(i_-2);
end

f_val=0;
for i_=1:L_
    f_val=f_val+Coeff_vals(i_)*Cheb_vals(i_);
end

end