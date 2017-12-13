%% Capacitive deionization 2D simulation
%
% Developed by: Xia Shang
% Advisors: Kyle Smith, Roland Cusick
% Copyright: Xia Shang, Roland Cusick, Kyle Smith
% University of Illinois at Urbana-Champaign
% All rights reserved.
%
% Funded by: US National Science Foundation Award No. 1605290 entitiled
% "SusChEM: Increasing Access to Sustainable Freshwater Resources with
% Membrane Capacitve Deionization", and Joint Center for Energy Storage
% Research, an Energy Innovation Hub funded by the U.S. Department of
% Energy, Office of Science, Basic Energy Sciences.
%
% Redistribution and use in source and binary forms, with or
% without modification, are permitted provided that the following
% conditions are met:
%
%     *   Redistributions of source code must retain the above copyright notice,
%         this list of conditions and the following disclaimer.
%     *   Redistributions in binary form must reproduce the above
%         copyright notice, this list of conditions and the following
%         disclaimer in the documentation and/or other materials provided
%         with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
% OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [T_r, V_cell, C_eff] = CDI_2D_Demo(C0, Q_FC,...
            J, Q_fix_pos, Q_fix_neg,...
            V_max, V_min, Cycle_limit, L_FC, T_pos, T_neg, ...
            T_FC, R_FC, R_Macro, R_Micro, dt_record)
%clc; clear; close all;
tic


%% cell operating parameters
%C0=30;                         % mM, influent concentration
%V_max=200;                    % mV, cell voltage
%V_min = 0;
%J= 1;                          % mA/cm2, projected area
%Cycle_limit = 2;               % -, total cycle number
%Q_FC = 0.0736;                 % mL/min, flow rate suggested range (0.05 - 1 mL/min)

%% Cell design
T_CC1 = 0;                % um, thickness of current collector 1
T_pos = 450;              % um, thickness of Cathod
T_FC = 250;               % um, thickness of Flow channel
T_neg = 450;              % um, thickness of Anode
T_CC2 = 0;                % um, thickness of current collector 2
L_FC = 9;                 % cm, Length of flow channel
W_FC = 2;                 % cm, Width of flow channel
V_internal = 0.455;       % mL, volume of a well mixed internal tank, dead-volume

%% Electrode Porosities
M_e = 0.4214;                    % g/cm3; electrode desity
C_Measured = 49;                 % F/g; Stern-layer capacitance, this value can be initially estimated based on the whole cell capacitance measured in CV at low scan rate (e.g., < 1mV/s) in 1 M NaCl
%Q_fix_pos = -4;                  % C/cm3, immobile charge density
%Q_fix_neg = 0;                   % C/cm3, immobile charge density
uatt=0;                          % -, attracitve chemical force using in the mD theory
R_Macro = 0.4;                   % Macropore ratio
R_Micro = 0.2;                   % Micropore ratio
R_FC = 0.7;                      % Average porosity of flow channel
R_e= 1-(1-R_Macro)*(1-R_Micro);  % Total porosity of electrode

%% numerical parameters
dt = 0.008;               % s, dt
dx = 25;                  % um
dy=200*L_FC;              % um

%% data storage to reduce membrane usage
tmin = 0;                        % s, starting time
tmax =round(800*Cycle_limit*1.2/J);   % s, ending time, this value can be adjusted
dt_record=1;            % s, time step to record data
Nt=round((tmax-tmin)/dt);   % total number of nodes in time
Nt_r = round((tmax-tmin)/dt_record);      % no. of stored nodes in time

% other related parameters
E_Dielectric = 80.4;                % dielectric constant of water
E_Permittivity = 8.85*10^-12;       % F/m; permittivity of free space
r_Na = 0.16;                        % nm, radius of Na+ ion
d_H2O = 0.29;                       % nm, diameter of a H2O molecule
d_Helmholtz = (r_Na + d_H2O)/10^9;  % m
c_Helmholtz = E_Dielectric*E_Permittivity/d_Helmholtz;  % F/m2
A_effective = C_Measured/c_Helmholtz;                   % m2/g; effective surface area
Temp = 298.16;                      % K, temperature
KB = 1.38*10^-23;                   % C V K-1; Boltzmann constant
e = 1.602*10^-19;                   % C;
Na = 6.02*10^23;                    % mol-1;
F = 96485.3365;                     % C mol-1; Faraday constant
M_NaCl=58.44;                       % g/mole
z_a =-1;                            % valence of anion
z_c = 1;                            % valence of cation
q_a = 1;                            % charge of anion
q_c = 1;                            % charge of cation
t_Na = 0.5;                         % transference number of Na+
t_Cl = 1-t_Na;                      % transference number of Cl-

% Diffusion coefficient
D_NaCl = 1.61*10^-5;                % cm2/s; diffusion coeff. in spacer
D_e=D_NaCl*10^8*R_Macro^0.5;        % um2/s, effective diffusion coefficient in electrode
D_FC = D_NaCl*10^8*R_FC^0.5;        % um2/s, effective diffusion coefficient
D_mid = 2*D_FC*D_e/(D_FC+D_e);

% resistance
R_contact = 9;                  % ohm cm2
R_electrode = 0.1625;               % ohm*m

% flow rate and mixing in the dead-volume
v_FC = Q_FC/T_FC*10^4/W_FC/60*10^4; %  um/s, flow velocity, assume uniform flow profile in the cross section
k_mix = Q_FC/60/V_internal;         % /s, factor used in the equation accounting for the effect of mixing in the dead-volume

%% Leakage current sub-models
% B-V model parameter
a_ele = 11000;                      % cm2/cm3 electrode volume
I0_C = 2.5*10^-7;                   % mA/cm2
I0_O2 = 5.5 *10^-7;                 % mA/cm2
E0_C = 0.207;                       % V, vs SHE
E0_O2 = 0.81;                       % V, vs SHE
RT_F = 0.0592;                      % V
a_C = 0.5;
a_O2 = 0.5;
E0_An = 0.5419;                     % V, vs SHE, measured at equlibrium condition when two electrodes were short circuited
E0_Ca = 0.5419;                     % V, vs SHE, measured at equlibrium condition when two electrodes were short circuited

% imperical limiting current model parameter
E_half=0.124;                       % V, vs SHE
IL_O2=3.558;                        % mA/cm3-electrode

%% Indexing
% Current collector
xmin_CC1 = 0;                       % um, minimum value of x
xmax_CC1 = (xmin_CC1+T_CC1)/dx;     % um, maximum value of x

%Postive electrode
xmin_pos = xmax_CC1;                % left boundary node point of the positive electrode
xmax_pos = (xmin_pos+T_pos)/dx;     % right boundary node point of the positive electrode
L_pos = L_FC;                       % cm, length of the positive electrode
W_pos = W_FC;                       % cm, width of the positive electrode
A_edge = 0;                         % cm2, area of cathod in the triangular spot
A_pos = L_pos*W_pos+A_edge;         % cm2, area of projected area of the positive electrode
V_pos = T_pos/10000 * A_pos;        % cm^3, volume of the positive electrode

%Flow channel
xmin_FC = xmax_pos;                 % left boundary node point of FC
xmax_FC = xmin_FC+T_FC/dx;          % right boundary node point of FC
A_edge_FC = 0;                      % cm2, area of flow channel in the triangular spot
A_FC = L_FC*W_FC +A_edge_FC;        % cm2, projected area of flow channel
Vol_FC = T_FC/10^4*A_FC/1000*R_FC;  % L, Volume of flow channel

%Negative electrode
xmin_neg = xmax_FC;                 % left boundary node point of the negtive electrode
xmax_neg = xmin_neg+T_neg/dx;       % right boundary node point of the negtive electrode

%Current collector 2
xmin_CC2 = xmax_neg;                % um, minimum value of x
xmax_CC2 = xmin_CC2+T_CC2/dx;       % um, maximum value of x

% indexing of vectors and matrix
nx = xmax_CC2-xmin_CC1;             % total number of compartments in the cell
ny=L_FC*10000/dy;                   % ny
N=T_pos/dx;                         % total number of compartments in the electrode
N_FC=T_FC/dx;                       % total number of compartments in the FC
V_i = W_pos*dy/10000*dx/10000;      % cm3, volume per compartment

x_ind=[1:xmax_pos, xmin_neg+1:xmax_neg];    % indexing for the electrodes
FC_ind=[xmax_pos+1:xmin_neg];               % indexing for the flow channel

FC_ind_A = zeros(numel(FC_ind)*ny,1);       % indexing for the flow channel in 2D
for i=1:ny
    FC_ind_A(1+(i-1)*numel(FC_ind):i*numel(FC_ind)) = FC_ind + (i-1)*nx;
end

%% Initiation
Vt=zeros(Nt_r, 4);                  % mV, 1st column: electrode voltage, 2nd column: IEM resistance voltage; 3rd column: spacer resistance voltage; 4th total voltage
T0=0;                               % s previous time step
T1=0;                               % s current time step
T_r=zeros(Nt_r,1);                  % s, reduced time domain for data storage
Cycle_index=zeros(100,1);           % index of cycles
Discharge_index=zeros(100,1);       % index of discharge cycles
VV_0 = zeros(nx,ny);                % mV, Voltage matrix
V_N = zeros (nx, ny);               % mV, Voltage matrix of ionic resistance in the electrode
II_0 = zeros(nx,ny);                % mA, Current matrix: current passing each capacitor
II_1 = zeros(nx,ny);                % mA, Current matrix: current passing each capacitor
Ri_M =zeros(nx,ny);                 % ohm, individual ionic resistor matrix
Re_M =zeros(nx,ny);                 % ohm, individual electronic resistor matrix
Cond_M = zeros(nx,ny);              % Conductivity matrix
IL = zeros(nx,ny);                  % mA, Current matrix: current passing each resistor in parallel to the capacitor
I_Ri = zeros (nx,ny);               % mA, Current matrix : current passing each ionic resistor
I_Re = zeros (nx,ny);               % mA, Current matrix: current passing each electronic resistor
C_0 = zeros (nx,ny);                % mM, Concentration matrix coresponding to the capacitor matrix
C_1 = zeros (nx,ny);                % mM, Concentration matrix coresponding to the capacitor matrix
C1 = zeros(nx,ny);                  % mM, concentration matrix
C_re =zeros(nx,ny,Nt_r);            % mM, recorded concentration tensor
II_re =zeros(nx,ny,Nt_r);           % mA, recorded current tensor
Q_mi_re = zeros(nx,ny,Nt_r);        % mV, recorded charge density tensor
IL_re = zeros(nx,ny,Nt_r);          % mA, recorded leakage density tensor
Qt_0 = zeros (nx,ny);               % C, charge density matrix
Qt_1 = zeros (nx,ny);               % C, charge density matrix
Qt_e=Qt_0;                          % C, charge density matrix
Qt_mi=-Qt_0;                        % C, charge density matrix
V_D=zeros(nx,ny);                   % mV, Donnan potential matrix
C_mi_0=zeros(nx,ny);                % mmol-L micropores, ionic charges matrix inside micropores
C_mi_1=zeros(nx,ny);                % mmol-L micropores, ionic charges matrix inside micropores
C_eff = zeros(Nt_r,1);              % mM avg effulent concentration per batch
C_eff_mix = C_eff;                  % mM, true effluent considering the effect of well mixed condition at the end of effluent

% Initial concentration and polarization
C_0(:,:)=C0;                        % mM, Initial concentration profile in the electrode phase
C_1(:,:)=C0;                        % mM, Initial concentration profile in the electrode phase
VV_0(:,:)=0;                        % mV, Initial polarization profile in the electrode phase

% determine the inital current distribution across the cell
I = J*A_FC;                         % mA, total current applied
indx=(1:nx)';
indx_1=circshift(indx,1);
indx_0=circshift(indx,-1);
indy=(1:ny)';
indy_1=circshift(indy,1);
indy_0=circshift(indy,-1);

% calculate ionic conductivity S/cm
Cond_M(:,:)=2*D_NaCl*C_0(:,:)/KB/Temp/1000000*F*e*R_e^1.5;
Cond_M(xmax_pos+1:xmin_neg,:)=2*D_NaCl*C_0(xmax_pos+1:xmin_neg,:)/KB/Temp/1000000*F*e*R_FC^1.5;
Cond_M(1:xmax_pos-1,:)=2.*Cond_M(1:xmax_pos-1,:).*Cond_M(indx_0(1:xmax_pos-1),:)./(Cond_M(1:xmax_pos-1,:)+Cond_M(indx_0(1:xmax_pos-1),:));
Cond_M(xmin_neg+1:xmax_neg-1,:)=2.*Cond_M(xmin_neg+1:xmax_neg-1,:).*Cond_M(indx_0(xmin_neg+1:xmax_neg-1),:)./(Cond_M(xmin_neg+1:xmax_neg-1,:)+Cond_M(indx_0(xmin_neg+1:xmax_neg-1),:));
Ri_M(:,:) = 1./Cond_M(:,:)*dx/dy/W_FC;
Ri_S=zeros(ny,1);
Ri_S(:)=sum(Ri_M(xmax_pos:xmin_neg,:));

% calculate the current distribution in the electrode
M_I = zeros(ny,ny);              % Coefficient matrix for current distribution in the electrode
M_I(1:(ny+1):(end-ny))=Ri_S(1:(end-1));
M_I((ny+1):(ny+1):end)=-Ri_S(2:end);
M_I(ny,:)=1;
S_I=sparse(M_I);
B_I =zeros(ny,1);               % vector for current distribution in the electrode
B_I(1:end-1)=VV_0(xmax_pos,indy_0(1:(end-1)))-VV_0(xmin_neg+1,indy_0(1:(end-1)))-VV_0(xmax_pos,1:(end-1))+VV_0(xmin_neg+1,1:(end-1));
B_I(end)=I;
I_y=S_I\B_I;                    % solve the current distribution in the electrode in the y-direction

% calculate the leakage current in the positive electrode
Oe = 0;                         % set the electronic resistance
%IL(1:xmax_pos,1)=a_ele*I0_C*V_i*exp(a_C/RT_F.*(VV(1:xmax_pos,1)/1000+E0_An-E0_C));
IL(1:xmax_pos,:)=0;
% Update the current distribution in the positive electrode
II_0(2:xmax_pos-1, :)=-IL(2:xmax_pos-1,:)+(VV_0(indx_0(2:xmax_pos-1),:)-VV_0(indx(2:xmax_pos-1),:))./(Oe+Ri_M(2:xmax_pos-1,:))-(VV_0(indx(2:xmax_pos-1),:)-VV_0(indx_1(2:xmax_pos-1),:))./(Oe+Ri_M(indx_1(2:xmax_pos-1),:));
II_0(1,:)=-IL(1,:)+(VV_0(2,:)-VV_0(1,:))/(Oe+Ri_M(1,:));
II_0(xmax_pos,:)=I_y-IL(xmax_pos,1)-sum(II_0(1:xmax_pos-1));

% calculate the leakage current in the negative electrode
IL(xmin_neg+1:xmax_neg,:)=0;
%IL(xmin_neg+1:xmax_neg,1)=a_ele*I0_O2*V_i*(exp(a_O2/RT_F.*(VV(xmin_neg+1:xmax_neg,1)/1000+E0_Ca-E0_O2))-exp(-a_O2/RT_F.*(VV(xmin_neg+1:xmax_neg,1)/1000+E0_Ca-E0_O2))) ;
%IL(xmin_neg+1:xmax_neg,1)=-IL_O2*V_i./(1+exp((VV_0(xmin_neg+1:xmax_neg,:)/1000+E0_Ca-E_half)/RT_F));

% update the current distribution in the negative electrode
II_0(xmin_neg+2:end-1,:)=-IL(xmin_neg+2:end-1,:)+(VV_0(indx_1(xmin_neg+2:end-1),:)-VV_0(indx(xmin_neg+2:end-1),:))./(Oe+Ri_M(indx_1(xmin_neg+2:end-1),:))-(VV_0(indx(xmin_neg+2:end-1),:)-VV_0(indx_0(xmin_neg+2:end-1),:))./(Oe+Ri_M(indx(xmin_neg+2:end-1),:));
II_0(end,:)=-IL(end,:)+(VV_0(end-1,:)-VV_0(end,:))./(Oe+Ri_M(end,:));
II_0(xmin_neg+1,:)=-I_y'-sum(II_0(xmin_neg+2:end,:))-sum(IL(xmin_neg+1:end,:));

% calculate the ionic current in the electrode
A = tril (ones(N));
A1=A';
A2=tril(ones(N_FC))';
I_Ri(1:xmax_pos,1)=A*(IL(1:xmax_pos,1)+II_0(1:xmax_pos,1));
I_Ri(xmax_pos+1:xmin_neg,1)=I_y(1);
I_Ri(xmin_neg+1:xmax_neg,1)=-A1*(IL(xmin_neg+1:xmax_neg,1)+II_0(xmin_neg+1:xmax_neg,1));


%B.C vector
bc=zeros(nx,ny);
bc(1,:)=0; bc(nx,:)=0;  %Neumann B.Cs
bc(:,1)=0; bc(:,ny)=0;
%bc(xmax_pos+1:xmin_neg,1)=C0/dy^2;  %Dirichlet B.Cs

%B.Cs at the corners:
bc(1,1)=0; bc(nx,1)=0;
bc(1,ny)=0; bc(nx,ny)=0;
bc(x_ind,:)=D_e*dt*bc(x_ind,:);
bc(FC_ind,:)=D_e*dt*bc(FC_ind,:);

%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx,1:nx-1,1,nx,nx);
Ax=Ex+Ex'-2*speye(nx);          %Dirichlet B.Cs
Ax(1,1)=-1; Ax(nx,nx)=-1;       %Neumann B.Cs
Ey=sparse(2:ny,1:ny-1,1,ny,ny);
Ay=Ey+Ey'-2*speye(ny);          %Dirichlet B.Cs
Ay(1,1)=-1; Ay(ny,ny)=-1;       %Neumann B.Cs
D=R_Macro*speye((nx)*(ny));
D(FC_ind_A,:)=D(FC_ind_A,:)/R_Macro*R_FC;
D_y = speye(nx)*D_e*dt;
D_y = kron(Ay/dy^2,D_y);
D_x = speye(ny)*D_e*dt;
D_x = kron(D_x,Ax/dx^2);
D_x (FC_ind_A,:) = D_x(FC_ind_A,:)/D_e*D_FC;    % update the diffusion coefficient in the flow channel
D_x (length(D_x)*(xmax_pos-1)+xmax_pos+1:length(D_x)*nx+nx:end) = D_x(2)/D_e*D_mid;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmax_pos)+xmax_pos:length(D_x)*nx+nx:end) = D_x(2)/D_e*D_mid;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmin_neg)+xmin_neg:length(D_x)*nx+nx:end) = D_x(2)/D_e*D_mid;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmin_neg-1)+xmin_neg+1:length(D_x)*nx+nx:end) = D_x(2)/D_e*D_mid;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmax_pos-1)+xmax_pos:length(D_x)*nx+nx:end) = -(D_mid+D_e)*dt/dx^2;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmax_pos)+xmax_pos+1:length(D_x)*nx+nx:end) = -(D_mid+D_FC)*dt/dx^2;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmin_neg-1)+xmin_neg:length(D_x)*nx+nx:end) = -(D_mid+D_FC)*dt/dx^2;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D_x (length(D_x)*(xmin_neg)+xmin_neg+1:length(D_x)*nx+nx:end) = -(D_mid+D_e)*dt/dx^2;    % update the diffusion coefficient at the interfaces between electrode and flow channel
D=D-D_x - D_y;
%D=D-D_x;
%D=D-D_e*dt*(kron(Ay/dy^2,speye(nx))+kron(speye(ny),Ax/dx^2));
%D((size(D,1)+1)*(xmax_pos)+1:(size(D,1)+1):1+(size(D,1)+1)*(xmin_neg-1))=R_FC-;
clear D_x D_y

% Other initial values
C_eff(1)=C0;
Qt_0(:,:)=0;                        %C/cm3-electrode; Initial charge density at time step n=1
Count_C = zeros(Nt,1);              % Iternation number for the concentration matrix for the last current matrix iteration at each time step
Count_Q = zeros(Nt,1);              % Iteration number for the current matrix at each time step
Error_Q = 0;                        % error of the charge iterative loop
flag_0=1;                           % flag for the current convergence loop
flag_1=1;                           % flag for the concentration convergence loop
flag_3=0;                           % flag for adjusting the concentration convergence
flag_4=1;                           % flag for the technique after charging
flag_5=1;                           % flag for the technique after discharging
flag_6=1;                           % flag for ocv study
flag_7=0;                           % flag for iterative time step
flag_8=0;                           % flag for iteration loops
Cycle_number=1;
n_r=1;
ind_B=(1:size(D,1))';
ind_B1=circshift(ind_B,1);
ind_B0=circshift(ind_B,-1);
err=zeros(500,1);
Vt_1=zeros(4,1);

% Calculate the initial charge storage and Donnan layer potential
Qt_mi(1:xmax_pos,:) = -Qt_e(1:xmax_pos,:)-Q_fix_pos;
Qt_mi(xmin_neg+1:end,:) = -Qt_e(xmin_neg+1:end,:) - Q_fix_neg;
V_D(x_ind,:)=-asinh(Qt_mi(x_ind,:)/R_Micro*10^6/F/2./(C0)/exp(uatt));   % non-dimentional Donnan layer potential
VV_0(x_ind,:)=V_D(x_ind,:)*25.7+Qt_e(x_ind,:)/C_Measured/M_e*1000;      % mV, electrode polarization
C_mi_1(x_ind,:)=2*(C0)*exp(uatt).*cosh(V_D(x_ind,:));                   % mmol-L micropores
C_eff_0 = C0;
C_mix_0 = C0;

% Calculate the inital voltage distribution
V_N (xmin_neg+1:xmax_neg,1)=  A1*(Ri_M(xmin_neg+1:xmax_neg,1).*I_Ri(xmin_neg+1:xmax_neg,1))+Oe*I_y(1)-VV_0(xmax_neg,1);              % liquid potential across one electrode when the cathod is grounded
V_N(xmax_pos+1:xmin_neg,1)=A2*Ri_M(xmax_pos+1:xmin_neg,1).*I_y(1)+V_N(xmin_neg+1,1);                    % liquid across Flow channel
V_N (1:xmax_pos,1)=A1*(Ri_M(1:xmax_pos,1).*I_Ri(1:xmax_pos,1))+V_N(xmax_pos+1,1) ;              % liquid across one electrode
Vt (1,1)=V_N(1,1)+I_y(1)*Oe+VV_0(1,1);              % mV, voltage in the solution
Vt (1,2)=R_contact*J;                               % mV, voltage due to the resistance of two IEMs and external resistance
Vt (1,3)=0;                                         % mV, membrane Donnan potential
Vt (1,4)=sum(Vt(1,1:3));                            % mV, Cell voltage

fprintf('Thank you for using the CDI simulation tool\n')
fprintf('***************************\n')
fprintf('*******CDI 2D Model********\n')
fprintf('Developed by Xia Shang at UIUC\n')
fprintf('Copyright: Xia Shang, Prof. Roland Cusick, Prof. Kyle Smith\n\n')
fprintf('*******Selected input******\n')
fprintf('Current Density = %4.2f mA/cm2\nVoltage Window= %4.2f mV\nInfluent Concentration = %4.2f mM.\n', J, V_max-V_min, C0)
time_stamp = datestr(now, 'HH:MM:SS');
date = datetime('today');
fprintf('Simulation starting time: %s', time_stamp)
fprintf('Date: %s', date)
fprintf('\n')
fprintf('*********Output***********\n')
fprintf('Time; Cell Voltage; Effluent Concentration\n')
%% Main time loop
for n=2:1:Nt            % time step n
    flag_7=0;           % reset flag_7
    C_0=C_1;            % update concentration field
    C_mi_0=C_mi_1;      % update micro-pores ionic concentration field
    Qt_0=Qt_1;          % update charge density field
    T0=T1;              % update time
    while flag_7==0;    % adaptive time loop
        flag_8=0;       % reset flag 8
        flag_0=1;       % reset flag_0 for the current while loop
        flag_3=0;       % reset flag_3
        T1=T0+dt;

        % guess electric charge density based on the information at the previous time step
        Qt_e(:,:)=Qt_0(:,:)+II_0(:,:)*dt/10^3/V_i; % C/cm3-electrode
        Count_Q(n)=0;
        while flag_0==1     % charge density iterative loop
            flag_1=1;       % reset flag_1 for the concentration while loop
            % guess Ce* based on the calculated charge efficiency at time step n-1
            C1=C_0;
            Count_C(n) = 0;
            while flag_1==1     % concentration iterative loop
                % Solve for total ionic charge density in the micropores
                Qt_mi(1:xmax_pos,:) = -Qt_e(1:xmax_pos,:)-Q_fix_pos;
                Qt_mi(xmin_neg+1:end,:) = -Qt_e(xmin_neg+1:end,:) - Q_fix_neg;
                V_D(x_ind,:)=-asinh(Qt_mi(x_ind,:)/R_Micro*10^6/F/2./(C1(x_ind,:))/exp(uatt));  % non-dimentional Donnan layer potential
                VV_0(x_ind,:)=V_D(x_ind,:)*25.7+Qt_e(x_ind,:)/C_Measured/M_e*1000;   % mV, electrode polarization
                C_mi_1(x_ind,:)=2*(C1(x_ind,:))*exp(uatt).*cosh(V_D(x_ind,:));  % mmol-L micropores

                % Solve for new Ce including the effect of diffusion (implicit, forward time, CN in diffusion, upwind explicit in aadvection)
                C_1=C1;
                C1(x_ind,:)=R_Macro*C_0(x_ind,:)-(C_mi_1(x_ind,:)-C_mi_0(x_ind,:))/2*R_Micro;
                C1(FC_ind,:)=R_FC*C_0(FC_ind,:);
                C1(FC_ind,2:end)=C1(FC_ind,2:end)+v_FC*dt/dy*(C_0(xmax_pos+1:xmin_neg,indy_1(2:end))-C_0(xmax_pos+1:xmin_neg,indy(2:end)));
                C1(FC_ind,1)=C1(FC_ind,1)+v_FC*dt/dy*(C0-C_0(xmax_pos+1:xmin_neg,1));
                C1=reshape(C1+bc,[],1);
                C1=D\C1;
                C1=reshape(C1,nx,ny);

                % check for convergence
                err(Count_C(n)+1,1)=max(max(abs(C1-C_1)./C_1));
                    if err(Count_C(n)+1,1)<10^-7;
                        flag_1=0;
                        flag_8=1;
                        C_1=C1;
                    end
                Count_C(n) = Count_C(n)+1;
            end

            % update the current and voltage distribution
            Cond_M(:,:)=2*D_NaCl*C_1(:,:)/KB/Temp/1000000*F*e*R_e^1.5;            % S/cm conductivity matrix
            Cond_M(xmax_pos+1:xmin_neg,:)=2*D_NaCl*C_1(xmax_pos+1:xmin_neg,:)/KB/Temp/1000000*F*e*R_FC^1.5;            % S/cm conductivity matrix
            Cond_M(1:xmax_pos-1,:)=2.*Cond_M(1:xmax_pos-1,:).*Cond_M(indx_0(1:xmax_pos-1),:)./(Cond_M(1:xmax_pos-1,:)+Cond_M(indx_0(1:xmax_pos-1),:));
            Cond_M(xmin_neg+1:xmax_neg-1,:)=2.*Cond_M(xmin_neg+1:xmax_neg-1,:).*Cond_M(indx_0(xmin_neg+1:xmax_neg-1),:)./(Cond_M(xmin_neg+1:xmax_neg-1,:)+Cond_M(indx_0(xmin_neg+1:xmax_neg-1),:));
            Ri_M(:,:) = 1./Cond_M(:,:)*dx/dy/W_FC; % ohm, resistance matrix
            Ri_S(:)=sum(Ri_M(xmax_pos:xmin_neg,:));

            % calculate the current distribution in the electrode
            M_I(1:(ny+1):(end-ny))=Ri_S(1:(end-1));
            M_I((ny+1):(ny+1):end)=-Ri_S(2:end);
            S_I=sparse(M_I);
            B_I(1:end-1)=VV_0(xmax_pos,indy_0(1:(end-1)))-VV_0(xmin_neg+1,indy_0(1:(end-1)))-VV_0(xmax_pos,1:(end-1))+VV_0(xmin_neg+1,1:(end-1));
            B_I(end)=I;
            I_y=S_I\B_I;                    % solve the current distribution in the electrode in the y-direction

            % calculate the leakage current
            IL(1:xmax_pos,:)=a_ele*I0_C*V_i*exp(a_C/RT_F.*(VV_0(1:xmax_pos,:)/1000+E0_An-E0_C)); % VB model
            %IL(1:xmax_pos,:)=0;    % no leakage current
            % Update the current distribution in the anode
            II_1(2:xmax_pos-1, :)=-IL(2:xmax_pos-1,:)+(VV_0(indx_0(2:xmax_pos-1),:)-VV_0(indx(2:xmax_pos-1),:))./(Oe+Ri_M(2:xmax_pos-1,:))-(VV_0(indx(2:xmax_pos-1),:)-VV_0(indx_1(2:xmax_pos-1),:))./(Oe+Ri_M(indx_1(2:xmax_pos-1),:));
            II_1(1,:)=-IL(1,:)+(VV_0(2,:)-VV_0(1,:))/(Oe+Ri_M(1,:));
            II_1(xmax_pos,:)=I_y'-sum(II_1(1:xmax_pos-1,:))-sum(IL(1:xmax_pos,:));

            % calculate the current distribution in the Cathode
            %IL(xmin_neg+1:xmax_neg,:)=0;   % no leakage curret
            IL(xmin_neg+1:xmax_neg,:)=-IL_O2*V_i./(1+exp((VV_0(xmin_neg+1:xmax_neg,:)/1000+E0_Ca-E_half)/RT_F));    % limiting current model
            %IL(xmin_neg+1:xmax_neg,:)=a_ele*I0_O2*V_i*(exp(a_O2/RT_F.*(VV_0(xmin_neg+1:xmax_neg,:)/1000+E0_Ca-E0_O2))-exp(-a_O2/RT_F.*(VV_0(xmin_neg+1:xmax_neg,:)/1000+E0_Ca-E0_O2)));
            %BV model

            % update the current distribution in the cathode
            II_1(xmin_neg+2:end-1,:)=-IL(xmin_neg+2:end-1,:)+(VV_0(indx_1(xmin_neg+2:end-1),:)-VV_0(indx(xmin_neg+2:end-1),:))./(Oe+Ri_M(indx_1(xmin_neg+2:end-1),:))-(VV_0(indx(xmin_neg+2:end-1),:)-VV_0(indx_0(xmin_neg+2:end-1),:))./(Oe+Ri_M(indx(xmin_neg+2:end-1),:));
            II_1(end,:)=-IL(end,:)+(VV_0(end-1,:)-VV_0(end,:))./(Oe+Ri_M(end,:));
            II_1(xmin_neg+1,:)=-I_y'-sum(II_1(xmin_neg+2:end,:))-sum(IL(xmin_neg+1:end,:));

            % Calculate the charge density
            Qt_1(:,:)=Qt_0(:,:)+(II_0(:,:)+II_1(:,:))/2/1000/V_i*dt;

            % Check for convergence
            Error_Q =max(max(abs(Qt_1-Qt_e)));
            if (Error_Q<10^-7||Count_Q(n)>200)&&flag_8==1;
                flag_0=0;
                flag_7=1;
                II_0=II_1;
            else
                Qt_e=Qt_1;
            end
            Count_Q(n) = Count_Q(n)+1;
        end
    end
    C_eff_1=mean(C_1(xmax_pos+1:xmin_neg,ny));
    C_mid = (C_eff_1+C_eff_0)/2;
    C_mix_1 = C_mid - (C_mid - C_mix_0)*exp(-k_mix*dt);
    C_eff_0 = C_eff_1;
    C_mix_0 = C_mix_1;

    % Calculate the  voltage distribution
    Vt_1 (1)=-VV_0(xmin_neg+1,end)+I_y(end)*Ri_S(end)+VV_0(xmax_pos,end);       % mV, voltage in the solution
    Vt_1 (2)=R_contact/A_FC*I;            % mV, voltage due to the resistance of two IEMs and external resistance
    Vt_1 (3)=0;                           % mV, Donnan potential
    Vt_1 (4)=sum(Vt_1(1:3));

    % check if cell voltage reaches the maximum voltage limit
    if Vt_1(4)>V_max&&flag_4==1;
        I=-I;
        Discharge_index(Cycle_number)=n_r;
        flag_4=0;
        flag_5=1;
    end

    % check if cell voltage reaches the minimum voltage limit
    if Vt_1(4)<0&&flag_5==1 && n_r>100
        I=-I;
        Cycle_number=Cycle_number+1
        Cycle_index(Cycle_number)=n_r;
        flag_5=0;
        flag_4=1;
    end

    % check for data storage
    if T1-T_r(n_r)>=dt_record
        n_r=n_r+1;
        T_r(n_r)=T1;
        C_eff(n_r)=C_eff_1;
        C_eff_mix(n_r)=C_mix_1;
        Vt(n_r,:)=Vt_1(:);
        C_re(:,:,n_r)=C_1(:,:);
        II_re(:,:,n_r)=II_1(:,:);
        Q_mi_re(:,:,n_r)=Qt_mi;
        IL_re(:,:,n_r) = IL;
        fprintf('T = %4.2f s, V_cell = %4.2f mV, C_eff = %4.2f mM.\n', T1, Vt_1(4), C_eff_1)
    end
    % check if the maximum cycle litmits is reached
    if Cycle_number==Cycle_limit+1;
        %I=0;
        flag_6=0;
        break
    end
end

T_r = T_r(1:n_r)';
V_cell = Vt(1:n_r, 4)';
C_eff = C_eff(1:n_r)';

%%%%%%post-treatment, save data in file
%filename = sprintf('CDI_2D_target_effluent_%icm_%imV_%imM_amph_D_%iA_%iFg_%iC_%s.mat',L_FC, V_max,C0,int32(J*10),C_Measured,Q_fix_pos,date);
%save(filename,'Vol_FC','V_max','n_r','T_r','II_re','C_re','C_eff','C_eff_mix','Vt','Q_mi_re','Cycle_index','Discharge_index','L_FC','Q_FC','J','A_FC','dx','dy','Batch_index','B_number','IL_re','C_Measured','Q_fix_pos','Q_fix_neg','k_mix')

toc
end

