%Programmed by: Ido Schwartz
%About: It calculates the dispersion of Valance and Conduction Energy Bands
%Gallium Arsenide semiconductor bands calculation in K*p model.
close all;
clear all;
clc;

%%%%%%%%%% data %%%%%%%%%%
E_g = 1.424;                            % gap energy [eV]
E_p = 28.8;                              % valence dipole energy [eV]
g_1 = 6.79;                               % luttinger parameter #1
g_2 = 1.924;                             % luttinger parameter #1
g_3 = 2.782;                             % luttinger parameter #1
delta = 0.341;                       % spin orbit coupling energy [eV]
m_c = 0.065;                             % conduction band mass [m_0]
E_v = 0;                                    % valence band initial energy [eV]
E_c = E_v+E_g;                    % conduction band initial energy [eV]
c = 3e8;                                     % speed of light [m/sec]
m_0 = 0.51e6/c^2;                    % free electron mass [eV*sec^2/m^2]
h = 4.135e-15;                           % Planck's constant [eV*sec]
h_bar = h/(2*pi);                   % h_bar [eV*sec]
h_m = h_bar^2/m_0;              % [eV*m^2]
a = 5.6532e-10;                         % latice const [m]

%%%%% initial calculation %%%%%
p = sqrt(h_m*E_p/2);            
A = 0.5*h_m/m_c-(p^2/E_g)*(E_g+(2/3)*delta)/(E_g+delta);
B = 0;
N = -3*h_m*g_3+p^2/E_g;
M = -0.5*h_m*(g_1-2*g_2);
L = -0.5*h_m*(g_1+4*g_2)+p^2/E_g;

%%%%%% k values %%%%%%%
k= (-0.1*2*pi/a):(0.01*2*pi/a):(0.1*2*pi/a);
%%%% (1,0,0) direction %%%%
k_x = k;
k_y = 0*k_x;
k_z = 0*k_x;
%%%% (1,1,0) direction %%%%
% k_x = 1/sqrt(2)*k;
% k_y = k_x;
% k_z = 0*k_x;
%%%% (1,1,0) direction %%%%
% k_x = 1/sqrt(3)*k;
% k_y = k_x;
% k_z = k_x;
%%%% (1,1,1) direction %%%%
% k_x = 1/sqrt(14)*k;
% k_y = 2*k_x;
% k_z = 3*k_x;

%%% the matrices %%%
d=zeros(8,length(k));
for r=1:length(k)
G_1 = [         E_c               i*p*k_x(r)             i*p*k_y(r)               i*p*k_z(r)    ;...
            -i*p*k_x(r)        E_v-delta/3                     0                                  0            ;...
            -i*p*k_y(r)                   0                      E_v-delta/3                      0            ;...
            -i*p*k_z(r)                   0                                  0                      E_v-delta/3];
        
G_2 = [       A*k(r)^2                          B*k_y(r)*k_z(r)                            B*k_x(r)*k_z(r)                              B*k_y(r)*k_x(r);...
           B*k_y(r)*k_z(r)   L*k_x(r)^2+M*(k_y(r)^2+k_z(r)^2)           N*k_x(r)*k_y(r)                              N*k_x(r)*k_z(r);...
           B*k_x(r)*k_z(r)                    N*k_x(r)*k_y(r)         L*k_y(r)^2+M*(k_x(r)^2+k_z(r)^2)               N*k_y(r)*k_z(r);...
           B*k_y(r)*k_x(r)                    N*k_x(r)*k_z(r)                            N*k_z(r)*k_y(r)       L*k_z(r)^2+M*(k_x(r)^2+k_y(r)^2)
           ];
G_so = -(delta/3)*[0   0    0    0;...
                                      0   0    i   0;...
                                      0  -i   0    0;...
                                      0   0    0    0];
                                
Gamma = -(delta/3)*[0   0   0    0;...
                                        0   0   0   -1;...
                                        0   0   0    i;...
                                        0   1  -i   0];
                                    
G = G_1+G_2+G_so;
H = zeros(8,8);
H(1:4,1:4) = G;
H(5:8,5:8) = conj(G);
H(1:4,5:8) = Gamma;
H(5:8,1:4) = -conj(Gamma);
d(:,r) = real(eig(H));
end
d=sort(real(d));
k=a/(2*pi)*k;
k_x=a/(2*pi)*k_x;
figure;
plot (k,d(7,:),'m','LineWidth',1.5);
hold on;
plot (k,d(5,:),'g','LineWidth',1.5);
hold on;
plot (k,d(3,:),'b','LineWidth',1.5);
hold on;
plot(k,d(1,:),'r','LineWidth',1.5);
xlim([ -0.13 0.13]);
title('\bfConduction and Valence band energy dispersion, k vector at [1 0 0] direction');
xlabel('\bfka /2\pi');
ylabel('\bfEnergy [eV]');
legend('electrons','heavy holes','light holes','SO holes');
grid;