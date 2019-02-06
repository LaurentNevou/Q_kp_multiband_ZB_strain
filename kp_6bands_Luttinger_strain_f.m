function[E]=kp_6bands_Luttinger_strain_f(k_list, Dso, g123, av, bv, dv, exx, ezz);

% Calvin Yi-Ping Chao and Shun Lien Chuang
% PHYSICAL REVIEW B VOLUME 46, NUMBER 7, (1992)
% "Spin-orbit-coupling effects on the valence-band structure of strained semiconductor quantum wells"

% Warren James Elder
% PhD thesis: "Semi-empirical modelling of SiGe hetero-structures" (2012)
% Chapter3: "The k.p method" page 68
% https://spiral.imperial.ac.uk:8443/handle/10044/1/9816

% Maxi Scheinert
% PhD thesis: "Optical Pumping: A Possible Approach towards a SiGe Quantum Cascade Laser"
% page 18
% https://core.ac.uk/download/pdf/20642119.pdf
% https://www.e-helvetica.nb.admin.ch/directAccess?callnumber=rero-004-104096

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;
Dso=Dso*e;
g1=g123(1);
g2=g123(2);
g3=g123(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
exy = 0;
ezx = 0;
eyz = 0;
ee  = exx+eyy+ezz;
av  = abs(av)*e;
bv  = abs(bv)*e;
dv  = abs(dv)*e;

Pe =  av * ee;
Qe =  -bv * (exx-ezz);
Re =  sqrt(3)/2 * bv*(exx-eyy) - 1i*dv*exy;
Se = -dv * (ezx - 1i*eyz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k+ =  kx + 1i*ky
%k- =  kx - 1i*ky

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);

Pk = -H0 * g1 * k^2;
Qk = -H0 * g2 *(kx^2 + ky^2 - 2*kz^2);
Rk =  H0 * sqrt(3)  * (g2*(kx^2-ky^2) - 2i*g3*kx*ky );
Sk =  H0 *2*sqrt(3) * g3*(kx-1i*ky)*kz;

P = Pk + Pe;
Q = Qk + Qe;
R = Rk + Re;
S = Sk + Se;

Hdiag = [ P+Q   P-Q   P-Q   P+Q   -Dso+P   -Dso+P  ];

%   HH     LH       LH        HH           SO              SO

H=[
    0      -S        R         0     -sqrt(1/2)*S       sqrt(2)  *R  % HH
    0       0        0         R     -sqrt(2)  *Q       sqrt(3/2)*S  % LH
    0       0        0         S      sqrt(3/2)*S'      sqrt(2)  *Q  % LH
    0       0        0         0     -sqrt(2)  *R'     -sqrt(1/2)*S' % HH
    0       0        0         0             0              0        % SO
    0       0        0         0             0              0        % SO
];

H=H'+H+diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end