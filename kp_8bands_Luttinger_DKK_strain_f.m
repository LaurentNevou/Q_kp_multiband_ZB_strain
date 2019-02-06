function[E]=kp_8bands_Luttinger_DKK_strain_f(k_list, Eg, EP, Dso, F, g123, ac, av, bv, dv, exx, ezz)

% DKK model: Dresselhaus, Kip and Kittel

% Calin Galeriu
% PhD thesis: "k.p Theory of semiconductor nanostructures" (2005)
% Chapter 3, page 26
% Download:
% https://web.wpi.edu/Pubs/ETD/Available/etd-120905-095359/unrestricted/cgaleriu.pdf

% Stefan Birner (Nextnano)
% PhD thesis: "Modeling of semiconductor nanostructures and semiconductor-electrolyte interfaces" (2011)
% Chapter3, page 36: "Multi-band k.p envelope function approximation"
% Download:
% https://mediatum.ub.tum.de/doc/1084806/1084806.pdf
% https://www.nextnano.com/downloads/publications/PhD_thesis_Stefan_Birner_TUM_2011_WSIBook.pdf

% Thomas B. Bahder,
% "Eight-band k.p model of strained zinc-blende crystals", PRB 41, 11992 (1990)
% https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.11992
% https://www.researchgate.net/publication/235532200_Eight-band_k_p_model_of_strained_zinc-blende_crystals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dso = Dso*e;
Eg  = Eg*e;
EP  = EP*e;
P   = sqrt(EP*hbar^2/(2*m0)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
exy = 0; eyx=0;
ezx = 0; exz=0;
eyz = 0; ezy=0;
ee  = exx+eyy+ezz;
ac  = -abs(ac)*e;
av  = +abs(av)*e;
bv  = +abs(bv)*e;
dv  = +abs(dv)*e;

l = av-2*bv;
m = av+bv;
bb = 0;
n  = sqrt(3)*dv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gc= 1+2*F + EP*(Eg+2*Dso/3) / (Eg*(Eg+Dso)) ;   % =1/mc  electron in CB eff mass

% renormalization of the paramter from 6x6kp to 8x8kp
% gc=gc-EP/3*( 2/Eg + 1/(Eg+Dso) );

gc = 1+2*F;
g1 = g123(1)-EP/(3*Eg);
g2 = g123(2)-EP/(6*Eg);
g3 = g123(3)-EP/(6*Eg);

L =  H0*(-1-g1-4*g2);
M =  H0*(-1-g1+2*g2);
N = -H0*6*g3;
B =  H0*0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonien %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k_list(:,1))
  
kx = k_list(i,1);
ky = k_list(i,2);
kz = k_list(i,3);

k=sqrt(kx.^2 + ky.^2 + kz.^2);
Hdiag = H0*k^2*[gc 1 1 1]  +  [Eg -Dso/3 -Dso/3 -Dso/3 ];

H4=[

  0     1i*P*kx   1i*P*ky   1i*P*kz
  0        0         0         0
  0        0         0         0
  0        0         0         0

];

HH4 = H4' + H4 + diag(Hdiag); 

HR=[
  0             B*ky*kz                 B*kx*kz                B*kx*ky
B*ky*kz   L*kx^2+M*(ky^2+kz^2)          N*kx*ky                N*kx*kz
B*kx*kz         N*kx*ky          L*ky^2+M*(kx^2+kz^2)          N*ky*kz
B*kx*ky         N*kx*kz                 N*ky*kz         L*kz^2+M*(kx^2+ky^2)
];


Hso=[

 0   0   0   0   0   0   0   0
 0   0   1   0   0   0   0   1i
 0  -1   0   0   0   0   0   1
 0   0   0   0   0  -1i -1   0
 0   0   0   0   0   0   0   0
 0   0   0  -1i  0   0  -1   0
 0   0   0   1   0   1   0   0
 0   1i -1   0   0   0   0   0

];



Hs1=[

ac*ee            bb*eyz              bb*exz            bb*exy
bb*eyz      l*exx+m*(eyy+ezz)        n*exy             n*exz
bb*exz           n*exy          l*eyy+m*(exx+ezz)      n*eyz
bb*exy           n*exz               n*eyz        l*ezz+m*(exx+eyy) 

];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here are the terms that make it slightly different from Pistl-1 and Fishman Hamiltonian
%% Removing those terms gives exactly the same results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hs2=[

0   -1i*P*(exx*kx+exy*ky+exz*kz)  -1i*P*(eyx*kx+eyy*ky+eyz*kz)  -1i*P*(ezx*kx+ezy*ky+ezz*kz)
0                0                             0                             0
0                0                             0                             0
0                0                             0                             0

];

Hs=Hs1;%+Hs2+Hs2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=[HH4  zeros(4,4) ; zeros(4,4)  HH4] + Hso*Dso/(3i) + [HR  zeros(4,4) ; zeros(4,4)  HR] + [Hs  zeros(4,4) ; zeros(4,4)  Hs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end