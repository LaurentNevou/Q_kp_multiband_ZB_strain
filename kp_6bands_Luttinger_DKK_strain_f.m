function[E]=kp_6bands_Luttinger_DKK_strain_f(k_list, Dso, g123, av, bv, dv, exx, ezz)

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
Dso = Dso*e;
g1=g123(1);
g2=g123(2);
g3=g123(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eyy = exx;
exy = 0; eyx=0;
ezx = 0; exz=0;
eyz = 0; ezy=0;
ee  = exx+eyy+ezz;
av  = abs(av)*e;
bv  = abs(bv)*e;
dv  = abs(dv)*e;

l = av-2*bv;
m = av+bv;
n  = sqrt(3)*dv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Hdiag = H0*k^2 -Dso/3*ones(1,6);

HR=[
   L*kx^2+M*(ky^2+kz^2)          N*kx*ky                N*kx*kz
         N*kx*ky          L*ky^2+M*(kx^2+kz^2)          N*ky*kz
         N*kx*kz                 N*ky*kz         L*kz^2+M*(kx^2+ky^2)
];

Hso=[
   0   1   0     0   0   1i
  -1   0   0     0   0   1
   0   0   0    -1i -1   0
   0   0  -1i    0  -1   0
   0   0   1     1   0   0
   1i -1   0     0   0   0
];

Hs=[
    l*exx+m*(eyy+ezz)        n*exy             n*exz
         n*exy          l*eyy+m*(exx+ezz)      n*eyz
         n*exz               n*eyz        l*ezz+m*(exx+eyy) 
];


H = diag(Hdiag) + Hso*Dso/(3i) + [HR  zeros(3,3) ; zeros(3,3)  HR] + [Hs  zeros(3,3) ; zeros(3,3)  Hs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e ;

end

end