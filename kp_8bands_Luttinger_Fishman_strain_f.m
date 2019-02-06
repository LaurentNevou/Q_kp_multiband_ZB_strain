function[E]=kp_8bands_Luttinger_Fishman_strain_f(k_list, Eg, EP, Dso, F, g123, ac, av, bv, dv, exx, ezz)

% Guy Fishman
% "Semi-Conducteurs: les Bases de la Theorie k.p " (2010)
% 3.3.4 L’hamiltonien projeté sur {G6 ; G8 ; G7} : l’hamiltonien H8 de Pidgeon-Brown
% page 169
% https://www.amazon.fr/Semi-Conducteurs-Bases-Theorie-K-P-Fishman/dp/2730214976/ref=sr_1_fkmr1_1?ie=UTF8&qid=1548234034&sr=8-1-fkmr1&keywords=guy+fishman+kp
% https://www.abebooks.fr/semi-conducteurs-bases-th%C3%A9orie-k.p-Fishman-ECOLE/30091636895/bd
% https://www.decitre.fr/livres/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971.html
% https://www.unitheque.com/Livre/ecole_polytechnique/Semi_conducteurs_les_bases_de_la_theorie_K.p-35055.html
% https://www.eyrolles.com/Sciences/Livre/semi-conducteurs-les-bases-de-la-theorie-k-p-9782730214971/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% charge de l electron [Coulomb]
m0=9.10938188E-31;              %% electron mass [kg]
H0=hbar^2/(2*m0) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dso = Dso*e;
Eg  = Eg*e;
EP  = EP*e;
P   = sqrt(EP*hbar^2/(2*m0)) ;

% gc= 1+2*F + EP*(Eg+2*Dso/3) / (Eg*(Eg+Dso)) ;   % =1/mc  electron in CB eff mass

% renormalization of the paramter from 6x6kp to 8x8kp
% gc=gc-EP/3*( 2/Eg + 1/(Eg+Dso) );

gc = 1+2*F;
g1 = g123(1)-EP/(3*Eg);
g2 = g123(2)-EP/(6*Eg);
g3 = g123(3)-EP/(6*Eg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fishman is using the Pidgeon Brown Hamiltonian...
% It uses 3 additionnal parameters
% If we set thoses parameter = Luttinger parameters, the results are exactly 
% the same as for the other models like Pistol or DKK

%gD1=gD1-EP/(3*(Eg+Dso));
%gD2=gD2-EP/12*(1/Eg+1/(Eg+Dso));
%gD3=gD3-EP/12*(1/Eg+1/(Eg+Dso));

gD1=g1;
gD2=g2;
gD3=g3;

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

Hs=[

 ac*ee   0           0                 0                 0                 0             0             0
   0   ac*ee         0                 0                 0                 0             0             0
   0     0   av*ee-bv*(exx-ezz)        0                 0                 0             0             0
   0     0           0         av*ee+bv*(exx-ezz)        0                 0    sqrt(2)*bv*(exx-ezz)   0
   0     0           0                 0         av*ee+bv*(exx-ezz)        0             0   -sqrt(2)*bv*(exx-ezz)
   0     0           0                 0                 0         av*ee-bv*(exx-ezz)    0             0
   0     0           0        sqrt(2)*bv*(exx-ezz)       0                 0           av*ee           0
   0     0           0                 0       -sqrt(2)*bv*(exx-ezz)       0             0           av*ee

];

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
kpp =  kx + 1i*ky;
kmm =  kx - 1i*ky;

AA = H0*g2*( 2*kz^2 - kx.^2 - ky.^2 );
BB = H0*2*sqrt(3)*g3*kz*(kx - 1i*ky) ;
CC = H0*sqrt(3)*(g2*(kx^2-ky^2)-2i*g3*kx*ky);

AA_D = H0*gD2*( 2*kz^2 - kx.^2 - ky.^2 );
BB_D = H0*2*sqrt(3)*gD3*kz*(kx - 1i*ky) ;
CC_D = H0*sqrt(3)*(gD2*(kx^2-ky^2)-2i*gD3*kx*ky);

Hdiag = -H0*k^2*[-gc -gc g1 g1 g1 g1 gD1 gD1] + [ +Eg +Eg +AA  -AA  -AA  +AA  -Dso  -Dso ];

% Ec- Ec+     HH+                 LH+               LH-              HH-            SO+              SO-
H=[
  0   0  -sqrt(1/2)*P*kpp   sqrt(2/3)*P*kz    sqrt(1/6)*P*kmm         0         sqrt(1/3)*P*kz    sqrt(1/3)*P*kmm  % Ec+
  0   0        0           -sqrt(1/6)*P*kpp   sqrt(2/3)*P*kz   sqrt(1/2)*P*kmm  sqrt(1/3)*P*kpp  -sqrt(1/3)*P*kz   % Ec-
  0   0        0                  BB                CC                0         sqrt(1/2)*BB_D    sqrt(2)  *CC_D   % HH+
  0   0        0                   0                 0               CC        -sqrt(2)  *AA_D   -sqrt(3/2)*BB_D   % LH+
  0   0        0                   0                 0              -BB        -sqrt(3/2)*BB_D'   sqrt(2)  *AA_D   % LH-
  0   0        0                   0                 0                0        -sqrt(2)  *CC_D'   sqrt(1/2)*BB_D'  % HH-
  0   0        0                   0                 0                0              0                0            % SO+
  0   0        0                   0                 0                0              0                0            % SO-
];

H=H'+H+diag(Hdiag);
H=H+Hs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end