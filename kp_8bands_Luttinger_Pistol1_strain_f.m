function[E]=kp_8bands_Luttinger_Pistol1_strain_f(k_list, Eg, EP, Dso, F, g123, ac, av, bv, dv, exx, ezz)

% Craig E. Pryor and M. E. Pistol
% "Atomistic k.p theory", Journal of Applied Physics 118, 225702 (2015); 
% https://doi.org/10.1063/1.4936170
% https://www.researchgate.net/publication/273067732_Atomistic_kp_theory

% Thomas B. Bahder,
% "Eight-band k.p model of strained zinc-blende crystals", PRB 41, 11992 (1990)
% https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.11992
% https://www.researchgate.net/publication/235532200_Eight-band_k_p_model_of_strained_zinc-blende_crystals

% Guy Fishman
% "Semi-Conducteurs: les Bases de la Theorie k.p " (2010)
% 4.5.4 L’hamiltonien de Bir-Pikus
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
P0  = sqrt(EP*hbar^2/(2*m0)) ; % Here I use "P0" instead of "P" because it uses "P" inside the H for something else

% gc= 1+2*F + EP*(Eg+2*Dso/3) / (Eg*(Eg+Dso)) ;   % =1/mc  electron in CB eff mass

% renormalization of the paramter from 6x6kp to 8x8kp
% gc=gc-EP/3*( 2/Eg + 1/(Eg+Dso) );

gc = 1+2*F;
g1 = g123(1)-EP/(3*Eg);
g2 = g123(2)-EP/(6*Eg);
g3 = g123(3)-EP/(6*Eg);

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
   0     0           0         av*ee+bv*(exx-ezz)        0                 0             0   -sqrt(2)*bv*(exx-ezz)
   0     0           0                 0         av*ee+bv*(exx-ezz)        0     sqrt(2)*bv*(exx-ezz)  0
   0     0           0                 0                 0         av*ee-bv*(exx-ezz)    0             0
   0     0           0                 0        sqrt(2)*bv*(exx-ezz)       0           av*ee           0
   0     0           0       -sqrt(2)*bv*(exx-ezz)       0                 0             0           av*ee

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

B =  0;
A =  Eg + gc*H0 * k^2;
U =  sqrt(1/3) * P0 * kz;
V =  sqrt(1/6) * P0 * kmm;
W =  1i*sqrt(1/3)*B * kx*ky;
T =  sqrt(1/6)*B * kz*kpp;
P = -g1*H0 * k^2;
Q = -g2*H0 *(kx^2 + ky^2 - 2*kz^2);
R = -H0 * sqrt(3)  * (g2*(kx^2-ky^2) - 2i*g3*kx*ky );
S =  H0 *2*sqrt(3) * g3*(kx-1i*ky)*kz;


% Ec- Ec+     HH+               LH+             LH-            HH-            SO-            SO+

Hdiag = [A A P+Q P-Q P-Q P+Q -Dso+P -Dso+P];

H=[
  0   0        0               T'+V'       sqrt(2)*(W-U)  -sqrt(3)*(T-V)     (W-U)     sqrt(2)*(T'+V')  % Ec-
  0   0  -sqrt(3)*(T'+V')  sqrt(2)*(W-U)       (T-V)            0       -sqrt(2)*(T-V)      W'+U        % Ec+
  0   0        0                -S               R              0       -sqrt(2)  *R     sqrt(1/2)*S    % HH+
  0   0        0                 0               0              R        sqrt(3/2)*S     sqrt(2)  *Q    % LH+
  0   0        0                 0               0              S       -sqrt(2)  *Q     sqrt(3/2)*S'   % LH-
  0   0        0                 0               0              0        sqrt(1/2)*S'    sqrt(2)  *R'   % HH-
  0   0        0                 0               0              0              0              0         % SO-
  0   0        0                 0               0              0              0              0         % SO+
];

H=H'+H+diag(Hdiag);
H=H+Hs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E(:,i) = eig(H)/e;

end

end