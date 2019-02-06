%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% last update 05Feb2019, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Here, you have to choose your material among the following %%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material='AlAs';
Material='GaAs';
%Material='InAs';
%Material='AlSb';
%Material='GaSb';
%Material='InSb';
%Material='AlP';
%Material='GaP';
%Material='InP';
%Material='AlN';
%Material='GaN';
%Material='InN';
%Material='Si';
%Material='Ge';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=300;                  % Temperature [Kelvin]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library
ExtractParameters

Nk=100;                   %% number of k-points for the dispersion
[k_list,k]=kZB_f(Nk,a);   %% function to compute the k-vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Activate the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kp_6x6_Luttinger_strain  = 1;
kp_8x8_Luttinger_strain  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add the strain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exx=0.01; %% exx = (a0-a)/a0
ezz=-2*c12/c11*exx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0;FS=20;
c=[
1 0 0
0 0 1
0 1 0
1 0 1
0 1 1
1 1 0
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kp_6x6_Luttinger_strain == 1;
  i=i+1;
  E{i}=kp_6bands_Luttinger_DKK_strain_f(k_list, Dso, g123, av, bv, dv, exx, ezz);
  %E{i}=kp_6bands_Luttinger_strain_f(k_list, Dso, g123, av, bv, dv, exx, ezz);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 6x6 bands');
end

if kp_8x8_Luttinger_strain == 1;
  i=i+1;
  E{i}=kp_8bands_Luttinger_DKK_strain_f(k_list, Eg, EP_L, Dso, F, g123, ac, av, bv, dv, exx, ezz);
  %E{i}=kp_8bands_Luttinger_Pistol1_strain_f(k_list, Eg, EP_L, Dso, F, g123, ac, av, bv, dv, exx, ezz);
  %E{i}=kp_8bands_Luttinger_Pistol2_strain_f(k_list, Eg, EP_L, Dso, F, g123, ac, av, bv, dv, exx, ezz);
  %E{i}=kp_8bands_Luttinger_Fishman_strain_f(k_list, Eg, EP_L, Dso, F, g123, ac, av, bv, dv, exx, ezz);
  s{i}=strcat('\fontsize{',num2str(FS),'}\color[rgb]{',num2str(c(i,:)),'}k.p 8x8 bands');
  
  EE=sort(E{i}(:,Nk));
  EEg=EE(8)-EE(6);
  s{i+1}=strcat('\fontsize{',num2str(FS),'}\color{black}Eg=',num2str(EEg,'%.3f'),'eV');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[50 100 1000 900]);

FS=20;
LW=2;

xscale=[-2 2];
yscale=[-1 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on;

xlabel('Wavevector (nm-1)')
ylabel('Energy (eV)')
xlim(xscale)
ylim(yscale)

text(0.8*(xscale(2)-xscale(1))+xscale(1),0.7,s);
text(-0.03*(xscale(2)-xscale(1))+xscale(1),-0.06*(yscale(2)-yscale(1))+yscale(1),{'kx:[100]'},'fontsize',FS);
text(0.92*(xscale(2)-xscale(1))+xscale(1),-0.06*(yscale(2)-yscale(1))+yscale(1),{'kz:[001]'},'fontsize',FS);

if exx==0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}unstrained'));
elseif exx>0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}exx=',num2str(exx*100),'% tensiled'));
elseif exx<0
  title(strcat(M{1},' bandstructure @T=',num2str(T),'K, \color{red}exx=',num2str(exx*100),'% compressed'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(E)

  plot(k*1e-9,E{j},'-', 'linewidth',LW,'color',c(j,:))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%