function[k_list,k]=kZB_f(Nk,a)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Brillouin zone vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% (100) or kx-direction
kx=linspace(-1,0,Nk);
k1 = [ kx'  0*kx'  0*kx' ];

%%% (001) or kz-direction
kz=linspace(0,1,Nk);
k2 = [ 0*kz'  0*kz'  1*kz' ];


k_list=[ k1 ; k2  ]*2*pi/a*0.2;
%k_ZB=[ k1 ; k2 ; k3 ; k4 ]*2*pi/a;

k=sqrt( k_list(:,1).^2 + k_list(:,2).^2 + k_list(:,3).^2 );
k=[ -k(1:Nk) ; k(Nk+1:end)  ];

%Nk=50;

%%%% L-valley
%kL=linspace(-0.5,0,Nk);
%k1 = [ kL'  kL'  kL' ];
%%%% X-valley
%kx=linspace(0,1,Nk);
%k2 = [ kx'  0*kx'  0*kx' ];
%%%% X point to U point
%ky=linspace(0,0.25,floor(Nk/2));
%kz=linspace(0,0.25,floor(Nk/2));
%k3 = [ ky'*0+1  ky'  kz' ];
%%%% U point to Gamma point
%kx=linspace(1,0,Nk);
%ky=linspace(0.25,0,Nk);
%k4 = [ kx'  ky'  ky' ];
%%%% here, it just to checked that it is isotrope... and make the last plot
%kx=linspace(0,1,Nk);
%k5 = [
%kx'*1  kx'*0  kx'*0
%kx'*0  kx'*1  kx'*0
%kx'*0  kx'*0  kx'*1
%];
%
%k_list=[ k1 ; k2 ; k3 ; k4 ; k5]*2*pi/a;
%%k_list=[ k4(1,:)]*2*pi/a;
%%k_list=[k_list(1,1) 0 0];
%%k=sqrt( k_list(:,1).^2 + k_list(:,2).^2 + k_list(:,3).^2 );

end