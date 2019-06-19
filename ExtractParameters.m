%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Extract general parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lattice parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'a');
  if idx==1
    a=M{2}(i-1)*1e-10;
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Kane model

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'EP_Kane');
%  if idx==1
%    EP_K = M{2}(i-1);
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EP from the Luttinger model

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'EP_Luttinger');
  if idx==1
    EP_L = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 7c

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg7c');
%  if idx==1
%    Eg7c = M{2}(i-1);
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 8c

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg8c');
%  if idx==1
%    Eg8c = M{2}(i-1);
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6v

%for i=1:length(DB.textdata(:,1))
%  idx=strcmp(DB.textdata{i,1},'Eg6v');
%  if idx==1
%    Eg6v = M{2}(i-1);
%    % break % removing the break makes it slower but more compatible between Matlab and Octave
%  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy level of the band 6c (the band gap in a direct band gap semiconductor)

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Eg6cG');
  if idx==1
    Eg6c = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "alpha" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'alphaG');
  if idx==1
    alphaG = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandgap Temperature dependency parameter "beta" at the Gamma point

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'betaG');
  if idx==1
    betaG = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eg   = Eg6c - (alphaG*T^2) ./ (T+betaG);   %Eg = Eg0 - (a*T.^2)./(T + b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split-off band energy Delta-so

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'Dso');
  if idx==1
    Dso = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter for the electrons

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'F');
  if idx==1
    F = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma1 

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma1');
  if idx==1
    g1 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma2

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma2');
  if idx==1
    g2 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luttinger parameter gamma3

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'gamma3');
  if idx==1
    g3 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g123 = [g1 g2 g3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ac Conduction band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'ac');
  if idx==1
    ac = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% av Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'av');
  if idx==1
    av = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bv Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'bv');
  if idx==1
    bv = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dv Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'dv');
  if idx==1
    dv = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c11 Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'c11');
  if idx==1
    c11 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c12 Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'c12');
  if idx==1
    c12 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c44 Valence band Strain parameter

for i=1:length(DB.textdata(:,1))
  idx=strcmp(DB.textdata{i,1},'c44');
  if idx==1
    c44 = M{2}(i-1);
    % break % removing the break makes it slower but more compatible between Matlab and Octave
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%