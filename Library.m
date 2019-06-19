%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Library load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB = importdata('materialDB.csv',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M{1}=Material;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here is a patch to be able to load the tables in Matlab AND Octave %%%%%%
% Matlab see the header in multiple cells while Octave see the header in one cell only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB.textdata(1,:))==1  %% Octave data load

    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);
    
    M{2} = DB.data(:,idxM);
else  %% Matlab data load

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M{1});
      if idx==1
        M{2}=DB.data(:,i-1);
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%