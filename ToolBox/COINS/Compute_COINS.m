%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SPECTRUM: A MATLAB Toolbox for Top-down Proteomics     %
%                           Version 2.0.0                          %
%        Copyright (c) Biomedical Informatics Research Laboratory, %
%          Lahore University of Management Sciences Lahore (LUMS), %
%                           Pakistan.                              %
%                (http://biolabs.lums.edu.pk/BIRL)                 %
%                    (safee.ullah@gmail.com)                       %
%                 Last Modified on: 25-October-2020                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FPA new modified 05042020

function FinalImp_Data = Compute_COINS(Imported_Data)
%COMPUTE_FPA creates complementary ions for flipped peak spectra

FPA_tolerance = getappdata(0,'FPA_Tol');  %user dependent tolerance for Flipped peak spectra
if (ischar(FPA_tolerance))
    FPA_tolerance = str2num(FPA_tolerance);
end
FPATolIndex = getappdata(0,'FPAtolUnit');
tolerance_for_Flipped_Peaks =[];
if(FPATolIndex == 1)
    FPAtolUnit = 'Da';
    tolerance_for_Flipped_Peaks = FPA_tolerance;
elseif (FPATolIndex == 2)
    FPAtolUnit = 'mmu';
    tolerance_for_Flipped_Peaks = FPA_tolerance;
elseif (FPATolIndex == 3)
    FPAtolUnit = 'ppm';
    tolerance_for_Flipped_Peaks = FPA_tolerance/1000000;
else
    FPAtolUnit = '%';
    tolerance_for_Flipped_Peaks = FPA_tolerance/1000000;
end

FPA_percentageallowed = getappdata(0,'FPA_percentage');  %user dependent tolerance for Flipped peak spectra
percentage_for_Flipped_Peaks =[];
if(FPA_percentageallowed == 1)
    percentage_for_Flipped_Peaks = 10;
elseif (FPA_percentageallowed == 2)
    percentage_for_Flipped_Peaks = 25;
elseif (FPA_percentageallowed == 3)
    percentage_for_Flipped_Peaks = 50;
elseif (FPA_percentageallowed == 4)
    percentage_for_Flipped_Peaks = 75;
else (FPA_percentageallowed == 5)
    percentage_for_Flipped_Peaks = 100;
end


com_ion = [Imported_Data(:,1)]; %separates column containing masses in imported data
com_int = [Imported_Data(:,2)]; %separates column containing intensities in imported data
comp_int = [];

Peaklist_size = size(Imported_Data,1); %get peaklist size
if Peaklist_size > 5
%%For extracting peaks from the center of the peaklist

    %% 10% flipped peaks from the center allowed only
    if percentage_for_Flipped_Peaks == 10
        Peaklist_size = size(Imported_Data,1); %get peaklist size
        Mid_45_percent= (45/100)*Peaklist_size; %find 37.5% of peaklist to get 25% of data from middle of peaklist
        Mid_45_percent_roundoff = round(Mid_45_percent); %roundoff to later find index to iterate through peaklist
        if Mid_45_percent_roundoff == 0
            Mid_45_percent_roundoff = 1;
        end
        Mid_55_percent= Peaklist_size-Mid_45_percent_roundoff; %find 75% of peaklist to get 50% of data from middle of peaklist
        
        Extract_New_Peaklist =[];
        for j=drange(Mid_45_percent_roundoff:Mid_55_percent)
            Extract_values = Imported_Data(j);
            Extract_ints = Imported_Data(j,2);
            Extract_New_Peaklist = [Extract_New_Peaklist; Extract_values,Extract_ints];
            Extract_values =[];
        end
        
        %% 25% flipped peaks from the center allowed only
    elseif percentage_for_Flipped_Peaks == 25
        Peaklist_size = size(Imported_Data,1); %get peaklist size
        Mid_375_percent= (37.5/100)*Peaklist_size; %find 37.5% of peaklist to get 25% of data from middle of peaklist
        Mid_375_percent_roundoff = round(Mid_375_percent); %roundoff to later find index to iterate through peaklist
        if Mid_375_percent_roundoff == 0
            Mid_375_percent_roundoff = 1;
        end
        Mid_625_percent= Peaklist_size-Mid_375_percent_roundoff; %find 75% of peaklist to get 50% of data from middle of peaklist
        
        Extract_New_Peaklist =[];
        for j=drange(Mid_375_percent_roundoff:Mid_625_percent)
            Extract_values = Imported_Data(j);
            Extract_ints = Imported_Data(j,2);
            Extract_New_Peaklist = [Extract_New_Peaklist; Extract_values,Extract_ints];
            Extract_values =[];
        end
        
        %50% flipped peaks from the center allowed only
    elseif percentage_for_Flipped_Peaks == 50
        Peaklist_size = size(Imported_Data,1); %get peaklist size
        Mid_025_percent= (25/100)*Peaklist_size; %find 25% of peaklist to get 50% of data from middle of peaklist
        Mid_025_percent_roundoff = round(Mid_025_percent); %roundoff to later find index to iterate through peaklist
        if Mid_025_percent_roundoff == 0
            Mid_025_percent_roundoff = 1;
        end
        Mid_075_percent= Peaklist_size-Mid_025_percent_roundoff; %find 75% of peaklist to get 50% of data from middle of peaklist
        
        Extract_New_Peaklist =[];
        for j=drange(Mid_025_percent_roundoff:Mid_075_percent)
            Extract_values = Imported_Data(j);
            Extract_ints = Imported_Data(j,2);
            Extract_New_Peaklist = [Extract_New_Peaklist; Extract_values,Extract_ints];
            Extract_values =[];
        end
        
        %75% flipped peaks from the center allowed only
    elseif percentage_for_Flipped_Peaks == 75
        Peaklist_size = size(Imported_Data,1); %get peaklist size
        Mid_0125_percent= (12.5/100)*Peaklist_size; %find 12.5% of peaklist to get 75% of data from middle of peaklist
        Mid_0125_percent_roundoff = round(Mid_0125_percent); %roundoff to later find index to iterate through peaklist
        if Mid_0125_percent_roundoff == 0
            Mid_0125_percent_roundoff = 1;
        end
        Mid_0875_percent= Peaklist_size-Mid_0125_percent_roundoff; %find 87.5% of peaklist to get 75% of data from middle of peaklist
        
        Extract_New_Peaklist =[];
        for j=drange(Mid_0125_percent_roundoff:Mid_0875_percent)
            Extract_values = Imported_Data(j,:);
            Extract_New_Peaklist = [Extract_New_Peaklist; Extract_values];
            Extract_values =[];
        end
        
        %100% all flipped peaks allowed
    else percentage_for_Flipped_Peaks == 100
        Peaklist_size = size(Imported_Data,1); %get peaklist size
        All_ions= (100/100)*Peaklist_size; %find 100% of peaklist to get all the data from the peaklist
        All_ions_roundoff = round(All_ions); %roundoff to later find index to iterate through peaklist
        if All_ions_roundoff == 0
            All_ions_roundoff = 1;
        end
        
        Extract_New_Peaklist =[];
        for j=2:All_ions_roundoff
            Extract_values = Imported_Data(j);
            Extract_ints = Imported_Data(j,2);
            Extract_New_Peaklist = [Extract_New_Peaklist; Extract_values,Extract_ints];
            Extract_values =[];
        end
    end
    
    %Using user set value for percentage of flipped peaks allowed to generate
    %complementary ions
    find_comp_ion_1_tol = [];
    find_comp_ion_2_tol = [];
    
    Fragmentation_Types_Selected = getappdata(0,'Fragmentation'); %find user-specified fragmentation method to compute complementary ion
    similar_compion = {'CID','IMD','BIRD','SID','HCD'}; %same kind of fragments for these techniques ie b and y
    similar_compion2 = {'ECD','ETD'}; %same kind of fragments for these techniques ie c and z
    similar_compion3 = {'EDD','NETD'}; %same kind of fragments for these techniques ie a and x
    data = getappdata(0,'HandleIon'); %find special ions
    match_fragtechnique = strcmpi(Fragmentation_Types_Selected,similar_compion); % compares the selected fragmentation technique with the set of similar techniques (ions wise)
    match_fragtechnique2 = strcmpi(Fragmentation_Types_Selected,similar_compion2); % compares the selected fragmentation technique with the set of similar techniques (ions wise)
    match_fragtechnique3 = strcmpi(Fragmentation_Types_Selected,similar_compion3); % compares the selected fragmentation technique with the set of similar techniques (ions wise)
    %masses to be added
    H = 1.007825035;
    C = 12.0000;
    N_term = 1.01;
    C_term = 17.01;
    O = 15.99491463;
    N = 14.003074;
    CO = C + O;
    NH = N + H;
    
    for i = 1 : size(Extract_New_Peaklist,1)
        
        %%find matches with for b and y ions
        if any(match_fragtechnique) %b and y ions
            % Formula: b_ion = Rf +<N>
            % y_ion = Rf + C + 2H
            Rf = Extract_New_Peaklist(i,1) - N_term ;
                    comp_y_ion = Rf + C_term + H + H ;
            
            floor_comp_y_ion = floor(comp_y_ion);
            
            find_comp_y_ion = ismember(floor_comp_y_ion-tolerance_for_Flipped_Peaks : floor_comp_y_ion+tolerance_for_Flipped_Peaks, com_ion); %find if computed ion already exists within user tolerance in the experimental peaklist
            if find_comp_y_ion == 1 %if computed ion already exists in the experimental peaklist
                comp_ion =[]; %do not append computed ion
                comp_int =[]; % do not append corresponding intensity
            elseif find_comp_y_ion == 0
                com_ion = [com_ion; comp_y_ion]; %append if computed ion does not exist in the experimental peaklist
                comp_int = [Extract_New_Peaklist(i,2)]; %replicates same intensity for complementary ion as its counterpart
                com_int = [com_int; comp_int];
            end
            comp_ion =[];
            comp_int =[];
            %% c and z ions
        else
            if any(match_fragtechnique2) %c and z ions
                % Formula: c_ion = Rf +<N> + N + 3H
                % z_ion = Rf + C - NH
                Rf = Extract_New_Peaklist(i,1) - N_term - N - (H + H+ H);
                comp_z_ion = Rf + C_term - (NH);
                floor_comp_z_ion = floor(comp_z_ion);
                
                find_comp_z_ion = ismember(floor_comp_z_ion-tolerance_for_Flipped_Peaks : floor_comp_z_ion+tolerance_for_Flipped_Peaks, com_ion); %find if computed ion already exists within user tolerance in the experimental peaklist
                if find_comp_z_ion == 1 %if computed ion already exists in the experimental peaklist
                    comp_ion =[]; %do not append computed ion
                    comp_int =[]; % do not append corresponding intensity
                elseif find_comp_z_ion == 0
                    com_ion = [com_ion; comp_z_ion]; %append if computed ion does not exist in the experimental peaklist
                    comp_int = [Extract_New_Peaklist(i,2)]; %replicates same intensity for complementary ion as its counterpart
                    com_int = [com_int; comp_int];
                end
                comp_ion =[];
                comp_int =[];
                
                %% a and x ions
            elseif any(match_fragtechnique3) %c and z ions
                % Formula: a_ion = Rf + <N> - CO
                % x_ion = Rf + C + CO
                Rf = Extract_New_Peaklist(i,1) - N_term + CO ;
                comp_x_ion = Rf + C_term + CO ;
                
                floor_comp_x_ion = floor(comp_x_ion);
                
                find_comp_x_ion = ismember(floor_comp_x_ion-tolerance_for_Flipped_Peaks : floor_comp_x_ion+tolerance_for_Flipped_Peaks, com_ion); %find if computed ion already exists within user tolerance in the experimental peaklist
                if find_comp_x_ion == 1 %if computed ion already exists in the experimental peaklist
                    comp_ion =[]; %do not append computed ion
                    comp_int =[]; % do not append corresponding intensity
                elseif find_comp_x_ion == 0
                    com_ion = [com_ion; comp_x_ion]; %append if computed ion does not exist in the experimental peaklist
                    comp_int = [Extract_New_Peaklist(i,2)]; %replicates same intensity for complementary ion as its counterpart
                    com_int = [com_int; comp_int];
                end
                comp_ion =[];
                comp_int =[];
            end
        end
    end
    
    %%
    %Append duplicated intensities of complementary ions to main peak list
    ImportedData = [];
    Imported_Data = [com_ion, com_int];
    Sortcomp = sortrows(Imported_Data, [1,2]); %sort new peaklist
    n=size(Sortcomp,1);
    MS1_data = Sortcomp(n,:); %keeping MS1 on top for later search
    MS2_data = Sortcomp(1:n-1,:);
    ImportedData = [MS1_data];
    ImportedData = [ImportedData; MS2_data]; %final imported data
    
    %replacing negative numbers by zeros
    rowtodelete =[];
    rowstodelete = [];
    for x=1:size(ImportedData,1)
        % for x=1:numel(ImportedData)
        if ImportedData(x)<0 || ImportedData(x)<50
            ImportedData(x) =0;
            if ImportedData(x) == 0
                ImportedData(x,:) = [];
            end
        end
        if x == size(ImportedData,1)
            break
        end
    end
    FinalImp_Data = ImportedData;

else
    FinalImp_Data = Imported_Data;
end

pwd = fileparts(mfilename('fullpath'));
[filepath, folder] = fileparts(pwd);
cd(strcat(filepath));
return;
end
