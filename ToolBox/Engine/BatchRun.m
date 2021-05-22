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
function BatchRun()
progressbar('SPECTRUM');
data = {};
csvData = [cellstr('File Name'),cellstr('Protein Header'),cellstr('Terminal Modification'),cellstr('Protein Sequence'),cellstr('Protein Tuncation'),cellstr('Truncation Position'),cellstr('Score'),cellstr('Molecular Weight'),cellstr('No of Modifications'),cellstr('No of Fragments Matched'),cellstr('Run Time'),cellstr('E-Value'),cellstr('PST Length:1'),cellstr('PST Length:2'),cellstr('PST Length:3'),cellstr('PST Length:4'),cellstr('PST Length:5'),cellstr('PST Length:6'),cellstr('PST Length:7'),cellstr('PST Length:8'),cellstr('PST Length:9'),cellstr('PST Length:10')]; %#ok<AGROW>
setappdata(0,'P_condotion',1);
Database_Path = getappdata(0,'Database_Path'); %#ok<*NASGU> % Database path of single and batch search modes
Selected_Database = getappdata(0,'Selected_Database'); % Concatenated: Database_Path + Database selected from menu_Database
Batch_Peaklist_Data = getappdata(0,'Batch_Peaklist_Data'); % Peaklist data for BATCH search mode
DirectoryContents=getappdata(0,'BatchFolderContents');
Path=getappdata(0,'path');

BlindPTMvar = getappdata(0,'BlindPTM');
MW_Tolerance = getappdata(0,'Molecular_Weight_Tol'); % Molecular Weight tolerance
MWTolIndex = getappdata(0,'MWunit');
if(MWTolIndex == 1)
    MWunit = 'Da';
elseif (MWTolIndex == 2)
    MWunit = 'mmu';
elseif (MWTolIndex == 3)
    MWunit = 'ppm';
else
    MWunit = '%';
end
Peptide_Tolerance = getappdata(0,'Peptide_Tol'); % Peptide tolerance
FilterPSTs = getappdata(0,'FilterPSTs');
PepTolIndex = getappdata(0,'PepUnit');
if(PepTolIndex == 1)
    PepUnit = 'Da';
elseif (PepTolIndex == 2)
    PepUnit = 'mmu';
elseif (PepTolIndex == 3)
    PepUnit = 'ppm';
else
    PepUnit = '%';
end

f = filesep;
initial_path=pwd;
setappdata(0,'inital',initial_path);

PTM_Tolerance = getappdata(0,'PTM_Tol'); % Peptide tolerance
HandleIon = getappdata(0,'HandleIon');
User_EST_parameters=getappdata(0,'User_EST_parameters');
Fixed_Modifications = getappdata(0,'Fixed_Modifications'); % Fixed Modifications
Variable_Modifications = getappdata(0,'Variable_Modifications'); % Variable Modifications
Fragmentation_Type = getappdata(0,'Fragmentation_type'); % Fragmentation type
Experimental_Protein_Mass = getappdata(0,'Experimental_Protein_Mass'); % Experimental Proein Mass (Protein mass)
Project_Title=getappdata(0,'Project_Title');
progressbar('SPECTRUM: Loading Database');

Candidate_ProteinsList_Full = LoadDatabaseBatch(Selected_Database);
progressbar('SPECTRUM: Performing Search');

for b= 1:size(DirectoryContents,1)     % for mgf and text files  (all files present in specified folder)
    clc;
    batchModeTimer = tic;
    progressbar(b/(size(DirectoryContents,1)+1));
    Batch_Peaklist_Data(b).name = Batch_Peaklist_Data(b).name(1:end-4);
    %% add path
    Save_Batch_File = strcat(getappdata(0,'Result_Folder'),f,Project_Title,'_',Batch_Peaklist_Data(b).name,'.results');
    Batch_Search_File = [Path,f,DirectoryContents(b).name];
    
    if(getappdata(0,'Type_file')==1)
        
        Imported_Data =  importdata(Batch_Search_File);
       
        m = max(Imported_Data(:,2));
        HandleIon = getappdata(0,'HandleIon');
        for tupleIndex = 2 : size(Imported_Data,1)
            %% Used for converting mono isotopic mass into mz value
            if (HandleIon{9,1} == 0)
                Imported_Data(tupleIndex,1) = Imported_Data(tupleIndex,1)+1.00727647;
            end
            Imported_Data(tupleIndex,2) = Imported_Data(tupleIndex,2)/m;
        end
        Imported_Data;
        
        
        %Is FPA checked?
        FPAvalue = getappdata(0,'FPAvalue');
        if isempty(FPAvalue)
            FPAvalue =0;
        end
        %% Generate complementary ions if FPA is checked
        if FPAvalue ==1
            cd(strcat(pwd,'\COINS'));
            FinalImp_Data = Compute_COINS(Imported_Data);
%             Imported_Data = FinalImp_Data;
            
            %write new data in complementary ions folder
            CurrentWorkingDirectory = fileparts(mfilename('fullpath'));
            idcs   = strfind(CurrentWorkingDirectory,'\Engine');
            newdir = CurrentWorkingDirectory(1:idcs(end)-1);
            file_path = strcat(fullfile(newdir, 'ComplementaryIons'));
            [pathstr,name,extension] = fileparts(Batch_Search_File) ;
            f= filesep;
            file = strcat(file_path, f, name, extension);
            compfile = fopen(file,'w');
            dlmwrite(file,FinalImp_Data,'delimiter','\t','newline', 'pc');
            fclose('all');
            
        elseif FPAvalue ==0 
            Imported_Data;
        end
        
        %added 07042020
        if FPAvalue == 1 %Is FPA checked?
            Comp_Peaks_File = FinalImp_Data;
            Sorted = sortrows(FinalImp_Data,1);
            PeakListMW_Comp = Sorted(:,1);
            Intensity_Comp = Sorted(:,2);
            setappdata(0,'Comp_Peaklist_Data',FinalImp_Data);
            setappdata(0,'Comp_Fragments_Masses',PeakListMW_Comp);  % MS2s (in ascending order) then, MS1 only 
            setappdata(0,'Comp_Int',Intensity_Comp);   
            if (isempty(FinalImp_Data))
                continue
            end
            Experimental_Protein_Mass = FinalImp_Data(1,1);

        else
            Peaks_File = Imported_Data;
            Sorted = sortrows(Imported_Data,1);
            PeakListMW = Sorted(:,1);
            Intensity = Sorted(:,2);
            setappdata(0,'Peaklist_Data',Imported_Data);
            setappdata(0,'Fragments_Masses',PeakListMW);
            setappdata(0,'Int',Intensity);
            Experimental_Protein_Mass = Imported_Data(1,1);
        end
        Peaks_File = Imported_Data;
        setappdata(0,'Peaklist_Data',Imported_Data);
        
    elseif(getappdata(0,'Type_file')==2)
        file_data = mzxmlread(Batch_Search_File);
        data_all=mzxml2peaks(file_data, 'Levels', 2);
        [dummy, Index] = sort(cellfun('size',data_all,1), 'descend'); %#ok<*ASGLU>
        data_all_sort=data_all(Index);
        Imported_Data =data_all_sort{1};
        m = max(Imported_Data(:,2));
        HandleIon = getappdata(0,'HandleIon');
        for tupleIndex = 1 : size(Imported_Data,1)
            %% Used for converting mono isotopic mass into mz value
            if (HandleIon{9,1} == 0)
                Imported_Data(tupleIndex,1) = Imported_Data(tupleIndex,1)+1.00727647;
            end
            Imported_Data(tupleIndex,2) = Imported_Data(tupleIndex,2)/m;
        end
        %added 26042021
        if FPAvalue == 1 %Is FPA checked?
            Comp_Peaks_File = FinalImp_Data;
            Sorted = sortrows(FinalImp_Data,1);
            PeakListMW = Sorted(:,1);
            Intensity = Sorted(:,2);
            setappdata(0,'Comp_Peaklist_Data',FinalImp_Data);
            setappdata(0,'Comp_Fragments_Masses',PeakListMW_Comp);  % MS2s (in ascending order) then, MS1 only 
            setappdata(0,'Comp_Int',Intensity_Comp);   
            if (isempty(FinalImp_Data))
                continue
            end
            Experimental_Protein_Mass = FinalImp_Data(1,1);
        else
            Peaks_File = Imported_Data;
            Sorted = sortrows(Imported_Data,1);
            PeakListMW = Sorted(:,1);
            Intensity = Sorted(:,2);
            setappdata(0,'Peaklist_Data',Imported_Data);
            setappdata(0,'Fragments_Masses',PeakListMW);
            setappdata(0,'Int',Intensity);
            [data_all,times]=mzxml2peaks(file_data, 'Levels', 1);
            [dummy, Index] = sort(cellfun('size',data_all,1), 'descend');
            Protein_Mass=data_all{Index(1)};
            Experimental_Protein_Mass=Protein_Mass(length(Protein_Mass));
        end
    else
        Imported_Data =  importdata(Batch_Search_File);
        m = max(Imported_Data(:,2));
        HandleIon = getappdata(0,'HandleIon');
        for tupleIndex = 2 : size(Imported_Data,1)
            %% Used for converting mono isotopic mass into mz value
            if (HandleIon{9,1} == 0)
                Imported_Data(tupleIndex,1) = Imported_Data(tupleIndex,1)+1.00727647;
            end
            Imported_Data(tupleIndex,2) = Imported_Data(tupleIndex,2)/m;
        end
    %Is FPA checked?
    FPAvalue = getappdata(0,'FPAvalue');
    if isempty(FPAvalue)
        FPAvalue =0;
    end
        %added 26042021
        if FPAvalue == 1 %Is FPA checked?
            Comp_Peaks_File = FinalImp_Data;
            Sorted = sortrows(FinalImp_Data,1);
            PeakListMW = Sorted(:,1);
            Intensity = Sorted(:,2);
            setappdata(0,'Comp_Peaklist_Data',FinalImp_Data);
            setappdata(0,'Comp_Fragments_Masses',PeakListMW_Comp);  % MS2s (in ascending order) then, MS1 only 
            setappdata(0,'Comp_Int',Intensity_Comp);   
            if (isempty(FinalImp_Data))
                continue
            end
            Experimental_Protein_Mass = FinalImp_Data(1,1);

        else
            Peaks_File = Imported_Data;
            Sorted = sortrows(Imported_Data,1);
            PeakListMW = Sorted(:,1);
            Intensity = Sorted(:,2);
            setappdata(0,'Peaklist_Data',Imported_Data);
            setappdata(0,'Fragments_Masses',PeakListMW);
            setappdata(0,'Int',Intensity);
            Experimental_Protein_Mass = Imported_Data(1,1);
        end
    end

     %%% Updated 20210409   %% BELOW
    % Adding here to cater all file formats
if FPAvalue == 1 %Is FPA checked?
    temp_Imported_Data = FinalImp_Data(2:end, :);
    Sorted = sortrows(temp_Imported_Data,1);
    %Sorted = sortrows(Imported_Data,1);
    PeakListMW_Comp = Sorted(:,1);
    Intensity_Comp = Sorted(:,2);
    PeakListMW_Comp(end+1) = FinalImp_Data(1,1);
    Intensity_Comp(end+1) = FinalImp_Data(1,2);
    setappdata(0,'Comp_Peaklist_Data',FinalImp_Data);  % MS1 then, MS2s with their intensities [Original Order]
    setappdata(0,'Comp_Fragments_Masses',PeakListMW_Comp);  % MS2s (in ascending order) then, MS1 only 
    setappdata(0,'Comp_Int',Intensity_Comp);                % Only intensities of corresponding to the PeakListMW (MS values)
    %%% Updated 20210409   %% ABOVE
else
     %%% Updated 20210409   %% BELOW
    % Adding here to cater all file formats
    temp_Imported_Data = Imported_Data(2:end, :);
    Sorted = sortrows(temp_Imported_Data,1);
    %Sorted = sortrows(Imported_Data,1);
    PeakListMW = Sorted(:,1);
    Intensity = Sorted(:,2);
    PeakListMW(end+1) = Imported_Data(1,1);
    Intensity(end+1) = Imported_Data(1,2);
    setappdata(0,'Peaklist_Data',Imported_Data);  % MS1 then, MS2s with their intensities [Original Order]
    setappdata(0,'Fragments_Masses',PeakListMW);  % MS2s (in ascending order) then, MS1 only 
    setappdata(0,'Int',Intensity);                % Only intensities of corresponding to the PeakListMW (MS values)
    %%% Updated 20210409   %% ABOVE
end
    
    try
        %% tuner
        if getappdata(0,'tuner') == 1 % Executes when Auto-Tune Checkobx is on
            Slider_Value=1;
            [Tuned_MolWt, ~,~ ,~,~] =  MassTuner(Slider_Value,Experimental_Protein_Mass,MW_Tolerance);
            if(abs(Tuned_MolWt - Experimental_Protein_Mass) < 3)
                Tuned_Mass = Tuned_MolWt;
                Experimental_Protein_Mass = Tuned_Mass;
               %Upated 20210410 BELOW
%                 set(handles.edit_TunedMass,'String',Tuned_Mass); % The fetched Tuned Mass assigned to the Tuned Mass in Main window
%                 set(handles.edit_TunedMass,'ForegroundColor',[0 0 0]); % Black Color
                %Upated 20210410 ABOVE
            end % if ((get(handles.checkbox_TunedMass, 'Value')) == 1)
        end
        
        %% Engine
        if FilterPSTs == 1
            Tags_Ladder = Prep_Score_PSTs(str2num(User_EST_parameters{1}),str2num(User_EST_parameters{2}),str2num(User_EST_parameters{3}),str2num(User_EST_parameters{4})); %#ok<ST2NM>
        else
            Tags_Ladder = {};
        end
        
        FilterDB=getappdata(0,'FilterDB');
        if FilterDB
            [Candidate_ProteinsList_Initial,Candidate_ProteinsList_Truncated] = ParseDatabaseBatch(Candidate_ProteinsList_Full,Experimental_Protein_Mass,MW_Tolerance,Fixed_Modifications,Variable_Modifications,Tags_Ladder);
        else
            Candidate_ProteinsList_Initial = Candidate_ProteinsList_Full;
            %% Obtain proteins greater then MolW.
            Candidate_ProteinsList_Truncated = Candidate_ProteinsList_Full;
        end
        
        Candidate_ProteinsList_MW=Score_Mol_Weight(Candidate_ProteinsList_Initial,Experimental_Protein_Mass);
        Candidate_ProteinsList = Updated_ParseDatabase(Experimental_Protein_Mass,Tags_Ladder,Candidate_ProteinsList_MW,PTM_Tolerance,Fixed_Modifications,Variable_Modifications );
        Candidate_ProteinsList = Insilico_frags_Generator_modifier(Candidate_ProteinsList,Fragmentation_Type,HandleIon);
        if BlindPTMvar == 1
            [sizeHopeInfo,Hop_Info_name, Hop_Info_AA, Hop_Info_end, Hop_Info_start] = BlindPTM_Extraction();
            Candidate_ProteinsListModified = BlindPTM(Candidate_ProteinsList,Peptide_Tolerance, 1,PepUnit,sizeHopeInfo,Hop_Info_name, Hop_Info_AA, Hop_Info_end, Hop_Info_start);
        else
            Candidate_ProteinsListModified = [];
        end
        Candidate_ProteinsList = [Candidate_ProteinsList; Candidate_ProteinsListModified];
        %Is FPA checked?
        FPAvalue = getappdata(0,'FPAvalue');
        if isempty(FPAvalue)
            FPAvalue =0;
        end
      
        %remove if else condition later because statement sre same for FPA
        %and nFPA
        if FPAvalue ==1
%         Matches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Comp_Peaklist_Data'),Peptide_Tolerance,PepUnit); %% Use mass diff to estimate blind ptms.
        Matches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Comp_Peaklist_Data'), getappdata(0,'Comp_Int'),Peptide_Tolerance,PepUnit); %Updated 20210410 %% Use mass diff to estimate blind ptms.
        else
        Matches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Fragments_Masses'), getappdata(0,'Int'),Peptide_Tolerance,PepUnit); %Updated 20210410 %% Use mass diff to estimate blind ptms.
%         Matches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Peaklist_Data'),Peptide_Tolerance,PepUnit); %% Use mass diff to estimate blind ptms.
        end
            
        Matches = BlindPTM_Localization(Matches,Experimental_Protein_Mass);
        
        %% Truncation Starts Here
        truncation = getappdata(0,'Truncation');
        if isempty(truncation)
            truncation = 0;
        end
        TruncatedMatches = [];
        if truncation == 1
            %% Candidate_ProteinsList_Truncated = ParseDatabaseBatch_Truncation(Candidate_ProteinsList_Full,Experimental_Protein_Mass,MW_Tolerance,Tags_Ladder);
            [Candidate_ProteinsList_Left,Candidate_ProteinsList_Right] = PreTruncation(Candidate_ProteinsList_Truncated,MW_Tolerance);
            
            %% Ordinary Truncation
            [Candidate_ProteinsList_Left, RemainingProteins_Left] = Truncation_Left(Candidate_ProteinsList_Left);
            [Candidate_ProteinsList_Right, RemainingProteins_Right] = Truncation_Right(Candidate_ProteinsList_Right);
            Candidate_ProteinsList_UnModified = [Candidate_ProteinsList_Left; Candidate_ProteinsList_Right];
            
            %% In silico Treatment
            Candidate_ProteinsList_UnModified = Insilico_frags_Generator_modifier(Candidate_ProteinsList_UnModified,Fragmentation_Type,HandleIon);
            if BlindPTMvar == 1
                RemainingProteins_Right = Insilico_frags_Generator_modifier(RemainingProteins_Right,Fragmentation_Type,HandleIon);
                RemainingProteins_Left = Insilico_frags_Generator_modifier(RemainingProteins_Left,Fragmentation_Type,HandleIon);
                
                %% Blind PTM localization
                [RemainingProteins_Right_Modified] = BlindPTM_Truncation_Right(RemainingProteins_Right, Peptide_Tolerance, 1, PepUnit, sizeHopeInfo,Hop_Info_name, Hop_Info_AA, Hop_Info_end, Hop_Info_start);
                [RemainingProteins_Left_Modified] = BlindPTM_Truncation_Left(RemainingProteins_Left, Peptide_Tolerance, 1, PepUnit, sizeHopeInfo,Hop_Info_name, Hop_Info_AA, Hop_Info_end, Hop_Info_start);
                
                %% Modified Truncation
                %20181029 - Adding additional parameters for pepUnits etc
                %[Candidate_ProteinsList_Modified_Right] = Truncation_Right_Modification(RemainingProteins_Right_Modified, Fragmentation_Type);
                [Candidate_ProteinsList_Modified_Right] = Truncation_Right_Modification(RemainingProteins_Right_Modified, Peptide_Tolerance, PepUnit, Fragmentation_Type);
                [Candidate_ProteinsList_Modified_Left] = Truncation_Left_Modification(RemainingProteins_Left_Modified, Peptide_Tolerance, PepUnit, Fragmentation_Type);
            else
                Candidate_ProteinsList_Modified_Right = [];
                Candidate_ProteinsList_Modified_Left = [];
            end
            %% Final Scoring
            Candidate_ProteinsList = [Candidate_ProteinsList_UnModified;Candidate_ProteinsList_Modified_Left;Candidate_ProteinsList_Modified_Right];
            Candidate_ProteinsList = FIlter_Truncated_Proteins(Candidate_ProteinsList, Tags_Ladder);
            
            
        %Is FPA checked?
        FPAvalue = getappdata(0,'FPAvalue');
        if isempty(FPAvalue)
            FPAvalue =0;
        end
      
        %remove if else condition later because same
        if FPAvalue ==1
%           TruncatedMatches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Comp_Peaklist_Data'),Peptide_Tolerance,PepUnit); %% Use mass diff to estimate blind ptms.
            TruncatedMatches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Comp_Peaklist_Data'), getappdata(0,'Comp_Int'),Peptide_Tolerance,PepUnit); %Updated 20210410 %% Use mass diff to estimate blind ptms.

        else
            TruncatedMatches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Fragments_Masses'), getappdata(0,'Int'),Peptide_Tolerance,PepUnit); %Updated 20210410 %% Use mass diff to estimate blind ptms.
%             TruncatedMatches = Insilico_Score(Candidate_ProteinsList,getappdata(0,'Peaklist_Data'),Peptide_Tolerance,PepUnit); %% Use mass diff to estimate blind ptms.
        end
        end
        Matches = [Matches; TruncatedMatches]; %#ok<AGROW>
        setappdata(0,'Matches',Matches);
        tElapsed = num2str(toc(batchModeTimer));
        
        %% Compute Final Score
        final_score;
        
        %% Compute E-Value
        Matches=getappdata(0,'Matches');
        if numel(Matches) ~= 0
            Matches =  PrSM_Evalue(Matches);
        end
        
        %% For table
        data = cell(numel(Matches),32);
        %sort result with respect to their scores
        for z=1:numel(Matches)
            % get values from structure in which result are stored
            data{z,1} = Matches{z}.Header;
            data{z,2} = Matches{z}.Name;
            data{z,3} = Matches{z}.MolW;
            data{z,4} = Matches{z}.Final_Score;%
            data{z,5} = Matches{z}.Sequence;
            data{z,6} = Matches{z}.PTM_score;
            data{z,7} = Matches{z}.PTM_name;
            data{z,8} = Matches{z}.EST_Score;
            data{z,9} = Matches{z}.PTM_seq_idx;
            data{z,10} = Matches{z}. PTM_site;
            data{z,11} = Matches{z}. MWScore;
            data{z,12} = Matches{z}.Matches_Score;
            data{z,13} = Matches{z}.LeftIons;
            data{z,14} = Matches{z}.RightIons;
            
            if BlindPTMvar == 1
                Blind = Matches{z}.Blind;
            else
                Blind.Start = -1;
                Blind.End = -1;
                Blind.Mass = -1;
            end
            
            data{z,15} = Blind.Start;
            data{z,16} = Blind.End;
            data{z,17} = Blind.Mass;
            data{z,18} = Matches{z}.Match;
            data{z,19} = Matches{z}.Terminal_Modification;
            
            if truncation == 1
                if strcmp(Matches{z}.Truncation,'Left')
                    TruncationMessage = 'Truncation at N-Terminal Side';
                elseif strcmp(Matches{z}.Truncation,'Right')
                    TruncationMessage = 'Truncation at C-Terminal Side';
                else
                    TruncationMessage = 'No Truncation';
                end
                TruncationLocation = Matches{z}.TruncationIndex;
            else
                TruncationMessage = 'No Truncation';
                TruncationLocation = -1;
            end
            data{z,20} = TruncationMessage;
            data{z,21} = TruncationLocation;
            data{z,22} = Matches{z}.Evalue;
        end
        
        %% Write File
        if numel(data) ~= 0
            data=sortrows(data,-4);
            menu_File = fopen(Save_Batch_File,'a');
            size_data=size(data);
            for out_put=1:size_data(1)
                
                if numel(Tags_Ladder)~= 0 % PST determination
                    pst = [];
                    a = [];
                    PSTend= [];
                    S2 = [];
                    S1 = [];
                    for idx=1:numel(Tags_Ladder)
                        PST = strfind(data(out_put,5), Tags_Ladder{1,idx}{1,1});
                        if (isempty(PST{:}))
                            continue %#ok<*AGROW>
                        else
                            c = cell2mat(PST); %to check if number of elements at start position are more than 1
                            if numel(c) == 1
                                for v=1:numel(c)
                                    S1 = cellfun(@(x) x(1,1),PST);
                                    a(v).start = S1;
                                    a(v).tag = Tags_Ladder{1,idx}{1,1};
                                    a(v).length = numel(a(v).tag);
                                    y = a.length;
                                    PSTend = S1+ y;
                                    a(v).end = PSTend;
                                end
                            elseif numel(c) > 1 %more than one start value for same tag
                                for k=1:numel(c)
                                    S2 = PST{1,1}(1,k); %tags catered if start position is diffeent but tag is same
                                    a(k).start = S2;
                                    a(k).tag = Tags_Ladder{1,idx}{1,1};
                                    a(k).length = numel(a(k).tag);
                                    x = a.length;
                                    PSTend = S2+ x;
                                    a(k).end = PSTend;
                                end
                            end
                            pst = [pst; a(:)];%#ok<*AGROW>   %psts with their seqs, lengths, start and end positions
                            a= [];
                        end
                    end
                    
                    %% for iteration converted to cell
                    cellnumstart= {};
                    cellnumend = {};
                    indexstart = [];
                    indexend = [];
                    tagsequence = {};
                    tagseq = {};
                    taglen = {};
                    taglength = {};
                    tags ={};
                    startpst = {};
                    endpst = {};
                    
                    % finding unique tags with different start value
                    if(~isempty(pst))
                        for idx = 1:numel(pst)
                            cellnumstart =[cellnumstart; pst(idx).start];
                            cellnumend = [cellnumend; pst(idx).end];
                        end
                        indexstart = unique(cell2mat(cellnumstart));  %unique start position to avoid repetition
                        indexend = unique(cell2mat(cellnumend));       %unique end position to avoid repetition
                        for j=1:numel(indexstart)
                            for k=1:numel(cellnumstart)
                                if indexstart(j) == cell2mat(cellnumstart(k))
                                    tagsequence = pst(k).tag;     %finding sequences
                                    taglength = pst(k).length;   %finding lengths
                                    startofpst = pst(k).start;  %finding start position
                                    endofpst = pst(k).end;    %finding end position
                                end
                            end
                            tagseq = [tagseq; tagsequence]; %list of sequences
                            taglen = [taglen; taglength];   %list of lengths
                            startpst = [startpst; startofpst];   %list of start positions
                            endpst = [endpst; endofpst];         %list of end positions
                        end
                        
                        %% getting unique pst end values
                        tags = {'sequence','length', 'start', 'end'};
                        tags = [tags; tagseq, taglen, startpst, endpst];   %unique tags with unique start positions
                        sameend = {};
                        sameendtags = {};
                        seqs ={};
                        lens = {};
                        len = {};
                        endings ={};
                        ending ={};
                        
                        %% looping through to find tags with same end position to isolate later
                        for j = 1:size(tags)
                            for k=2:size(tags)
                                if j == k   %for the same tag do not compare but continue the loop
                                    continue
                                end
                                if isequal(tags{j,4}, tags{k,4}) %if the end position of tags is same then get them and make a list
                                    seqs = tags{j};
                                    len = tags{j,2};
                                    ending = tags{j,4};
                                end
                            end
                            sameend = [sameend; seqs];   %list of tags with same end position
                            lens = [lens; len]; %list of lengths of the tags with same end position
                            endings = [endings; ending];
                            seqs =[];
                            len =[];
                            ending =[];
                        end
                        sameend = {sameend, lens, endings};
                        sameendtags =[sameend{:}];   %list of tags with same end position
                        
                        %% if the list of same end tags is empty that is tags have different/unique end values
                        uniqueendtags = {};
                        uniqueend =[];
                        if (isempty(sameendtags)) %getting tags with different/unique end values
                            for x=2:size(tags)
                                seqs = tags{x,1};    %then get all the tags in 'tags' list which already has unique start values
                                lens = tags{x,2};
                                uniqueend = {seqs, lens};
                                uniqueendtags = [uniqueendtags; uniqueend];
                            end
                            uniqueendtags = [uniqueendtags]; %list of tags with same end position
                        end
                        
                        %% if there are no tags with same ends then loop through unique tags
                        if (isempty(sameendtags))
                            findpst = [];
                            findtag = [];
                            tagswithoutmax = {};
                            taglistwithoutmax = {};
                            maxtag = [];
                            lenmat ={};
                            MaxValue = [];
                            findmax =[];
                            for i=2:size(tags)
                                for j = 1:size(uniqueendtags)    %loop through list of tags with same end
                                    findpst =  strcmp (tags{i,1}, uniqueendtags{j,1});    %compare with the tags with unique start values
                                    if (isempty(findpst))   %if nothing is returned break the loop
                                        break
                                    end
                                    if findpst == 1    %if there is a match then
                                        taglen = tags{i,2};    %find length of matched tag
                                        lenmat = cell2mat(uniqueendtags(:,2));   %convert this to a matrix to loop through easily
                                        for j = 1:size(lenmat) %matrix of lengths of matched tags
                                            if lenmat(j) >= lenmat(:,1)    %find max length
                                                MaxValue = lenmat(j);  %find maximum value in the list by looping through it
                                            end
                                        end
                                        if taglen == MaxValue   %if the tag length equals maximum value & lenmat(:) ~= MaxValue
                                            for x=1:size(lenmat)
                                                findmax =  strfind(lenmat(x,1), MaxValue);
                                                if findmax == 1
                                                    maxtag = tags{i,1};    %get the tag with max value
                                                    break
                                                end
                                            end
                                            %else findx = strcmp(maxtag(1,1), tags(i,1))
                                            %if findx == 0
                                        else
                                            tagswithoutmax = [tags{i,1}];    %if it is not a match then add it to separate list
                                        end
                                    end
                                end
                                taglistwithoutmax = [taglistwithoutmax; tagswithoutmax];   %list of tags with unique ends
                                tagswithoutmax = []; %empty to remove previous value from loop
                            end
                        end
                        
                        %% if tag list without max is empty check all lenghths then get all tags for there is no max
                        findmaxtag =[];
                        if (isempty(sameendtags))
                            if (isempty(taglistwithoutmax))
                                for i=2:size(tags)
                                    findmaxtag =  strfind(maxtag(1,1), tags{i,1});
                                    if findmaxtag == 1
                                        break
                                    else
                                        tagswithoutmax = [tags{i,1}];
                                    end
                                    taglistwithoutmax = [taglistwithoutmax; tagswithoutmax];
                                end
                                taglistwithoutmax = [taglistwithoutmax]; %list of tags with unique ends
                            end
                        end
                        
                        %% if eliminating max tag does not leave behind anything in the list
                        taglistnomaxintags ={};
                        nomaxintags =[];
                        if (~isempty(uniqueendtags))
                            if (isempty(taglistwithoutmax))%if the list of same end tags is empty that is tags have different/unique end values
                                for i=2:size(tags)
                                    seqs = tags{i,1};    %then get all the tags in 'tags' list which already has unique start values
                                    lens = tags{i,2};
                                    nomaxintags = {seqs, lens};
                                    taglistnomaxintags = [taglistnomaxintags; nomaxintags];
                                end
                                taglistnomaxintags = [taglistnomaxintags]; %list of tags with same end position
                            end
                        end
                        
                        %% Find sub psts with same end to eliminate
                        if (isempty(sameendtags))
                            pststags={};
                            tagsnew = {};
                            tagfindnew =[];
                            for i=2:size(tags)
                                for j=1:size(uniqueendtags)
                                    tagfindnew =  strcmp(tags(i,1) , uniqueendtags{j,1});   %compare tag list with updated list of pst with same end to eliminate
                                    if tagfindnew == 1
                                        pststags = [];
                                        break
                                    else
                                        pststags = [tags(i,:)];
                                    end
                                end
                                tagsnew = [tagsnew; pststags];
                                pststags = [];
                            end
                            
                            if (isempty(tagsnew))  %if this list is empty that is tags minus max value is empty
                                for i=2:size(tags)
                                    pststags = [tags(i,:)];    %then get all the tags in 'tags' list which exist for it is all that is
                                    tagsnew = [tagsnew; pststags];   %list of tags with same end position
                                end
                                tagsnew = [tagsnew];
                            end
                        end
                        
                        %% Getting tag of max length to avoid counting same tag multiple times for different lengths
                        if (~isempty(sameendtags))
                            findpst = [];
                            findtag = [];
                            tagswithoutmax = {};
                            taglistwithoutmax = {};
                            maxtag = {};
                            lenmat ={};
                            MaxValue = [];
                            g ={};
                            findmax =[];
                            for i=2:size(tags)
                                for j = 1:size(sameendtags)    %loop through list of tags with same end
                                    findpst =  strcmp (tags{i,1}, sameendtags{j,1});    %compare with the tags with unique start values
                                    if (isempty(findpst))   %if nothing is returned break the loop
                                        break
                                    end
                                    if findpst == 1    %if there is a match then
                                        taglen = tags{i,2};    %find length of matched tag
                                        lenmat = cell2mat(sameendtags(:,2));   %convert this to a matrix to loop through easily
                                        for j = 1:size(lenmat) %matrix of lengths of matched tags
                                            if lenmat(j) >= lenmat(:,1)    %find max length
                                                MaxValue = lenmat(j);  %find maximum value in the list by looping through it
                                            end
                                        end
                                        if taglen == MaxValue   %if the tag length equals maximum value & lenmat(:) ~= MaxValue
                                            for x=1:size(lenmat)
                                                findmax =  strfind(lenmat(x,1), MaxValue);
                                                if findmax == 1
                                                    maxtag = {tags{i,1}, tags{i,4}};    %get the tag with max value
                                                    break
                                                end
                                            end
                                            %else findx = strcmp(maxtag(1,1), tags(i,1))
                                            %if findx == 0
                                        else
                                            tagswithoutmax = {tags{i,1}, tags{i,4}};    %if it is not a match then add it to separate list
                                        end
                                    end
                                end
                                taglistwithoutmax = [taglistwithoutmax; tagswithoutmax];   %list of tags with unique ends
                                tagswithoutmax = []; %empty to remove previous value from loop
                            end
                        end
                        
                        %% if eliminating max tag does not leave behind anything in the list
                        if (isempty(taglistwithoutmax))%if the list of same end tags is empty that is tags have different/unique end values
                            for i=2:size(tags)
                                seqs = tags{i,1};    %then get all the tags in 'tags' list which already has unique start values
                                lens = tags{i,2};
                                tagswithoutmax = {seqs, lens};
                                taglistwithoutmax = [taglistwithoutmax; tagswithoutmax];
                            end
                            taglistwithoutmax = [taglistwithoutmax]; %list of tags with same end position
                        end
                        
                        %% Update list of pst with same end
                        if (~isempty(sameendtags))
                            tagsameend = {};   %subset of psts with same end
                            tagswithoutmax = {};
                            for j = 1:size(taglistwithoutmax)
                                for i=1:size(sameendtags)
                                    findtag =  strcmp(taglistwithoutmax(j,1) , sameendtags{i,1});   %compare updated taglist to get tags with same end
                                    if findtag == 1     %if match then get the matched pst
                                        findnum = strfind((cell2mat(taglistwithoutmax(j,2))) , sameendtags{i,3});
                                        if findnum == 1 %if match then see the end position
                                            tagswithoutmax = {taglistwithoutmax{j,1}, taglistwithoutmax{j,2}};
                                        end
                                    else    %if no match continue the loop
                                        continue
                                    end
                                end
                                tagsameend = [tagsameend; tagswithoutmax];    %subset of pst list updated with same end
                                tagswithoutmax = [];
                            end
                            
                            %% Find sub psts with same end to eliminate
                            pststags={};
                            tagsnew = {};
                            tagfindnew =[];
                            for i=2:size(tags)
                                for j=1:size(tagsameend)
                                    tagfindnew =  strcmp(tags(i,1) , tagsameend{j,1});   %compare tag list with updated list of pst with same end to eliminate
                                    if tagfindnew == 1   %if tags matched then see the end position
                                        findnumz = strfind((cell2mat(tags(i,4))) , tagsameend{j,2});
                                        if findnumz == 1 %if end position is same then skip the tag
                                            pststags = [];
                                        else    %if end position is different but tag is same then keep it
                                            pststags = [tags(i,:)];
                                        end
                                    else
                                        pststags = [tags(i,:)];
                                    end
                                end
                                tagsnew = [tagsnew; pststags];  %final tags list
                                pststags = [];
                            end
                            
                            if (isempty(tagsnew))  %if this list is empty that is tags minus max value is empty
                                pststags = [tags(i,:)];    %then get all the tags in 'tags' list which exist for it is all that is
                                tagsnew = [tagsnew; pststags];   %list of tags with same end position
                            end
                        end
                        
                        
                        %% Count psts for each protein according to their length and append PST Count in Result Files
                        count1 = 0;
                        count2 = 0;
                        count3 = 0;
                        count4 = 0;
                        count5 = 0;
                        count6 = 0;
                        count7 = 0;
                        count8 = 0;
                        count9 = 0;
                        count10 = 0;
                        
                        %PST determination
                        if(~isempty(tagsnew))
                            for idx = 1:size(tagsnew)
                                if tagsnew{idx,2} == 1
                                    count1 = count1+1;
                                elseif tagsnew{idx,2} == 2
                                    count2 = count2+1;
                                elseif tagsnew{idx,2} == 3
                                    count3 = count3+1;
                                elseif tagsnew{idx,2} == 4
                                    count4 = count4+1;
                                elseif tagsnew{idx,2} == 5
                                    count5 = count5+1;
                                elseif tagsnew{idx,2} == 6
                                    count6 = count6+1;
                                elseif tagsnew{idx,2} == 7
                                    count7 = count7+1;
                                elseif tagsnew{idx,2} == 8
                                    count8 = count8+1;
                                elseif tagsnew{idx,2} == 9
                                    count9 = count9+1;
                                else tagsnew{idx,2} == 10
                                    count10 = count10+1;
                                end
                            end
                        end
                        
                        %pst count of individual protein
                        TagsCount = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10; count1, count2, count3, count4, count5, count6, count7, count8, count9, count10};
                        
                        %% For table
                        data{out_put,23} = TagsCount{2,1};
                        data{out_put,24} = TagsCount{2,2};
                        data{out_put,25} = TagsCount{2,3};
                        data{out_put,26} = TagsCount{2,4};
                        data{out_put,27} = TagsCount{2,5};
                        data{out_put,28} = TagsCount{2,6};
                        data{out_put,29} = TagsCount{2,7};
                        data{out_put,30} = TagsCount{2,8};
                        data{out_put,31} = TagsCount{2,9};
                        data{out_put,32} = TagsCount{2,10};
                        
                        % results file with pst count of individual protein
                        formatSpec = '> %s | Score: %f | Molweight: %f | # Matched Fragments: %d | Terminal Modification: %s | E-Value: %g | PST Length 1= %d | PST Length 2= %d | PST Length 3= %d | PST Length 4= %d | PST Length 5= %d | PST Length 6= %d | PST Length 7= %d | PST Length 8= %d | PST Length 9= %d | PST Length 10= %d '; %#ok<*AGROW>
                        fprintf(menu_File,formatSpec,data{out_put,1},data{out_put,4},data{out_put,3},data{out_put,18},data{out_put,19},data{out_put,22},TagsCount{2,1},TagsCount{2,2},TagsCount{2,3},TagsCount{2,4},TagsCount{2,5},TagsCount{2,6},TagsCount{2,7},TagsCount{2,8},TagsCount{2,9},TagsCount{2,10}); %#ok<AGROW>
                        fprintf(menu_File,'\n');
                        fprintf(menu_File,data{out_put,5});
                        fprintf(menu_File,'\n');
                        for indi = 1 : numel(data{out_put,7})
                            %                     Modification = strcat('Modification Name: ',data{out_put,7}{indi},' | Modification Site: ',data{out_put,10}{indi},' | Site Index: ', num2str(data{out_put,9}(indi)) );
                            Modification = 'Modification Name: %s | Modification Site: %d | Site Index: %d';
                            fprintf(menu_File,Modification, data{out_put,7}{indi}, data{out_put,10}{indi}, data{out_put,9}(indi));
                            fprintf(menu_File,'\n');
                        end
                        if BlindPTMvar == 1
                            if data{out_put,15} ~= -1
                                %                         Modification = strcat('Modification Name: Unknown | Modification Weight: ', num2str(data{out_put,17}),' | Modification lies between index :',num2str(data{out_put,15}),'-',num2str(data{out_put,16}));
                                Modification = 'Modification Name: Unknown | Modification Weight: %f | Modification lies between index : %d - %d';
                                fprintf(menu_File, Modification', data{out_put,17}, data{out_put,15}, data{out_put,16});
                                fprintf(menu_File,'\n');
                            end
                        end
                    else %results file when no psts found
                        formatSpec = '> %s | Score: %f | Molweight: %f | # Matched Feagments: %d | Terminal Modification: %s | E-Value: %g '; %#ok<*AGROW>
                        fprintf(menu_File,formatSpec,data{out_put,1},data{out_put,4},data{out_put,3},data{out_put,18},data{out_put,19},data{out_put,22}); %#ok<AGROW>
                        fprintf(menu_File,'\n');
                        fprintf(menu_File,data{out_put,5});
                        fprintf(menu_File,'\n');
                        for indi = 1 : numel(data{out_put,7})
                            %                     Modification = strcat('Modification Name: ',data{out_put,7}{indi},' | Modification Site: ',data{out_put,10}{indi},' | Site Index: ', num2str(data{out_put,9}(indi)) );
                            Modification = 'Modification Name: %s | Modification Site: %d | Site Index: %d';
                            fprintf(menu_File,Modification, data{out_put,7}{indi}, data{out_put,10}{indi}, data{out_put,9}(indi));
                            fprintf(menu_File,'\n');
                        end
                        if BlindPTMvar == 1
                            if data{out_put,15} ~= -1
                                %                         Modification = strcat('Modification Name: Unknown | Modification Weight: ', num2str(data{out_put,17}),' | Modification lies between index :',num2str(data{out_put,15}),'-',num2str(data{out_put,16}));
                                Modification = 'Modification Name: Unknown | Modification Weight: %f | Modification lies between index : %d - %d';
                                fprintf(menu_File, Modification', data{out_put,17}, data{out_put,15}, data{out_put,16});
                                fprintf(menu_File,'\n');
                            end
                        end
                    end
                else %results file when no psts allowed
                    formatSpec = '> %s | Score: %f | Molweight: %f | # Matched Feagments: %d | Terminal Modification: %s | E-Value: %g'; %#ok<*AGROW>
                    fprintf(menu_File,formatSpec,data{out_put,1},data{out_put,4},data{out_put,3},data{out_put,18},data{out_put,19},data{out_put,22}); %#ok<AGROW>
                    fprintf(menu_File,'\n');
                    fprintf(menu_File,data{out_put,5});
                    fprintf(menu_File,'\n');
                    for indi = 1 : numel(data{out_put,7})
                        %                     Modification = strcat('Modification Name: ',data{out_put,7}{indi},' | Modification Site: ',data{out_put,10}{indi},' | Site Index: ', num2str(data{out_put,9}(indi)) );
                        Modification = 'Modification Name: %s | Modification Site: %d | Site Index: %d';
                        fprintf(menu_File,Modification, data{out_put,7}{indi}, data{out_put,10}{indi}, data{out_put,9}(indi));
                        fprintf(menu_File,'\n');
                    end
                    if BlindPTMvar == 1
                        if data{out_put,15} ~= -1
                            %                         Modification = strcat('Modification Name: Unknown | Modification Weight: ', num2str(data{out_put,17}),' | Modification lies between index :',num2str(data{out_put,15}),'-',num2str(data{out_put,16}));
                            Modification = 'Modification Name: Unknown | Modification Weight: %f | Modification lies between index : %d - %d';
                            fprintf(menu_File, Modification', data{out_put,17}, data{out_put,15}, data{out_put,16});
                            fprintf(menu_File,'\n');
                        end
                    end
                end
            end %print file
                fclose(menu_File);
                headerRegex =  data{1,1};
                headerRegex(regexp(headerRegex,'[,]'))=[];
                headerRegex = regexprep(headerRegex,'[\n\r]+','');
                csvData = [csvData; {DirectoryContents(b).name},headerRegex,{data{1,19}},{data{1,5}},{data{1,20}},{data{1,21}},data{1,4},{data{1,3}},numel(data{1,7}),{data{1,18}},tElapsed,{data{1,22}},{data{1,23}},{data{1,24}},{data{1,25}},{data{1,26}},{data{1,27}},{data{1,28}},{data{1,29}},{data{1,30}},{data{1,31}},{data{1,32}}]; %#ok<AGROW>
                cd(initial_path);
        else
            cd(initial_path);
            cd(getappdata(0,'Result_Folder'));
            menu_File = fopen(Save_Batch_File,'a');
            fprintf(menu_File,'\n');
            fprintf(menu_File,'No Result Found Please search with another set of parameters');
            fclose( menu_File);
            cd(initial_path)
        end
     isResultAvailable = 1;
        
    catch Error %results file if error or no results found 
        menu_File = fopen(Save_Batch_File,'a');
        fprintf(menu_File,'\n');
        fprintf(menu_File,'No Result Found Please search with another set of parameters');
        fclose( menu_File);
        setappdata(0,'P_condotion',isResultAvailable); % setappdata(0,'P_condotion',0);  %Updated 20210426
        continue ;
    end
    
end 
%% write excel file
cd(getappdata(0,'Result_Folder'));
fid = fopen(strcat(Project_Title,'.csv'),'w');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\r\n',csvData{1,:});
for csvIndex = 2:numel(csvData)/22
    fprintf(fid,'%s, %s, %s, %s, %s, %d, %f, %f, %d, %d, %s, %g, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\r\n',csvData{csvIndex,:});
end
fclose(fid);
cd(initial_path)
progressbar(1/1);
setappdata(0,'P_condotion',isResultAvailable);