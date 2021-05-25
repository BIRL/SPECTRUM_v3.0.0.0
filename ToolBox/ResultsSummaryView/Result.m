%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SPECTRUM: A MATLAB Toolbox for Top-down Proteomics     %
%                           Version 3.0.0.0                        %
%        Copyright (c) Biomedical Informatics Research Laboratory, %
%          Lahore University of Management Sciences Lahore (LUMS), %
%                           Pakistan.                              %
%                (http://biolabs.lums.edu.pk/BIRL)                 %
%                    (safee.ullah@gmail.com)                       %
%                 Last Modified on: 25-May-2021                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = Result(varargin)
% UNTITLED1 MATLAB code for untitled1.fig
%      UNTITLED1, by itself, creates a new UNTITLED1 or raises the existing
%      singleton*.
%
%      H = UNTITLED1 returns the handle to a new UNTITLED1 or the handle to
%      the existing singleton*.
%
%      UNTITLED1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED1.M with the given input arguments.
%
%      UNTITLED1('Property','Value',...) creates a new UNTITLED1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled1

% Last Modified by GUIDE v2.5 28-Mar-2019 11:10:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @untitled1_OpeningFcn, ...
    'gui_OutputFcn',  @untitled1_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled1 is made visible.
function untitled1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled1 (see VARARGIN)

% Choose default command line output for untitled1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled1 wait for user response (see UIRESUME)
% uiwait(handles.Figure_Result);


% --- Outputs from this function are returned to the command line.
function varargout = untitled1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
progressbar('Please wait while SPECTRUM compiles results...');
Matches=getappdata(0,'Matches');
initial_path=pwd;
row={};

%Attach a control to GUI for inputting user-provided protein count
% output= getappdata(0,'Proteins_Output');
output = 1000;

if(numel(Matches)< output)
    num_output=numel(Matches);
else
    num_output= output;
end
BlindPTMvar = getappdata(0,'BlindPTM');
truncation = getappdata(0,'Truncation');

variable_progress = numel(Matches)*3;

for z=1:numel(Matches)
    progressbar(z/variable_progress);
    
    % get values from structure in which result are stored
    Name=Matches{z}.Name;
    Header = Matches{z}.Header;
    %     ID=Matches{z}.Id;
    Score=Matches{z}.Final_Score;%
    MolW=Matches{z}.MolW;
    PTM_score=Matches{z}.PTM_score;
    PTM_name=Matches{z}.PTM_name;
    Sequence=Matches{z}.Sequence;
    EST_Score=Matches{z}.EST_Score;
    PTM_seq_idx=Matches{z}.PTM_seq_idx;
    PTM_site=Matches{z}. PTM_site;
    MWScore=Matches{z}. MWScore;
    Matches_Score=Matches{z}.Matches_Score;
    LeftIons=Matches{z}.LeftIons;
    RightIons=Matches{z}.RightIons;
    Blind = Matches{z}.Blind;
    Terminal=Matches{z}.Terminal_Modification;
    Evalue = Matches{z}.Evalue;
    No_of_Matches = Matches{z}.Match;
    % store in a cell array
    
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
        TruncationLocation = '-';
    end
    
    if(z==1)
        data={Header,Name,num2str(MolW),num2str(Score),Sequence,PTM_score,PTM_name,EST_Score,PTM_seq_idx,PTM_site,MWScore,...
            Matches_Score,LeftIons,RightIons,Blind.Start,Blind.End,Blind.Mass,No_of_Matches,Terminal,TruncationMessage,TruncationLocation,Evalue,z};
    else
        protein_data={Header,Name,num2str(MolW),num2str(Score),Sequence,PTM_score,PTM_name,EST_Score,PTM_seq_idx,PTM_site,MWScore,...
            Matches_Score,LeftIons,RightIons,Blind.Start,Blind.End,Blind.Mass,No_of_Matches,Terminal,TruncationMessage,TruncationLocation,Evalue,z};
        data = cat(1,data,protein_data);
    end
end

data=sortrows(data,-4);
setappdata(0,'data',data);

for z = 1:num_output
    progressbar(z+(variable_progress/3)/variable_progress);
    Name=data{z,1};
    ID=data{z,2};
    Score=data{z,4};
    MolW=data{z,3};
    Terminal=data{z,19};
    row=cat(2,row,num2str(z));
    set(handles.uitable2,'RowName',row);
    if(z==1)
        protein_data= {Name,ID,num2str(MolW),num2str(Score),Terminal,'...view...'};
        set(handles.uitable2,'data',protein_data);
    else
        data_1=get(handles.uitable2,'Data');
        protein_data= {Name,ID,num2str(MolW),num2str(Score),Terminal,'...view...'};
        newRowdata = cat(1,data_1,protein_data);
        set(handles.uitable2,'data',newRowdata);
    end
end

cd(getappdata(0,'Result_Folder'));
file = getappdata(0,'Single_Peaklist_Data');
temp = strsplit(file,'\');
Peaklist_Data = temp(numel(temp));
Project_Title=getappdata(0,'Project_Title');
Save_File = strcat(Project_Title,'_',Peaklist_Data,'.results');
Tags_Ladder=getappdata(0,'Tags_Ladder');
%Is FPA checked?
FPAvalue = getappdata(0,'FPAvalue');
if isempty(FPAvalue)
    FPAvalue =0;
end

if FPAvalue ==1
    Fragments_Peaklist_Data = getappdata(0,'Comp_Peaklist_Data');
else
    Fragments_Peaklist_Data = getappdata(0,'Peaklist_Data');
end


%% looping for every protein in file
if numel(data) ~= 0
    menu_File = fopen(Save_File{1,1},'a');
    size_data=size(data);
    for out_put=1:size(data)
        progressbar((out_put + (2*variable_progress/3))/variable_progress);
        
        % PST determination
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
            
            %% PST count of individual protein
            TagCount = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10; count1, count2, count3, count4, count5, count6, count7, count8, count9, count10};
            
            %% Results file with pst count of individual protein
            header=strcat('>',data{out_put,1},' | Score:',data{out_put,4},' | Molweight:',data{out_put,3},' | # Matched Fragments:',num2str(data{out_put,18}),' | Terminal Modification:',data{out_put,19},' | E-Value:',num2str(data{out_put,22}),' | PST Length 1= ', num2str(TagCount{2,1}),' | PST Length 2= ', num2str(TagCount{2,2}),' | PST Length 3= ', num2str(TagCount{2,3}),' | PST Length 4= ', num2str(TagCount{2,4}),' | PST Length 5= ', num2str(TagCount{2,5}),' | PST Length 6= ', num2str(TagCount{2,6}),' | PST Length 7= ', num2str(TagCount{2,7}),' | PST Length 8= ', num2str(TagCount{2,8}),' | PST Length 9= ', num2str(TagCount{2,9}), ' | PST Length 10= ', num2str(TagCount{2,10})); %#ok<*AGROW>
            fprintf(menu_File,header);
            fprintf(menu_File,'\n');
            fprintf(menu_File,data{out_put,5});
            fprintf(menu_File,'\n');
            for indi = 1 : numel(data{out_put,7})
                Modification = strcat('Modification Name: ',data{out_put,7}{indi},' | Modification Site: ',data{out_put,10}{indi},' | Site Index: ', num2str(data{out_put,9}(indi)) );
                fprintf(menu_File,Modification);
                fprintf(menu_File,'\n');
            end
            
            if BlindPTMvar == 1
                if data{out_put,15} ~= -1
                    Modification = strcat('Modification Name: Unknown | Modification Weight: ', num2str(data{out_put,17}),' | Modification lies between index :',num2str(data{out_put,15}),'-',num2str(data{out_put,16}));
                    fprintf(menu_File,Modification);
                    fprintf(menu_File,'\n');
                end
            end
        end
    end
    fclose( menu_File);
else
    cd(initial_path);
    cd(getappdata(0,'Result_Folder'));
    menu_File = fopen(Save_File,'a');
    fprintf(menu_File,'\n');
    fprintf(menu_File,'No Result Found Please search with another set of parameters');
    fclose( menu_File);
    cd(initial_path)
end
cd(initial_path)
varargout{1} = handles.output;

% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
a=eventdata.Indices;
if a(2)==6
    data=getappdata(0,'data');
    
    Matches=getappdata(0,'Matches');
    setappdata(0,'Res',a(1));
    % Matches=getappdata(0,'Matches');
    Rank=data(a(1),23);
    setappdata(0,'Rank',Rank{1,1});
    try
        addpath(strcat(pwd,'/ResultsDetailedView'))
        %        close(Detail_view);
        Detail_view;
    catch exception
        errordlg(getReport(exception,'basic','hyperlinks','off'));
    end
end


% --- Executes on button press in pushbutton_Finish.
function pushbutton_Finish_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.Figure_Result)


% --- Executes on scroll wheel click while the figure is in focus.
function Figure_Result_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to Figure_Result (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on uitable2 and none of its controls.
function uitable2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_back.
function pushbutton_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = filesep;
close(handles.Figure_Result);
addpath(strcat(pwd,f,'Score'));
Final_Scoring;
setappdata(0,'w1',getappdata(0,'w1'));
setappdata(0,'w2',getappdata(0,'w2'));
setappdata(0,'w3',getappdata(0,'w3'));
final_score;

result=getappdata(0,'Matches');
setappdata(0,'Matches',result);

addpath(strcat(pwd,f,'ResultsSummaryView'));
%  der=Tag_sizeone(Tags_Ladder);
%   setappdata(0,'Tags_Ladder',Tag_ladder);
%  setappdata(0,'Proteins_Output',Proteins_Output);
if(getappdata(0,'FScoringWindowClosed') == false)
    Result;
    
end


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run ('PSTcount.m');
