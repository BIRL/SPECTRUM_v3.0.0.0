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
function varargout = PSTcount(varargin)
% PSTCOUNT MATLAB code for PSTcount.fig
%      PSTCOUNT, by itself, creates a new PSTCOUNT or raises the existing
%      singleton*.
%
%      H = PSTCOUNT returns the handle to a new PSTCOUNT or the handle to
%      the existing singleton*.
%
%      PSTCOUNT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSTCOUNT.M with the given input arguments.
%
%      PSTCOUNT('Property','Value',...) creates a new PSTCOUNT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PSTcount_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PSTcount_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PSTcount

% Last Modified by GUIDE v2.5 02-Apr-2019 11:32:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PSTcount_OpeningFcn, ...
                   'gui_OutputFcn',  @PSTcount_OutputFcn, ...
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


% --- Executes just before PSTcount is made visible.
function PSTcount_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PSTcount (see VARARGIN)

% Choose default command line output for PSTcount
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PSTcount wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PSTcount_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;
Tags_Ladder=getappdata(0,'Tags_Ladder');
Tag = Tags_Ladder(1,:);
TagSeqs = [];
TagLengths = [];
for i= 1:numel(Tag)
    TagSeq = Tag{i}(1);
    TagLength = Tag{i}(2);
    TagSeqs = [TagSeqs; TagSeq];
    TagLengths = [TagLengths; TagLength];
end
Length= cell2mat(TagLengths);

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
for i=1:size(TagLengths)
    if TagLengths{i,1} == 1
        count1 = count1+1;
    elseif TagLengths{i,1} == 2
        count2 = count2+1;
    elseif TagLengths{i,1} == 3
        count3 = count3+1;
    elseif TagLengths{i,1} == 4
        count4 = count4+1;
    elseif TagLengths{i,1} == 5
        count5 = count5+1;
    elseif TagLengths{i,1} == 6
        count6 = count6+1;
    elseif TagLengths{i,1} == 7
        count7 = count7+1;
    elseif TagLengths{i,1} == 8
        count8 = count8+1;
    elseif TagLengths{i,1} == 9
        count9 = count9+1;
    else TagLengths{i,1} == 10
        count10 = count10+1;
    end
end
TagCount = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10; count1, count2, count3, count4, count5, count6, count7, count8, count9, count10};
hold on
MatData = cell2mat(TagCount);

cla reset; 

FPAvalue = getappdata(0,'FPAvalue');
if isempty(FPAvalue)
    FPAvalue =0;
end
if FPAvalue == 0
    bar(MatData(1,:), MatData(2,:), 'facecolor', [0 0.45 0.74]); % plot PST length and count
elseif FPAvalue == 1
        bar(MatData(1,:), MatData(2,:), 'facecolor', [0 0.5 0]); % plot PST length and count
end
% bar(MatData(1,:), MatData(2,:)); % plot PST length and count
hold on

xlim([0 10]);
hold on
xlabel('PST Length');
yline= max(MatData(2,:)+10);
graph_axes = size(yline+2);
ylim([0 yline]);
hold on
if yline < 20
set(gca, 'YTick', 0: yline+2)
end
ylabel('PST Count');
set(gca, 'TickDir','out');
set(gcf,'menu', 'none');
set(gcf,'toolbar', 'figure');
datacursormode on;

varargout{1} = handles.output;
