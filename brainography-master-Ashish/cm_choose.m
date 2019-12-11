function varargout = cm_choose(varargin)
% CM_CHOOSE MATLAB code for cm_choose.fig
%      CM_CHOOSE, by itself, creates a new CM_CHOOSE or raises the existing
%      singleton*.
%
%      H = CM_CHOOSE returns the handle to a new CM_CHOOSE or the handle to
%      the existing singleton*.
%
%      CM_CHOOSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CM_CHOOSE.M with the given input arguments.
%
%      CM_CHOOSE('Property','Value',...) creates a new CM_CHOOSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cm_choose_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cm_choose_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cm_choose

% Last Modified by GUIDE v2.5 08-May-2013 12:13:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cm_choose_OpeningFcn, ...
                   'gui_OutputFcn',  @cm_choose_OutputFcn, ...
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


% --- Executes just before cm_choose is made visible.
function cm_choose_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cm_choose (see VARARGIN)

% Choose default command line output for cm_choose
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

mainHandle = getappdata(0,'mainHandle');
setappdata(handles.output,'mainHandle',mainHandle);
tmp = guidata(mainHandle);

if ~isempty(tmp(tmp(end).currentVol).connectivityMatrix)
    CM = tmp(tmp(end).currentVol).connectivityMatrix;
    axes(handles.axes1);
    imagesc(CM);
    setappdata(handles.output,'CM',CM);
end
% UIWAIT makes cm_choose wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cm_choose_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% setappdata(handles.output,'handlestruct',handles.output);

% --- Executes on button press in pushbutton1. OK button
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainHandle = getappdata(handles.output,'mainHandle');
CM = getappdata(handles.output,'CM');
tmp = guidata(mainHandle);
tmp(tmp(end).currentVol).connectivityMatrix = CM;
guidata(mainHandle, tmp);
close(handles.output);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.output);

% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3. the browse button
function pushbutton3_Callback(hObject, eventdata, handles)  %BROWSE button
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    CMvar = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    matIn = load([pathname filesep filename]);
    eval(['CM = matIn.' CMvar.name ';']);
    axes(handles.axes1);
    imagesc(CM);
    setappdata(handles.output,'CM',CM);
end
    
% --- Executes on button press in pushbutton4. Import CM from workspace
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
importCM = inputdlg({'Please enter the variable to import:'});
if ~isempty(importCM)
    try
        CM = evalin('base', importCM{1});
    catch
        CM=[];
        disp(['Variable name ''' importCM{1} ''' not found.']);
    end
    if ~isempty(CM)
        axes(handles.axes1);
        imagesc(CM);
        setappdata(handles.output,'CM',CM);
    end
end
    
% --- Executes on button press in pushbutton5. The CLEAR button
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(getappdata(handles.output,'CM'))
    rmappdata(handles.output,'CM');
    cla(handles.axes1);
end
