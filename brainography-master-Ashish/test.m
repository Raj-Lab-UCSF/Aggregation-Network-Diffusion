function varargout = test(varargin)
% TEST MATLAB code for test.fig
%      TEST, by itself, creates a new TEST or raises the existing
%      singleton*.
%
%      H = TEST returns the handle to a new TEST or the handle to
%      the existing singleton*.
%
%      TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST.M with the given input arguments.
%
%      TEST('Property','Value',...) creates a new TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test

% Last Modified by GUIDE v2.5 13-May-2013 23:00:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_OpeningFcn, ...
                   'gui_OutputFcn',  @test_OutputFcn, ...
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


% --- Executes just before test is made visible.
function test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test (see VARARGIN)

% Choose default command line output for test
handles.output = hObject;
handles.mainHandle = getappdata(0,'mainHandle');

% Update handles structure
guidata(hObject, handles);
setappdata(0,'tempHandle',hObject);

N = guidata(handles.mainHandle);
currentVol = N(end).currentVol;
M = N(currentVol);
tmp = get(handles.uitable1,'Data');

if isempty(M.regionvalues) && ~isempty(M.brain_at)
    U = unique(M.brain_at(M.brain_at~=0));
    T = cell(size(U,1), 5);
    T(:,1) = num2cell(U);
elseif ~isempty(M.regionvalues)
    T = M.regionvalues;
end
set(handles.uitable1,'Data',T);

if isempty(M.brain_colormaprange)
    minVal = min(cell2mat(T(:,2)));
    maxVal = max(cell2mat(T(:,2)));
else
    minVal = M.brain_colormaprange(1);
    maxVal = M.brain_colormaprange(2);
end

try
    set(handles.pushbutton12,'BackgroundColor', M.singleColor);
catch
    set(handles.pushbutton12,'BackgroundColor',M.singleColor/255);
end
if M.singleColorFlag
    set(handles.checkbox4,'Value',1);
    set(handles.popupmenu1,'Enable','off');
end

cmString = get(handles.popupmenu1,'String');
if ~isempty(M.brain_colormapidx)
    if (M.brain_colormapidx > 0) && (M.brain_colormapidx < size(cmString,1) )
        set(handles.popupmenu1,'Value',M.brain_colormapidx);
        eval(['rawMap = ' M.brain_colormap '(200);']);
        setappdata(handles.output,'rawMap',rawMap);
    else
        rawMap = M.custom_colormap;
        set(handles.popupmenu1,'Value',M.brain_colormapidx);
        setappdata(handles.output,'rawMap',rawMap);
    end
end
set(handles.slider1,'Min',-500000,'Max',500000,'Value',minVal,'SliderStep',[1/1000000 0.1]);
set(handles.slider2,'Min',-500000,'Max',500000,'Value',maxVal,'SliderStep',[1/1000000 0.1]);
set(handles.edit1,'String',num2str(minVal));
set(handles.edit2,'String',num2str(maxVal));
    
set(handles.popupmenu1,'Value',M.brain_colormapidx);
% UIWAIT makes test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    regionvar = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    matIn = load([pathname filesep filename]);
    eval(['RV = matIn.' regionvar.name ';']);
    if length(RV) ~= size(T,1)
        disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(RV))]);
    else
        if size(RV,2) == size(T,1)
            RV = RV';
        end
        T(:,2) = num2cell(RV);
        set(handles.uitable1,'Data',T);
        set(handles.edit1,'String',num2str(min(RV)));
        set(handles.slider1,'Value',min(RV));
        set(handles.edit2,'String',num2str(max(RV)));
        set(handles.slider2,'Value',max(RV));
        if ~get(handles.checkbox4,'Value')
            setRVRGB(handles,getappdata(handles.output,'rawMap'));
        end
    end
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');

importRV = inputdlg({'Please enter the variable to import:'});
if ~isempty(importRV)
    try
        RV = evalin('base', importRV{1});
    catch
        RV=[];
        disp(['Variable name ''' importRV{1} ''' not found.']);
    end
    if ~isempty(RV)
        if length(RV) ~= size(T,1)
            disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(RV))]);
        else
            if size(RV,2) == size(T,1)
                RV = RV';
            end
            T(:,2) = num2cell(RV);
            set(handles.uitable1,'Data',T);
            set(handles.edit1,'String',num2str(min(RV)));
            set(handles.slider1,'Value',min(RV));
            set(handles.edit2,'String',num2str(max(RV)));
            set(handles.slider2,'Value',max(RV));
            if ~get(handles.checkbox4,'Value')
                setRVRGB(handles,getappdata(handles.output,'rawMap'));
            end
        end
    end
end
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    regionvar = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    if size(regionvar,1) < 2
        matIn = load([pathname filesep filename]);
        eval(['RV = matIn.' regionvar.name ';']);
        if length(RV) ~= size(T,1)
            disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(RV))]);
        else
            if size(RV,2) == size(T,1)
                RV = RV';
            end
            T(:,3:5) = num2cell(RV);
            
            set(handles.uitable1,'Data',T);
        end
    else
        disp('Please specify file containing single variable.');
    end
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');

importRV = inputdlg({'Please enter the variable to import:'});
if ~isempty(importRV)
    try
        RV = evalin('base', importRV{1});
    catch
        RV=[];
        disp(['Variable name ''' importRV{1} ''' not found.']);
    end
    if ~isempty(RV)
        if length(RV) ~= size(T,1)
            disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(RV))]);
        else
            if size(RV,2) == size(T,1)
                RV = RV';
            end
            T(:,3:5) = num2cell(RV);
            set(handles.uitable1,'Data',T);
        end
    end
end

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N = guidata(handles.mainHandle);
currentVol = N(end).currentVol;
M = N(currentVol);
M.regionvalues = get(handles.uitable1,'Data');
M.singleColorFlag = get(handles.checkbox4,'Value');
disp(['M.singleColorFlag ' num2str(M.singleColorFlag)]);
M.singleColor = get(handles.pushbutton12,'BackgroundColor');
M.brain_colormapidx = get(handles.popupmenu1, 'Value');
M.brain_colormaprange(1) = str2num(get(handles.edit1,'String'));
M.brain_colormaprange(2) = str2num(get(handles.edit2,'String'));

N(currentVol) = M;
guidata(handles.mainHandle,N);
close(handles.output);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');
T(:,3:5) = cell(size(T,1),3);
set(handles.uitable1,'Data',T);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');
T(:,2) = num2cell([1:size(T,1)]');
set(handles.uitable1,'Data',T);
set(handles.edit1,'String','1');
set(handles.edit2,'String',num2str(size(T,1)));
set(handles.slider1,'Value',1);
set(handles.slider2,'Value',size(T,1));
if ~get(handles.checkbox4,'Value')
    setRVRGB(handles,getappdata(handles.output,'rawMap'));
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.output);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
cmVal=get(hObject,'Value');
cmString=get(hObject,'String');
customOpt = size(cmString,1);
if cmVal == customOpt
    waitfor(import_colormap);
    H = guidata(handles.mainHandle);
    rawMap = H(H(end).currentVol).custom_colormap;
    
else
    eval(['rawMap = ' cmString{cmVal} '(200);']);
    
end
setappdata(handles.output,'rawMap',rawMap);
setRVRGB(handles,rawMap);
guidata(hObject,handles);

function setRVRGB(handles,rawMap)
minVal = get(handles.slider1,'Value');
maxVal = get(handles.slider2,'Value');
rvTable = get(handles.uitable1,'Data');
entriez = cell2mat(rvTable(:,2));
rgbMapping = getRGBTriple(rawMap,minVal,maxVal,entriez);
rvTable(:,3:5) = num2cell(rgbMapping);
set(handles.uitable1,'Data',rvTable);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
N = guidata(handles.mainHandle);
currentVol = N(end).currentVol;
M = N(currentVol);

if get(hObject,'Value')
%     brain_at = M.brain_at;
    
    set(handles.pushbutton12,'Enable','on');
    singleColor = get(handles.pushbutton12,'BackgroundColor');

%     uniqueRegions = unique(brain_at(find(brain_at~=0)));
%     regionvalues = cell(size(uniqueRegions,1),5);
%     regionvalues(:,1) = num2cell(uniqueRegions);
%     regionvalues(:,2) = num2cell(ones(size(uniqueRegions,1),1));
    T = get(handles.uitable1,'Data');
    T(:,3:5) = num2cell(repmat(singleColor,size(T,1),1));

    set(handles.uitable1,'Data',T);
    set(handles.popupmenu1,'Enable','off');
    set(handles.edit1,'Enable','off');
    set(handles.edit2,'Enable','off');
    set(handles.slider1,'Enable','off');
    set(handles.slider2,'Enable','off');
else
	set(handles.pushbutton12,'Enable','off');
    set(handles.popupmenu1,'Enable','on');
    set(handles.edit1,'Enable','on');
    set(handles.edit2,'Enable','on');
    set(handles.slider1,'Enable','on');
    set(handles.slider2,'Enable','on');
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = uisetcolor(get(hObject,'BackgroundColor'));
set(hObject,'BackgroundColor',V);

regionvalues = get(handles.uitable1,'Data');
numberROI = size(regionvalues,1);
regionvalues(:,3:5) = num2cell(repmat(V,numberROI,1));
% regionvalues(:,2) = num2cell(ones(1,numberROI));
set(handles.uitable1,'Data',regionvalues);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
minVal = str2double(get(hObject, 'String'));
maxVal = str2double(get(handles.edit2, 'String'));
if minVal > maxVal
    minVal = maxVal;
    set(hObject, 'String', num2str(minVal));
end
set(handles.slider1, 'Value', minVal);
setRVRGB(handles,getappdata(handles.output,'rawMap'));


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderMin = get(hObject,'Value');
sliderMax = get(hObject,'Value');

if sliderMin > sliderMax
    sliderMin = sliderMax;
    set(hObject,'Value',sliderMin);
end
set(handles.edit1,'String',num2str(sliderMin));
setRVRGB(handles,getappdata(handles.output,'rawMap'));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
minVal = str2double(get(handles.edit1, 'String'));
maxVal = str2double(get(hObject, 'String'));
if minVal > maxVal
    maxVal = minVal;
    set(hObject, 'String', num2str(maxVal));
end
set(handles.slider2, 'Value', maxVal);
setRVRGB(handles,getappdata(handles.output,'rawMap'));

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderMin = get(hObject,'Value');
sliderMax = get(hObject,'Value');

if sliderMin > sliderMax
    sliderMax = sliderMin;
    set(hObject, 'Value', sliderMax);
end

set(handles.edit2, 'String', num2str(sliderMax));

setRVRGB(handles, getappdata(handles.output, 'rawMap'));


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
N = guidata(handles.mainHandle);
M = N(N(end).currentVol);
if get(hObject,'Value') %LUT mode on
%     set(handles.pushbutton3,'Enable','Off');
%     set(handles.pushbutton4,'Enable','Off');
%     set(handles.pushbutton9,'Enable','Off');
    set(handles.checkbox4,'Enable','Off');
    set(handles.pushbutton12,'Enable','Off');
    set(handles.popupmenu1,'Enable','Off');
    set(handles.edit1,'Enable','Off');
    set(handles.slider1,'Enable','Off');
    set(handles.edit1,'Enable','Off');
    set(handles.slider2,'Enable','Off');
    
    set(handles.pushbutton6,'Enable','On');
    set(handles.pushbutton7,'Enable','On');
    set(handles.pushbutton10,'Enable','On');
else %LUT mode off
%     set(handles.pushbutton3,'Enable','On');
%     set(handles.pushbutton4,'Enable','On');
%     set(handles.pushbutton9,'Enable','On');
    set(handles.checkbox4,'Enable','On');
    set(handles.pushbutton12,'Enable','On');
    set(handles.popupmenu1,'Enable','On');
    set(handles.edit1,'Enable','On');
    set(handles.slider1,'Enable','On');
    set(handles.edit1,'Enable','On');
    set(handles.slider2,'Enable','On');
    
    set(handles.pushbutton6,'Enable','Off');
    set(handles.pushbutton7,'Enable','Off');
    set(handles.pushbutton10,'Enable','Off');
end
