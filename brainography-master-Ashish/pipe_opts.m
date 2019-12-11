function varargout = pipe_opts(varargin)
% PIPE_OPTS MATLAB code for pipe_opts.fig
%      PIPE_OPTS, by itself, creates a new PIPE_OPTS or raises the existing
%      singleton*.
%
%      H = PIPE_OPTS returns the handle to a new PIPE_OPTS or the handle to
%      the existing singleton*.
%
%      PIPE_OPTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIPE_OPTS.M with the given input arguments.
%
%      PIPE_OPTS('Property','Value',...) creates a new PIPE_OPTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pipe_opts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pipe_opts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pipe_opts

% Last Modified by GUIDE v2.5 13-May-2013 17:25:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pipe_opts_OpeningFcn, ...
                   'gui_OutputFcn',  @pipe_opts_OutputFcn, ...
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


% --- Executes just before pipe_opts is made visible.
function pipe_opts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pipe_opts (see VARARGIN)

% Choose default command line output for pipe_opts
handles.output = hObject;
handles.mainHandle = getappdata(0, 'mainHandle');
% Update handles structure
guidata(hObject, handles);

N = guidata(handles.mainHandle);
M = N(N(end).currentVol);

if ~isempty(M.pipeScale)
    set(handles.slider1,'Value', M.pipeScale);
    setappdata(handles.output,'sliderValue', M.pipeScale);
else
    set(handles.slider1,'Value',1.5);
    setappdata(handles.output,'sliderValue', 1.5);
end

if(M.pipeUniform)
    set(handles.checkbox1,'Value',1);
else
    set(handles.checkbox1,'Value',0);
end

set(handles.popupmenu4, 'Value', M.pipeScheme);
setPipeScheme(handles, M.pipeScheme);

setappdata(handles.output,'cMapString',M.pipeColorMap);

if ~isempty(M.connectivityMatrix)
    setappdata(handles.output,'CM',M.connectivityMatrix);
    cmMin = min(min(M.connectivityMatrix));
    cmMax = max(max(M.connectivityMatrix));
    set(handles.slider3,'Min',cmMin,'Max',cmMax);
    cmDim = size(M.connectivityMatrix);
    if isempty(M.pipeColorHyperCube)
        pipeColorHyperCube = zeros([cmDim 3 9]);
    else
        pipeColorHyperCube = M.pipeColorHyperCube;
    end
    pipeCoupletThreshold = M.pipeCoupletThreshold;
    set(handles.slider3,'Value',pipeCoupletThreshold);
    set(handles.edit2,'String',num2str(pipeCoupletThreshold));
    setappdata(handles.output,'pipeCoupletThreshold',pipeCoupletThreshold);
    set(handles.slider1,'Value',M.pipeScale,'Min',0.01,'Max',20.0);
    set(handles.edit1,'String',num2str(M.pipeScale));
    set(handles.popupmenu3,'Value',M.pipeStyle);
    axes(handles.axes2);
    switch M.pipeScheme
        case 1
            singletCube = pipeColorHyperCube(:,:,:,1);
            if ~any(any(sum(singletCube,3))) || length(unique(sum(singletCube,3))) > 1
                singlet = rand(1,3);
            else
                singlet = squeeze(singletCube(1,1,:))';
            end
            pipeColorHyperCube = setSingletCube(handles,singlet,pipeColorHyperCube);
            setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
        case 2
            eval(['cMap = ' M.pipeColorMap '(200);']);
            set(gcf,'colormap', cMap);
            setappdata(handles.output,'cMapString',M.pipeColorMap);
            pipeColorHyperCube = setCmapCube(handles,cMap,pipeColorHyperCube);
            setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
            cmString = get(handles.popupmenu1,'String');
            set(handles.popupmenu1,'Value',find(ismember(cmString,M.pipeColorMap)));
        case 3
            if M.pipeCoupletThreshold > cmMax || M.pipeCoupletThreshold < cmMin
                pipeCoupletThreshold = cmMin + (cmMax - cmMin)/2;
                setappdata(handles.output,'pipeCoupletThreshold',pipeCoupletThreshold);
            else
                pipeCoupletThreshold = M.pipeCoupletThreshold;
            end
            couplet = M.pipeCouplet;
            pipeColorHyperCube = setCoupletCube(handles,couplet,pipeCoupletThreshold,pipeColorHyperCube);
    end
    setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
    setappdata(handles.output,'pipeCoupletThreshold',pipeCoupletThreshold);
    setappdata(handles.output,'pipeCouplet',M.pipeCouplet);
    setappdata(handles.output,'cMapString',M.pipeColorMap);
    imagesc(M.connectivityMatrix);
else
    setappdata(handles.output,'CM',[]);
end

function setHueBox(RGB,boxHandle)
set(boxHandle,'BackgroundColor',RGB);
set(boxHandle,'ForegroundColor',RGB);
% if ~isempty(M.pipeColorCube)
%     
    

function setPipeScheme(handles, pipeScheme)
% Set GUI for choosing proper options
switch pipeScheme
    case 1 %single color
        set(handles.popupmenu1,'Enable','off');
        set(handles.pushbutton5,'Enable','on');
        set(handles.pushbutton6,'Enable','off');
        set(handles.slider3,'Enable','off');
        set(handles.edit2,'Enable','off');
    case 2 %color map
        set(handles.popupmenu1,'Enable','on');
        set(handles.pushbutton5,'Enable','off');
        set(handles.pushbutton6,'Enable','off');
        set(handles.slider3,'Enable','off');
        set(handles.edit2,'Enable','off');
    case 3 %two-color
        set(handles.popupmenu1,'Enable','off');
        set(handles.pushbutton5,'Enable','on');
        set(handles.pushbutton6,'Enable','on');
        set(handles.slider3,'Enable','on');
        set(handles.edit2,'Enable','on');
end

function pch = setSingletCube(handles,singlet,pch)
CM = getappdata(handles.output,'CM');
if ~isempty(CM)
    cmDim = size(CM);
    singletCube = cat(3,singlet(1)*ones(cmDim),singlet(2)*ones(cmDim),singlet(3)*ones(cmDim));
    pch(:,:,:,1) = singletCube;
    set(gcf,'ColorMap',repmat(singlet,[5 1]));
    setHueBox(singlet,handles.pushbutton5);
%     setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
else
    disp('Please set connectivity');
end

function pch = setCmapCube(handles,cMap,pch)
CM = getappdata(handles.output,'CM');
if ~isempty(CM)
   rgbentriez = getRGBTriple(cMap,min(min(CM)),max(max(CM)),CM(:));
   rgbCM = cat(3,reshape(rgbentriez(:,1),size(CM)),reshape(rgbentriez(:,2),size(CM)),reshape(rgbentriez(:,3),size(CM)));
   pch(:,:,:,2) = rgbCM;
else
    disp('Please set connectivity.');
end


function pipeColorHyperCube = setCoupletCube(handles,couplet,pipeCoupletThreshold,pipeColorHyperCube)
CM = getappdata(handles.output,'CM');
if ~isempty(CM)
    cmMin = min(min(CM)); cmMax = max(max(CM));
    bigun = find(CM > pipeCoupletThreshold);
    lilun = find(CM < pipeCoupletThreshold);
    coupletCube = zeros([size(CM) 3]);
    for i=1:3
        tmp = coupletCube(:,:,i);
        tmp(bigun) = couplet(2,i);
        tmp(lilun) = couplet(1,i);
        coupletCube(:,:,i) = tmp;
    end
    pipeColorHyperCube(:,:,:,3) = coupletCube;
    threshfracidx = round((pipeCoupletThreshold - cmMin)/(cmMax - cmMin) * 150);
    set(gcf, 'ColorMap', [repmat(couplet(1,:),[threshfracidx 1]); repmat(couplet(2,:),[150-threshfracidx 1])]);
    setHueBox(couplet(1,:),handles.pushbutton5);
    setHueBox(couplet(2,:),handles.pushbutton6);
else
    disp('Please set connectivity');
end


% UIWAIT makes pipe_opts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pipe_opts_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.output);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N = guidata(handles.mainHandle);
M = N(N(end).currentVol);

M.pipeCouplet = getappdata(handles.output,'pipeCouplet');
M.pipeScale = get(handles.slider1,'Value');
M.pipeStyle = get(handles.popupmenu3,'Value');
M.pipeCoupletThreshold = get(handles.slider3,'Value');
M.pipeColorMap = getappdata(handles.output,'cMapString');
M.pipeColorHyperCube = getappdata(handles.output,'pipeColorHyperCube');
M.pipeScheme = get(handles.popupmenu4,'Value');
M.pipeUniform = get(handles.checkbox1,'Value');

N(N(end).currentVol) = M;
guidata(handles.mainHandle,N);

close(handles.output);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

ps = get(hObject,'Value');
set(handles.edit1,'String',num2str(ps));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
cMapValue = get(hObject,'Value');
allCMap = get(hObject,'String');
cMapString = allCMap{cMapValue};

eval(['cMap = ' cMapString '(200);']);
set(gcf,'ColorMap',cMap);
pipeColorHyperCube = getappdata(handles.output,'pipeColorHyperCube');
pipeColorHyperCube = setCmapCube(handles,cMap,pipeColorHyperCube);
setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
setappdata(handles.output,'cMapString',cMapString);


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% CM = cm_choose;
waitfor(cm_choose);
N = guidata(handles.mainHandle);
M = N(N(end).currentVol);
setappdata(handles.output,'CM',M.connectivityMatrix);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
maxVal = get(handles.slider1,'Max');
posit = str2num(get(handles.edit1,'String'));

if posit < maxVal && posit > 0
    set(handles.slider1, 'Value', posit);
elseif posit > maxVal
    set(handles.slider1, 'Value', maxVal);
    set(handles.edit1, 'String', num2str(maxVal));
else
    set(handles.slider1, 'Value', 0);
    set(handles.edit1, 'String', num2str(0));
end


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


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

pipeScheme = get(hObject,'Value');
setPipeScheme(handles, pipeScheme);
CM = getappdata(handles.output,'CM');
pipeColorHyperCube = getappdata(handles.output,'pipeColorHyperCube');
if ~isempty(CM)
    switch pipeScheme
        case 1
            singletCube = pipeColorHyperCube(:,:,:,1);
            if ~any(any(sum(singletCube,3))) || length(unique(sum(singletCube,3))) > 1
                singlet = rand(1,3);
            else
                singlet = squeeze(singletCube(1,1,:))';
            end
            pipeColorHyperCube = setSingletCube(handles,singlet,pipeColorHyperCube);
            setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
        case 2
            cMapString = getappdata(handles.output,'cMapString');
            eval(['cMap = ' cMapString '(200);']);
            set(gcf,'colormap', cMap);
            setappdata(handles.output,'cMapString',cMapString);
            pipeColorHyperCube = setCmapCube(handles,cMap,pipeColorHyperCube);
            setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
            cmString = get(handles.popupmenu1,'String');
            set(handles.popupmenu1,'Value',find(ismember(cmString,cMapString)));
        case 3
            cmMax = max(max(CM)); cmMin = min(min(CM));
            pct = getappdata(handles.output,'pipeCoupletThreshold');
            couplet = getappdata(handles.output,'pipeCouplet');
            if pct > cmMax || pct < cmMin
                pipeCoupletThreshold = cmMin + (cmMax - cmMin)/2;
            else
                pipeCoupletThreshold = pct;
            end
            setappdata(handles.output, 'pipeCoupletThreshold', pipeCoupletThreshold);
            pipeColorHyperCube = setCoupletCube(handles, couplet, pipeCoupletThreshold, pipeColorHyperCube);
            setappdata(handles.output,'pipeColorHyperCube',pipeColorHyperCube);
            set(handles.slider3,'Value',pipeCoupletThreshold);
            set(handles.edit2,'String',num2str(pipeCoupletThreshold));
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = uisetcolor(get(handles.pushbutton5,'BackgroundColor'));
setHueBox(V,handles.pushbutton5);
pch = getappdata(handles.output,'pipeColorHyperCube');
if get(handles.popupmenu4,'Value') == 1
    pch = setSingletCube(handles, V, pch);
elseif get(handles.popupmenu4,'Value') == 3
    couplet = getappdata(handles.output,'pipeCouplet');
    couplet(1,:) = V;
    setappdata(handles.output,'pipeCouplet',couplet);
    pch = setCoupletCube(handles, couplet, get(handles.slider3,'Value'), pch);
end
setappdata(handles.output,'pipeColorHyperCube',pch);
    

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = uisetcolor(get(handles.pushbutton6,'BackgroundColor'));
setHueBox(V,handles.pushbutton6);
pch = getappdata(handles.output,'pipeColorHyperCube');
couplet = getappdata(handles.output,'pipeCouplet');
couplet(2,:) = V;
setappdata(handles.output,'pipeCouplet',couplet);
pch = setCoupletCube(handles, couplet, get(handles.slider3,'Value'), pch);
setappdata(handles.output,'pipeColorHyperCube',pch);

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pch = getappdata(handles.output,'pipeColorHyperCube');
couplet = getappdata(handles.output,'pipeCouplet');
pch = setCoupletCube(handles, couplet, get(handles.slider3,'Value'), pch);
setappdata(handles.output,'pipeColorHyperCube',pch);
set(handles.edit2,'String',num2str(get(handles.slider3,'Value')));


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
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
retVal = str2num(get(hObject,'String'));
CM = getappdata(handles.output,'CM');
cmMin = min(min(CM)); cmMax = max(max(CM));
if retVal < cmMin
    set(handles.slider3, 'Value', cmMin);
    set(hObject,'String',num2str(cmMin));
    retVal = cmMin;
elseif retVal > cmMax
    set(handles.slider3, 'Value', cmMax);
    set(hObject,'String',num2str(cmMax));
    retVal = cmMax;
else
    set(handles.slider3,'Value',retVal);
end
pch = getappdata(handles.output,'pipeColorHyperCube');
couplet = getappdata(handles.output,'pipeCouplet');
pch = setCoupletCube(handles, couplet, retVal, pch);
setappdata(handles.output,'pipeColorHyperCube',pch);


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
