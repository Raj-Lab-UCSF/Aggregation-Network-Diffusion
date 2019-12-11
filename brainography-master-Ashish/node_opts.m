function varargout = node_opts(varargin)
% NODE_OPTS MATLAB code for node_opts.fig
%      NODE_OPTS, by itself, creates a new NODE_OPTS or raises the existing
%      singleton*.
%
%      H = NODE_OPTS returns the handle to a new NODE_OPTS or the handle to
%      the existing singleton*.
%
%      NODE_OPTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NODE_OPTS.M with the given input arguments.
%
%      NODE_OPTS('Property','Value',...) creates a new NODE_OPTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before node_opts_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to node_opts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help node_opts

% Last Modified by GUIDE v2.5 09-May-2013 17:07:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @node_opts_OpeningFcn, ...
                   'gui_OutputFcn',  @node_opts_OutputFcn, ...
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


% --- Executes just before node_opts is made visible.
function node_opts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to node_opts (see VARARGIN)

% Choose default command line output for node_opts
handles.output = hObject;
mainHandle = getappdata(0,'mainHandle');
handles.mainHandle=mainHandle;
% Update handles structure
guidata(hObject, handles);

N = guidata(handles.mainHandle);
M = N(N(end).currentVol);

% if ~isempty(M.nodeScale)
disp(M.nodeScale);
set(handles.slider1, 'Min', 0, 'Max', 20.0, 'Value', M.nodeScale);
set(handles.edit4, 'String', num2str(M.nodeScale));
% else
%     set(handles.slider1,'Value',1.5);
%     set(handles.edit4,'String','1.5');
% end

if isempty(M.nodeProps) && ~isempty(M.brain_at)
    U = unique(M.brain_at(M.brain_at~=0));
    T = cell(size(U,1), 3);
    T(:,1) = num2cell(U);
    set(handles.uitable1, 'Data', T);
elseif ~isempty(M.nodeProps) && ~isempty(M.brain_at)
    set(handles.uitable1, 'Data', M.nodeProps);
end

if ~isempty(M.nodeSchema)
    setSchema(M.nodeSchema,handles);
else
    setSchema(defaultSchema,handles);
end

set(handles.popupmenu3, 'Value', M.nodeStyle);
    
% UIWAIT makes node_opts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = node_opts_OutputFcn(hObject, eventdata, handles) 
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
sV = get(handles.slider1, 'Value');
N = guidata(handles.mainHandle);
M = N(N(end).currentVol);
M.nodeScale = sV;
M.nodeProps = get(handles.uitable1, 'Data');
M.nodeStyle = get(handles.popupmenu3, 'Value');
M.nodeSchema = getappdata(handles.output, 'nodeSchema');

if size(M.nodeSchema,1) < max(cell2mat(M.nodeProps(:,3)))
    diffSch = max(cell2mat(M.nodeProps(:,3))) - size(M.nodeSchema,1);
    appSch = rand(diffSch,3);
    M.nodeSchema = [M.nodeSchema; appSch];
    disp('Correcting node schema');
end
    

N(N(end).currentVol) = M;
guidata(handles.mainHandle, N);

close(handles.output);

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% setappdata(handles.output,'sliderValue',get(hObject,'Value'));
set(handles.edit4,'String',num2str(get(hObject,'Value')));

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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
maxVal = get(handles.slider1,'Max');
posit = str2num(get(handles.edit4,'String'));

if posit < maxVal && posit > 0
    set(handles.slider1, 'Value', posit);
elseif posit > maxVal
    set(handles.slider1, 'Value', maxVal);
    set(handles.edit4, 'String', num2str(maxVal));
else
    set(handles.slider1, 'Value', 0);
    set(handles.edit4, 'String', num2str(0));
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

retVal = round(str2double(get(hObject,'String')));

if ~isnan(retVal) && retVal > 0
    set(hObject,'String',num2str(retVal));
    set(handles.popupmenu4,'Value',1);
    set(handles.popupmenu4,'String',cellstr(num2str([1:retVal]')));
    newSchema = rand(retVal,3);
else
    set(hObject,'String','1');
    set(handles.popupmenu4,'Value',1);
    newSchema = rand(1,3);
end

setappdata(handles.output,'nodeSchema', newSchema);
setSchema(newSchema,handles);
setHueBox(newSchema(1,:),handles.pushbutton5);



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
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
sc = getappdata(handles.output,'nodeSchema');
idx = get(hObject,'Value');
setHueBox(sc(idx,:),handles.pushbutton5);

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
setHueBox(V,hObject);
nob = get(handles.popupmenu4,'Value');
sc = getappdata(handles.output,'nodeSchema');
sc(nob,:) = V;
setappdata(handles.output,'nodeSchema',sc);

function setHueBox(RGB,boxHandle)
set(boxHandle,'BackgroundColor',RGB);
set(boxHandle,'ForegroundColor',RGB);

function setSchema(nodeSchema,handles)
set(handles.edit5,'String',num2str(size(nodeSchema,1)));
set(handles.popupmenu4,'String',cellstr(num2str([1:size(nodeSchema,1)]')));
setappdata(handles.output,'nodeSchema',nodeSchema);
setHueBox(nodeSchema(1,:),handles.pushbutton5);

function lobeSchema = defaultSchema
% colors = {'blue';  'magenta'; 'green'; 'red'; 'cyan'; 'yellow'};
lobeSchema = zeros(6,3);
lobeSchema(1,:) = [0 0 1];
lobeSchema(2,:) = [1 0 1];
lobeSchema(3,:) = [0 1 0];
lobeSchema(4,:) = [1 0 0];
lobeSchema(5,:) = [0 1 1];
lobeSchema(6,:) = [1 1 0];

        


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    T = get(handles.uitable1,'Data');
    lobe_ids = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    if size(lobe_ids,1) < 2
        matIn = load([pathname filesep filename]);
        eval(['LUT = matIn.' lobe_ids.name ';']);
        if length(LUT) ~= size(T,1)
            disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(LUT))]);
        else
            if size(LUT,2) == size(T,1)
                LUT = LUT';
            end
            T(:,3) = num2cell(LUT);
            set(handles.uitable1,'Data',T);
        end
    else
        disp('Please specify file containing single variable.');
    end
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
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
            T(:,3) = num2cell(RV);
            set(handles.uitable1,'Data',T);
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
        eval(['LUT = matIn.' regionvar.name ';']);
        if length(RV) ~= size(T,1)
            disp(['Imported mat wrong size: expecting ' num2str(size(T,1)) ', found ' num2str(length(RV))]);
        else
            if size(RV,2) == size(T,1)
                RV = RV';
            end
            T(:,2) = num2cell(RV);
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
            T(:,2) = num2cell(RV);
            set(handles.uitable1,'Data',T);
        end
    end
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');
T(:,3) = cell(length(T),1);
set(handles.uitable1,'Data',T);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T = get(handles.uitable1,'Data');
T(:,2) = cell(length(T),1);
set(handles.uitable1,'Data',T);
