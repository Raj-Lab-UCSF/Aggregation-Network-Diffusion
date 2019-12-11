function varargout = brainography(varargin)
dbstop if error
% BRAINOGRAPHY MATLAB code for brainography.fig
%      BRAINOGRAPHY, by itself, creates a new BRAINOGRAPHY or raises the existing
%      singleton*.
%
%      H = BRAINOGRAPHY returns the handle to a new BRAINOGRAPHY or the handle to
%      the existing singleton*.
%
%      BRAINOGRAPHY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINOGRAPHY.M with the given input arguments.
%
%      BRAINOGRAPHY('Property','Value',...) creates a new BRAINOGRAPHY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before brainography_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to brainography_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help brainography

% Last Modified by GUIDE v2.5 09-Jun-2013 18:42:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brainography_OpeningFcn, ...
                   'gui_OutputFcn',  @brainography_OutputFcn, ...
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


% --- Executes just before brainography is made visible.
function brainography_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brainography (see VARARGIN)

% Choose default command line output for brainography
handles.output = hObject;
% setappdata(handles.output,'mainHandle',handles.output);
setappdata(0,'mainHandle',hObject);
% Update handles structure
% guidata(hObject, handles);
% Reserve initial guidata with "Settings" for rendering that can later be
% modified by the user and won't be deleted if the
% guidata(hObject,struct('volString','Settings','brain_at',[],'renderRes',2,'currentVol',1,'saveImages',0,'saveMovie',0));

initStruct = ui_initialize(handles);
guidata(hObject,initStruct);


function  initStruct = ui_initialize(handles)
initStruct = newRenderStruct;
initStruct.volString = 'Settings';
initStruct.renderRes = 2;
initStruct.currentVol = 1;
initStruct.saveImages = 0;
initStruct.saveMovie = 0;
initStruct.mainHandle = handles;
initStruct.figstr = 'brainography1'; 
set(handles.edit2, 'String', initStruct.figstr);
set(handles.popupmenu1,'Value',1);
set(handles.popupmenu1,'String',{'+ Add New Volume'});
defaultGUI(handles);

% set(handles.edit2,'String',initStruct.figstr);
% set(handles.checkbox1,'Enable','Off');
% set(handles.checkbox4,'Enable','Off');
% set(handles.checkbox6,'Enable','Off');
% set(handles.pushbutton14,'Enable','Off');
% set(handles.pushbutton8,'Enable','Off');
% set(handles.pushbutton10,'Enable','Off');

% UIWAIT makes brainography wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = brainography_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version ostore struct in guif MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1}=hObject;
% varargout{2}=handles; %If I can find some way to output volume struct to
% command line
% disp(hObject);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if length(handles) > 1
    disp(guihandles(hObject));
    disp(handles);
    % EXECUTE FULL-SCALE IMAGE GENERATION USING GUIDATA(1:end-1)
    figure;
    BrainographyRender(handles,gca,1);
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles(handles(end).currentVol).opacity = str2num(get(hObject,'String'));
guidata(hObject,handles);
disp(get(hObject,'String'));

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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


H = guihandles(hObject);
volVal = get(H.popupmenu1,'Value');
volString = get(H.popupmenu1,'String');

if volVal ~= size(volString,1)  %add-volume option should always be last option in menu
    % clear/delete popupmenu entry
    volString(volVal)=[];
    set(H.popupmenu1,'String',volString); %update popupmenu1
    set(H.popupmenu1,'Value',1);
    handles(volVal)=[]; % clear the removed vol from guidata/handles
    guidata(hObject,handles);
    if length(handles) == 1
        defaultGUI(H);
    else
        populateGUI(H,handles(1))
    end
elseif length(handles) == 1
    defaultGUI(H);    
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
H=guihandles(hObject);
volVal=get(H.popupmenu1,'Value');
volString=get(H.popupmenu1,'String');

if volVal == size(volString,1)
    % open file chooser, then load the volume via SPM, set/check .dim
    % property, throw error if dim does not agree and size(volString,1)~=1,
    % user sets name for vol entry, add name to top of popupmenu1 value 'String'
    % Ordering: newest vols go on top
%     if length(handles) > 1
%         handles(handles(end).currentVol).opacity=get(H.edit1,'String');
%     end
    [filename, pathname] = uigetfile({'*.nii;*.img;*.hdr','NIfTI/ANALYZE Files (*.nii,*.hdr,*.img)'}, 'Select a Volume:');
    %     disp(filterindex);
    if isequal(filename,0)
        disp('User selected Cancel')
    elseif size(filename,1) > 1
            disp('Please choose only one file at a time'); % Need better error handling
    else
        V=load_untouch_nii([pathname filesep filename]);
        % if multiple volumes per img, prompt user to choose
        if V.hdr.dime.dim(5) > 1
            imgChoice=[];
            while isempty(imgChoice) || imgChoice > V.hdr.dime.dim(5) || imgChoice < 1
                imgChoice = floor(str2num(inputdlg({['Please choose volume number (1-' num2str(length(V)) '):']},'Volume Select',1,{'1'})));
            end
        else
            imgChoice=1;
        end
        vox=flipdim(V.img(:,:,:,imgChoice),1);
        Vmat = [V.hdr.hist.srow_x; V.hdr.hist.srow_y; V.hdr.hist.srow_z; 0 0 0 1];
        Vdim = V.hdr.dime.dim(2:4);
        %if first volume chosen, this sets precedent
        if length(handles) == 1
%             disp('First Volume');
            handles(end).dim = Vdim;
            handles(end).mat = Vmat;
        end
        if ~isequal(handles(end).dim, Vdim) || ~isequal(handles(end).mat, Vmat)
            %         end
            disp('Volume dimensions or orientation do not agree with previously loaded images.');
        else
            % User prompt for name: build struct before prepending it to handles
            newVolString = inputdlg({'Enter identifier for this volume:'},'Volume Import',1,{['Volume ' num2str(length(handles))]});
            tmp=newRenderStruct;
            tmp.brain_at=vox;
            tmp.volString=newVolString{1};
            tmp.numberROI = unique(vox(find(vox~=0)));
            tmp.regionvalues = cell(size(tmp.numberROI,1),5);
            tmp.regionvalues(:,1) = num2cell(tmp.numberROI);
            tmp.regionvalues(:,2) = num2cell([1:size(tmp.numberROI,1)]');
            tmp.regionvalues(:,3:5) = num2cell(getRGBTriple(bone(1.5*size(tmp.numberROI,1)),min(tmp.numberROI),max(tmp.numberROI),tmp.numberROI));
            
            handles=[tmp handles];
            handles(end).currentVol=1;
            guidata(hObject,handles); disp('guidata set');
            volString=[{tmp.volString}; volString];
            set(hObject,'String',volString);
            set(hObject,'Value',1);
            defaultGUI(H);
        end
    end
elseif volVal ~= handles(end).currentVol % User chose to view another loaded volume, not to import new volume
%     handles(handles(end).currentVol).opacity=get(H.edit1,'String');
    handles(end).currentVol=volVal;
    guidata(hObject,handles);
    
    % populate GUI window with guidata(volVal) properties
    popVol=handles(volVal);                
    populateGUI(H,popVol);
end



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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% regionvalues;
if length(handles) > 1
    waitfor(test);  %regionvalues
    H = guihandles(hObject);
    handles = guidata(hObject);
    scValue = handles(handles(end).currentVol).singleColorFlag;
    set(H.checkbox6,'Value',scValue);
    set(H.pushbutton14,'BackgroundColor',handles(handles(end).currentVol).singleColor);
    if scValue
        set(H.pushbutton14,'Enable','On');
    else
        set(H.pushbutton14,'Enable','Off');
    end
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, ~, handles) % Launch CM chooser
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if length(handles) > 1
    waitfor(cm_choose(hObject));
    H = guihandles(hObject);
    handles = guidata(hObject);
    if ~isempty(handles(handles(end).currentVol).connectivityMatrix)
        set(H.checkbox4,'Enable','On');
        set(H.checkbox4,'Value',0);
        set(H.pushbutton10,'Enable','Off');% How to get connectivity matrix back to this GUI? With handle?
    else
        set(H.checkbox4,'Enable','Off');
        set(H.checkbox4,'Value',0);
        set(H.pushbutton10,'Enable','Off');
    end
end

% --- Executes on button press in checkbox1. Nodes on/off option
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
currentVol = handles(end).currentVol;
H = guihandles(hObject);
handles(currentVol).nodes=get(hObject, 'Value');
nodeProps = handles(currentVol).nodeProps;

if handles(currentVol).nodes
    if isempty(nodeProps);
        at = handles(currentVol).brain_at;
        numberROI = unique(at(find(at~=0)));
        nodeProps = cell(size(numberROI,1),3);
        nodeProps(:,1) = num2cell(numberROI);
        if ~isempty(handles(currentVol).connectivityMatrix)
            nodeProps(:,2) = num2cell(sum(handles(currentVol).connectivityMatrix,2));
        else
            nodeProps(:,2) = num2cell(numberROI);
        end
        nodeProps(:,3) = num2cell(ones(size(numberROI)));
        handles(currentVol).nodeProps = nodeProps;
    end
    if isempty(handles(currentVol).nodeSchema)
        handles(currentVol).nodeSchema = rand(1,3);
    end
    set(H.pushbutton8,'Enable','On');
else
    set(H.pushbutton8,'Enable','Off');
end
guidata(hObject,handles);


% --- Executes on button press in checkbox2.  Save Images Option
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
handles(end).saveImages=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox3. Save Movie Option
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
handles(end).saveMovie=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
currentVol = handles(end).currentVol;
handles(currentVol).pipes=get(hObject,'Value');
H = guihandles(hObject);

if handles(currentVol).pipes
    set(H.pushbutton10,'Enable','On');
    pipeColorHyperCube = handles(currentVol).pipeColorHyperCube;
    CM = handles(currentVol).connectivityMatrix;
    numberROI = size(CM,1);
    if isempty(pipeColorHyperCube) || numberROI ~= size(pipeColorHyperCube,1)
        handles(currentVol).pipeScheme = 1;
        pipeColorHyperCube = zeros(numberROI, numberROI, 3, 3);
        tmp = zeros(numberROI, numberROI);
        for i=1:3
            tmp(:) = rand;
            pipeColorHyperCube(:,:,i,1) = tmp;
        end
        handles(currentVol).pipeColorHyperCube = pipeColorHyperCube;
    end
else
    set(H.pushbutton10,'Enable','Off');
end
guidata(hObject,handles);

function resetMyAxes(H)
ac = allchild(H.axes1);
for i = 1:size(ac,1)
    delete(ac(i));
end
view([1 0 0]);
% clmo(handlem('light'));
if get(H.checkbox5,'Value')
    set(H.checkbox5,'Value',0);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

H = guihandles(hObject);
handles(end).saveMovie = 0;
handles(end).saveImages = 0;
%cla reset;
resetMyAxes(H);


if length(handles) > 1
    BrainographyRender(handles,H.axes1,3);
%     view([1 1 0]);
end

function populateGUI(H, popVal)
set(H.checkbox6,'Enable','On');
set(H.checkbox1,'Enable','On');


if ~isempty(popVal.opacity)
    set(H.edit1,'String',popVal.opacity);
else
    set(H.edit1,'String','1.0');
end

if popVal.nodes
    set(H.checkbox1,'Value',1);
    set(H.pushbutton8,'Enable','On');
else
    set(H.checkbox1,'Value',0);
    set(H.pushbutton8,'Enable','Off');
end

if popVal.pipes && ~isempty(popVal.connectivityMatrix)
    set(H.checkbox4,'Enable','On');
    set(H.checkbox4,'Value',1);
    set(H.pushbutton10,'Enable','On');
else
    set(H.checkbox4,'Value',0);
    set(H.pushbutton10,'Enable','Off');
    if isempty(popVal.connectivityMatrix)
        set(H.checkbox4,'Enable','Off');
    else
        set(H.checkbox4,'Enable','On');
    end
end

if popVal.singleColorFlag
    set(H.checkbox6,'Value',1);
    set(H.pushbutton14,'Enable','On');
else
    set(H.checkbox6,'Value',0);
    set(H.pushbutton14,'Enable','Off');
end

set(H.checkbox1,'Enable','On');
set(H.pushbutton14,'BackgroundColor',popVal.singleColor);



function defaultGUI(H)

set(H.edit1,'String','1.0');

set(H.checkbox1,'Enable','On');
set(H.checkbox1,'Value',0);

set(H.checkbox4,'Enable','Off');
set(H.checkbox4,'Value',0);

set(H.pushbutton10, 'Enable', 'Off');

set(H.checkbox6,'Enable','On');
set(H.checkbox6,'Value',0);

set(H.pushbutton14,'BackgroundColor',[236/255 214/255 214/255]);
set(H.pushbutton14,'Enable','Off');

resetMyAxes(H);


% function newStruct = newRenderStruct
% %Initialize new renderstruct
% newStruct=struct('volString','','brain_at',[],'dim',[],'mat',[],'opacity',1.0, ...
%     'regionvalues',[],'singleColorFlag',0,'singleColor',[236/255 214/255 214/255],'brain_colormap','bone','custom_colormap',[], ...
%     'brain_colormapidx',1,'nodes',0,'pipes',0, 'connectivityMatrix',[], ...
%     'nodeScale',2.5,'nodeSchema',[],'nodeProps',[],'nodeStyle',1,'pipeScale',1.5, ...
%     'pipeScheme',1,'pipeColorHyperCube',[],'pipeCouplet',[rand(1,3);rand(1,3)], ...
%     'pipeColorMap','jet','pipeCoupletThreshold',50,'pipeStyle',1,'pipeUniform',0,'renderRes',[], ...
%     'currentVol',[],'saveImages',0,'saveMovie',0,'figstr','','mainHandle',[],'numberROI',[]);

% --- Executes on button press in pushbutton8. Launch Nodes Advanced
% Settings
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
waitfor(node_opts);

% --- Executes on button press in pushbutton10. Launch Pipes Advanced
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles(handles(end).currentVol).connectivityMatrix) 
    waitfor(pipe_opts);
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([1 0 0]);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([0 0 1]);

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view([0 1 0]);


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
clval = get(hObject, 'Value');
H = guihandles(hObject);
% clmo(handlem('light'))
if clval
    K=camlight('left');
    setappdata(H.figure1,'clhandle',K);
else
    set(getappdata(H.figure1,'clhandle'),'Visible','Off');
end
    

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles(end).figstr = get(hObject,'String');
guidata(hObject,handles);



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


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if length(handles) > 1
    currentVol = handles(end).currentVol;
    mainHandles = guihandles(hObject);
    
    if get(hObject,'Value')
        set(mainHandles.pushbutton14,'Enable','on');
        singleColor = get(mainHandles.pushbutton14,'BackgroundColor');
        
        if isempty(handles(currentVol).regionvalues)
            brain_at = handles(currentVol).brain_at;
            uniqueRegions = unique(brain_at(find(brain_at~=0)));
            regionvalues = cell(size(uniqueRegions,1),5);
            regionvalues(:,1) = num2cell(uniqueRegions);
            regionvalues(:,2) = num2cell(ones(size(uniqueRegions,1),1));
        else
            regionvalues = handles(currentVol).regionvalues;
        end
        
        regionvalues(:,3:5) = num2cell(repmat(singleColor,size(regionvalues,1),1));
        handles(currentVol).regionvalues = regionvalues;
        handles(currentVol).singleColorFlag = 1;
    else
        handles(currentVol).singleColorFlag = 0;
        set(mainHandles.pushbutton14,'Enable','off');
    end
    guidata(hObject,handles);
end
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = uisetcolor(get(hObject,'BackgroundColor'));
set(hObject,'BackgroundColor',V);

currentVol = handles(end).currentVol;
handles(currentVol).singleColor = V;
regionvalues = handles(currentVol).regionvalues;
numberROI = size(regionvalues,1);
regionvalues(:,3:5) = num2cell(repmat(V,numberROI,1));
handles(currentVol).regionvalues = regionvalues;
guidata(hObject,handles);


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    renderstruct = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    if length(renderstruct) == 1
        matIn = load([pathname filesep filename]);
        eval(['structIn = matIn.' renderstruct.name ';']);
        [tf,tmpStruct] = convertRenderStruct(structIn);
        
        if tf% && (length(structIn > 1))
            H = guihandles(hObject);
            structIn = tmpStruct;
            guidata(hObject,structIn);
            currentVol = structIn(end).currentVol;
            populateGUI(H, structIn(currentVol));
            volNameList = {'+ Add New Volume'};
            for i=length(structIn)-1:-1:1
                volNameList = [structIn(i).volString; volNameList];
            end
            set(H.popupmenu1,'String',volNameList);
            set(H.popupmenu1,'Value',currentVol);
        else
            disp('Not a proper Brainography scene struct.');
        end
    end
end

function rsTF = isRenderStruct(structIn)
rsTF = 1;
reqFields = fields(newRenderStruct);
incFields = fields(structIn);
for i=1:length(reqFields)
    if ~ismember(reqFields{i},incFields)
        rsTF = 0;
        break;
    end
end


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.mat'});
renderStruct = handles;
save([pathname filesep filename],'renderStruct');
disp(['Scene file saved to ' pathname filesep filename]);


% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H = guihandles(hObject);
initStruct = ui_initialize(H);
guidata(hObject, initStruct);
set(H.popupmenu1,'Value',1);
set(H.popupmenu1,'String',{'+ Add New Volume'});
defaultGUI(H);
