function H = BrainographyADNI116(startupstruct)

%startupstruct fields:
%baseStruct the Brainography struct filename created for this run, which is
%'brainography-master/sample_files/ashish_116_region_base.mat'
%plottype: 'surface' or 'glassbrain'
%nodevalues: 1x116 size of nodes
%nodescale 'Small','Normal','Large','XL'
%pipescale (same)
%nodeflag 1/0
%pipeflag 1/0, must supply connectivity
%connectivity: 116x116 symmetric; values > 1 determine pipe radius
%surfacevalues: 1x116 for opaque surface plot
%surfaceColorMap: string of MATLAB LUT

Smallnodescale = 1.5;
Largenodescale = 3.5;
XLnodescale = 5.5;


load(startupstruct.baseStruct);

if strcmp(startupstruct.plottype,'surface')
    renderStruct(1).opacity = 1;
    renderStruct(1).singleColorFlag = 0;
    renderStruct(1).regionvalues(:,2) = num2cell(startupstruct.surfacevalues);
    renderStruct(1).pipes = 0;
    renderStruct(1).nodes = 0;
    if isempty(startupstruct.surfaceColorMap)
        colormapname = 'jet';
    else
        colormapname = startupstruct.surfaceColorMap;
    end
    rgbVals = getRGBTriple(eval([colormapname '(150)']),min(startupstruct.surfacevalues),max(startupstruct.surfacevalues),startupstruct.surfacevalues);
    renderStruct(1).regionvalues(:,3:5) = num2cell(rgbVals);
    
elseif strcmp(startupstruct.plottype,'glassbrain')
    renderStruct(1).opacity = 0.08;
    renderStruct(1).singleColorFlag = 1;
    renderStruct(1).regionvalues(:,2) = num2cell(ones(size(startupstruct.nodevalues)));
    renderStruct(1).nodeProps(:,2) = num2cell(startupstruct.nodevalues); lnth = size(renderStruct(1).nodeProps,1);
    renderStruct(1).regionvalues(:,3:5) = num2cell(repmat(renderStruct(1).singleColor,lnth,1));
    
    if ~strcmp(startupstruct.nodescale,'Normal')
        switch startupstruct.nodescale
            case 'Small'
                renderStruct(1).nodeScale = Smallnodescale; 
            case 'Large'
                renderStruct(1).nodeScale = Largenodescale;
            case 'XL'
                renderStruct(1).nodeScale = XLnodescale;
        end
    end
    
    if ~strcmp(startupstruct.pipescale,'Normal') && startupstruct.pipeflag
        switch startupstruct.pipescale
            case 'Small'
                renderStruct(1).pipeScale = 1.5;
            case 'Large'
                renderStruct(1).pipeScale = 3.5;
            case 'XL'
                renderStruct(1).pipeScale = 4.5;
        end
    else
        renderStruct(1).pipeScale = 2.5;
    end
    
    if startupstruct.pipeflag
        renderStruct(1).connectivityMatrix = startupstruct.connectivity;
        renderStruct(1).pipes = 1;
        renderStruct(1).pipeStyle = 1;
        renderStruct(1).pipeScheme = 1;
        renderStruct(1).pipeColorHyperCube = zeros([size(startupstruct.connectivity),3,3]);
    end
end

if startupstruct.saveImages
    renderStruct(2).saveImages = 1;
    renderStruct(2).figstr = startupstruct.savefilename;
else
    renderStruct(2).saveImages = 0;
end
H = figure; BrainographyRender(renderStruct,gca);