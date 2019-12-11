function H = brainomatic(baseStructFileName,atlas_type, gbplottype,figuresavename,regionvalues,connectivityMatrix,pipesonoff,customLUT,pipesize,nodesize)

% baseStructFileName:
%    'brainography-master/sample_files/ashish_116_region_base.mat'
% gbplottype: 'surface' or 'glassbrain' for node/pipe
% figuresavename
% regionvalues: can be for either node sizes or surface values, determined
%    by gbplottype
% connectivityMatrix: 116x116 (optional)
% pipesonoff 1/0 (opt)
% customLUT: for surface plots (opt) 'jet','spring','winter',etc
% pipesize: 'Small','Normal','Large','XL'
% nodesize: same
% added by Raj 2/5/14: removed rescaling of regionvalues by max(abs(regionvalues))
% its the user's responsibility to rescale regionvalues to [0,1] or
% something predictable

if nargin < 8
    pipesize = 'Normal';
    nodesize = 'Normal';
    if nargin < 5
        connectivityMatrix = '';
        pipesonoff = 0;
        customLUT = 'jet';
    end
end
    
startupstruct.baseStruct = baseStructFileName;
startupstruct.plottype = gbplottype;

if ~isempty(figuresavename)
    startupstruct.savefilename = figuresavename;
    startupstruct.saveImages = 1;
else
    startupstruct.saveImages = 0;
end

if strcmp('surface',startupstruct.plottype)
    startupstruct.surfacevalues = regionvalues;
    if ~isempty(customLUT)
        startupstruct.surfaceColorMap = customLUT;
    else
        startupstruct.surfaceColorMap = '';
    end
else
    startupstruct.nodevalues = regionvalues;
    startupstruct.nodeflag = 1;
end

if ~isempty(connectivityMatrix) && pipesonoff
    startupstruct.connectivity = connectivityMatrix;
    startupstruct.pipeflag = 1;
else
    startupstruct.pipeflag = 0;
end

startupstruct.pipescale = pipesize;
startupstruct.nodescale = nodesize;

if strcmp(atlas_type, 'MS86')
    H = BrainographyMS86(startupstruct);
elseif strcmp(atlas_type, 'ADNI116')
    H = BrainographyADNI116(startupstruct);
end