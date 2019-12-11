T = random('norm',1,1,116);
regionvalues = sum(T);
connectivityMatrix = T;
customLUT='';
pipesize = 'Normal'; % 'Small','Large','XL'
nodesize = 'Large';
baseStructFileName = 'brainography-master\sample_files\ashish_116_region_base.mat'; %Need to be in local path or set absolute path in this string

%Glass brain network plot
gbplottype = 'glassbrain';
pipesonoff = 1;
figuresavename = 'gb_network';
H = brainomatic(baseStructFileName,gbplottype,figuresavename,regionvalues,connectivityMatrix,pipesonoff,customLUT,pipesize,nodesize);

%Glass brain no pipes
gbplottype = 'glassbrain';
pipesonoff = 0;
connectivityMatrix = [];
figuresavename = 'gb_nodes_only';
H = brainomatic(baseStructFileName,gbplottype,figuresavename,regionvalues,connectivityMatrix,pipesonoff,customLUT,pipesize,nodesize);

%Surface plot
gbplottype = 'surface';
pipesonoff = 0;
connectivityMatrix = [];
customLUT = 'winter';
figuresavename = 'surface_demo';
H = brainomatic(baseStructFileName,gbplottype,figuresavename,regionvalues,connectivityMatrix,pipesonoff,customLUT,pipesize,nodesize);
