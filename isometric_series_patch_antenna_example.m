% Generate a patch antenna array that is 4 parallel elements with 5 series elements

close all
clear
clc

physical_constants;
unit = 1e-3;

as.substrate.thickness = 1.6;
as.substrate.length = 100;
as.substrate.width = 50;
as.substrate.cells = 4;
as.substrate.epsR = 4.4;
as.substrate.kappa = 1;

as.pad.length = 0.97;
as.pad.width = 1.5;
as.pad.interconnect_length = 1.18;
as.pad.interconnect_width = 0.1;

as.feed.R = 50;
as.feed.length = 1;
as.feed.width = 0.1;

as.series_count = 8;
as.parallel_count = 8;
as.parallel_spacing = 2.5;
as.fc = 79e9;
as.cu_thickness = 0.035;
as.start_pos = [0 0 0];

SimBox = [as.substrate.width*2 as.substrate.length*2 200];

CSX = InitCSX();
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];

[CSX ports] = isometric_series_patch_antenna_generator(CSX,as);

patch_mesh = DetectEdges(CSX, [], 'SetProperty', 'patch');
mesh.x = [mesh.x SmoothMeshLines(patch_mesh.x, 0.5)];
mesh.y = [mesh.y SmoothMeshLines(patch_mesh.y, 0.5)];

mesh = DetectEdges(CSX,mesh);
mesh = SmoothMesh(mesh, c0/as.fc/unit/20);
CSX = DefineRectGrid(CSX, unit, mesh);

FDTD = InitFDTD();

Sim_Path =  'tmp_generator';
Sim_CSX = 'generated.xml';

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

CSXGeomPlot([Sim_Path '/' Sim_CSX]);