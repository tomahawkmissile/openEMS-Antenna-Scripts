% Generate a patch antenna array that is 4 parallel elements with 5 series elements

close all
clear
clc

physical_constants;
unit = 1e-3;
as.unit=unit;

as.series_count = 8;
as.parallel_count = 8;
as.parallel_spacing = 2.5;
as.f_start = 77e9;
as.fc = 79e9;
as.f_end = 81e9;
as.cu_thickness = 0.035;
as.start_pos = [0 0 0];

as.substrate.thickness = 0.1;
as.substrate.length = 120;
as.substrate.width = 100;
as.substrate.cells = 4;
as.substrate.epsR = 4.4;
as.substrate.kappa = 0.0027*2*pi*as.fc * EPS0*as.substrate.epsR;

as.pad.length = 0.97;
as.pad.width = 1.5;
as.pad.interconnect_length = 1.18;
as.pad.interconnect_width = 0.1;

as.feed.R = 50;
as.feed.length = 1;
as.feed.width = 0.1;

FDTD = InitFDTD('NrTs',10000);
FDTD = SetGaussExcite(FDTD,as.f_start,as.f_end);
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'};
FDTD = SetBoundaryCond(FDTD,BC);

[CSX, port, mesh, SimBox] = isometric_series_patch_antenna_generator(as);

%start = [mesh.x(4) mesh.y(4) mesh.z(4)];
%stop = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', -SimBox/2, SimBox/2);

Sim_Path =  'C:\Users\tomah\OneDrive\Documents\MATLAB\tmp';
Sim_CSX = 'generated.xml';

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

CSXGeomPlot([Sim_Path '/' Sim_CSX]);

RunOpenEMS(Sim_Path, Sim_CSX);

