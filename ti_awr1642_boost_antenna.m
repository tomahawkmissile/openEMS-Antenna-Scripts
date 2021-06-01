close all
clear
clc

addpath('E:\openEMS');

physical_constants;
unit = 1e-3;

% Antenna geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _________________________________________________________________
%|
%|
%|
%|
%|              pad.length
%|               +-------+     +-------+     +-------+     +-------+
%| <feed.length> |       |     |       |     |       |     |       |
%|       \/      |       |     |       |     |       |     |       |
%| --------------+       +-----+       +-----+       +-----+       | 
%|                                                                 | pad.width
%| --------------+       +-----+       +-----+       +-----+       |
%|       /\      |       |     |       |     |       |     |       |
%|    feed.width |       |     |       |     |       |     |       |
%|               +-------+     +-------+     +-------+     +-------+
%|                     pad.spacing
%|
%|
%|
%|_________________________________________________________________| |
%\__________________________________________________________________\|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All units in mm
substrate.thickness = 0.1;
substrate.length = 30;
substrate.width = 30;
substrate.cells = 4;

pad.length = 0.97;
pad.width = 1.5;
pad.spacing = 1.18;
feed.width = 0.1;
feed.length = 1;

f_start = 77e9;
f0 = 79e9;
f_end = 81e9;
fc = 85e9;

substrate.epsR = 3.48;
substrate.kappa = 0.0027*2*pi*f0 * EPS0*substrate.epsR;

feed.R = 50;

SimBox = [substrate.width*2 substrate.length*2 150];

FDTD = InitFDTD('NrTs',10000);
FDTD = SetGaussExcite(FDTD,f0,fc);
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'};
FDTD = SetBoundaryCond(FDTD,BC);

CSX = InitCSX();

mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/2 SimBox(3)/2];

CSX = AddMaterial(CSX, 'substrate');
CSX = SetMaterialProperty(CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa);

start = [-substrate.width/2 -substrate.length/2 0];
stop = [substrate.width/2 substrate.length/2 substrate.thickness];

CSX = AddBox(CSX, 'substrate', 1, start, stop);
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

CSX = AddMetal(CSX, 'groundplane');
start = [-substrate.width/2 -substrate.length/2 substrate.thickness];
stop = [substrate.width/2 substrate.length/2 substrate.thickness];
CSX = AddBox(CSX, 'groundplane', 10, start, stop);

translation = [0 0 0]; %translate the antenna in space if needed
CSX = AddMetal(CSX, 'pads');

start = [0 -0.05 0]+translation;
stop = start+[1 0.1 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %feed element

start = [1 -0.75 0]+translation;
stop = start+[0.97 1.5 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %first patch

start = [1+0.97 -0.05 0]+translation;
stop = start+[1.18 0.1 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %first connector

start = [1+0.97+1.18 -0.75 0]+translation;
stop = start+[0.97 1.5 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %second patch

start = [1+0.97+1.18+0.97 -0.05 0]+translation;
stop = start+[1.18 0.1 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %second connector

start = [1+0.97+1.18+0.97+1.18 -0.75 0]+translation;
stop = start+[0.97 1.5 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %third patch

start = [1+0.97+1.18+0.97+1.18+0.97 -0.05 0]+translation;
stop = start+[1.18 0.1 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %third connector

start = [1+0.97+1.18+0.97+1.18+0.97+1.18 -0.75 0]+translation;
stop = start+[0.97 1.5 0];
CSX = AddBox(CSX, 'pads', 10, start, stop); %fourth patch

patch_mesh = DetectEdges(CSX, [], 'SetProperty', 'patch');
mesh.x = [mesh.x SmoothMeshLines(patch_mesh.x, 0.5)];
mesh.y = [mesh.y SmoothMeshLines(patch_mesh.y, 0.5)];

start = [0 -0.05 0]+translation;
stop = start+[0.01 0.1 0];
[CSX, port] = AddLumpedPort(CSX, 5, 1, feed.R, start, stop, [1 0 0], true);

mesh = DetectEdges(CSX,mesh);
mesh = SmoothMesh(mesh, c0/f0/unit/20);
CSX = DefineRectGrid(CSX, unit, mesh);

start = [mesh.x(4) mesh.y(4) mesh.z(4)];
stop = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];

[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

Sim_Path = 'tmp_patch';
Sim_CSX = 'patch.xml';

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

CSXGeomPlot([Sim_Path '/' Sim_CSX]);

RunOpenEMS( Sim_Path, Sim_CSX);



%% postprocessing & do the plots
freq = linspace(f0, f_end, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = real(0.5 * port.uf.tot .* conj( port.if.tot )); % antenna feed power

% plot feed point impedance
figure
plot( freq/1e9, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e9, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5','Mode',1);

plotFF3D(nf2ff)

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,1e-3);