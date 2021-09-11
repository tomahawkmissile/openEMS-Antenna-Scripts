% Note: this script needs a lot of RAM!!

close all
clear
clc 

physical_constants;
unit = 1e-3;

series_count = 4;
parallel_count = 1;
parallel_spacing = 2.5;

f_start = 77e9;
fc = 79e9;
f_end = 81e9;
cu_thickness = 0.035;

translation = [0 0 0];

substrate.thickness = 0.1;
substrate.length = 30;
substrate.width = 30;
substrate.cells = 4;
substrate.epsR = 3.48;
substrate.kappa = 0.0027*2*pi*fc * EPS0*substrate.epsR;

pad.length = 0.97;
pad.width = 1.5;
pad.interconnect.length = 1.18;
pad.interconnect.width = 0.1;

feed.R = 50;
feed.length = 1;
feed.width = 0.1;

FDTD = InitFDTD('NrTs',50000);
FDTD = SetGaussExcite(FDTD,f_start,f_end);
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'};
FDTD = SetBoundaryCond(FDTD,BC);

SimBox = [substrate.length*2 substrate.width*2 25];
    
max_res = c0/f_end/unit/32;
CSX = InitCSX();

mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [0 SimBox(3)];
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];
mesh.z = [mesh.z linspace(cu_thickness+substrate.thickness,SimBox(3),(SimBox(3)-cu_thickness+substrate.thickness)/max_res)];

CSX = AddMaterial(CSX, 'substrate');
CSX = SetMaterialProperty(CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa);

start = [-substrate.width/2 -substrate.length/2 0];
stop = [substrate.width/2 substrate.length/2 substrate.thickness];
CSX = AddBox(CSX, 'substrate', 1, start, stop);

CSX = AddMetal(CSX, 'groundplane');
start = [-substrate.width/2 -substrate.length/2 0];
stop = [substrate.width/2 substrate.length/2 0];
CSX = AddBox(CSX, 'groundplane', 10, start, stop);

CSX = AddMetal(CSX, 'pads');

for p = 0:parallel_count-1
    % Feed geometry
    start = [0 p*parallel_spacing substrate.thickness] + translation;
    stop = start + [feed.length feed.width cu_thickness];
    CSX = AddBox(CSX, 'pads', 10, start, stop);
    port_length = 0.01;
    [CSX port{p+1}] = AddLumpedPort(CSX, 999, p+1, feed.R, start, start+[port_length feed.width cu_thickness], [0 0 1], true);
    mesh.x = [mesh.x 0 port_length]; %Add cells to capture port
    
    mesh.x = [mesh.x SmoothMeshLines([-SimBox(1)/2,feed.length], max_res, 1.4)]; %Smooth mesh from start of x domain up to end of feed
    
    for s = 0:series_count-1
        % Interconnect
        start = [feed.length+((s)*pad.length)+(s*pad.interconnect.length) (feed.width/2-pad.interconnect.width/2)+p*parallel_spacing substrate.thickness]+translation;
        stop = start+[pad.interconnect.length pad.interconnect.width cu_thickness];
        mesh.x = [mesh.x linspace(start(1)+max_res*2/3, stop(1)-max_res*2/3, 10)];
        mesh.y = [mesh.y linspace(start(2), stop(2), 4)];
        
        CSX = AddBox(CSX, 'pads', 10, start, stop);

        % Pad element
        start = [feed.length+(s*pad.length+(s+1)*pad.interconnect.length) (feed.width/2-pad.width/2)+p*parallel_spacing substrate.thickness]+translation;
        stop = start+[pad.length pad.width cu_thickness];
        mesh.x = [mesh.x start(1)-max_res*2/3 stop(1)+max_res*2/3]; %Generate lines parallel to pad boundaries obeying 1/3 rule
        mesh.y = [mesh.y start(2)-max_res*2/3 stop(2)+max_res*2/3]; %Generate lines parallel to pad boundaries obeying 1/3 rule
        mesh.x = [mesh.x linspace(start(1)+max_res*1/3 , stop(1)-max_res*1/3 , (stop(1)-start(1))/max_res)]; %Add cells in between x
        mesh.y = [mesh.y linspace(start(2)+max_res*1/3 , stop(2)-max_res*1/3 , (stop(2)-start(2))/max_res)]; %Add cells in between y
        CSX = AddBox(CSX, 'pads', 10, start, stop);
    end
    mesh.x = [mesh.x SmoothMeshLines([stop(1)+max_res*2/3,SimBox(1)/2], max_res, 1.4)]; %Smooth x domain from last cell of pad to end of sim box
    mesh.y = [mesh.y SmoothMeshLines([-SimBox(2)/2 start(2)-max_res*2/3], max_res, 1.4)]; %Smooth y domain from start of sim box to first cell touching a pad
    mesh.y = [mesh.y SmoothMeshLines([stop(2)+max_res*2/3,SimBox(2)/2], max_res, 1.4)]; %Smooth y domain from last pad cell to end of sim box
end

%patch_mesh = DetectEdges(CSX, [], 'SetProperty', 'pads');
%mesh.x = [mesh.x SmoothMeshLines([-SimBox(1)/2 SimBox(1)/2], max_res, 1.4)];
%mesh.y = [mesh.y SmoothMeshLines([-SimBox(2)/2 SimBox(2)/2], max_res, 1.4)];

%mesh = DetectEdges(CSX, mesh);
%mesh = SmoothMesh(mesh, max_res);
CSX = DefineRectGrid(CSX, unit, mesh);

start = [mesh.x(4) mesh.y(4) mesh.z(4)];
stop = [mesh.x(end-3) mesh.y(end-3) mesh.z(end-3)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

Sim_Path =  'C:\Users\tomah\OneDrive\Documents\MATLAB\out\tmp';
Sim_CSX = 'generated.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

CSXGeomPlot([Sim_Path '/' Sim_CSX]);

RunOpenEMS(Sim_Path, Sim_CSX);

%% postprocessing & do the plots
freq = linspace(f_start, f_end, 501 );
U = ReadUI( {'port_ut1','et'}, 'tmp/', freq ); % time domain/freq domain voltage
I = ReadUI( 'port_it1', 'tmp/', freq ); % time domain/freq domain current (half time step is corrected)

% plot time domain voltage
figure
[ax,h1,h2] = plotyy( U.TD{1}.t/1e-9, U.TD{1}.val, U.TD{2}.t/1e-9, U.TD{2}.val );
set( h1, 'Linewidth', 2 );
set( h1, 'Color', [1 0 0] );
set( h2, 'Linewidth', 2 );
set( h2, 'Color', [0 0 0] );
grid on
title( 'time domain voltage' );
xlabel( 'time t / ns' );
ylabel( ax(1), 'voltage ut1 / V' );
ylabel( ax(2), 'voltage et / V' );
% now make the y-axis symmetric to y=0 (align zeros of y1 and y2)
y1 = ylim(ax(1));
y2 = ylim(ax(2));
ylim( ax(1), [-max(abs(y1)) max(abs(y1))] );
ylim( ax(2), [-max(abs(y2)) max(abs(y2))] );

% plot feed point impedance
figure
Zin = U.FD{1}.val ./ I.FD{1}.val;
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
uf_inc = 0.5*(U.FD{1}.val + I.FD{1}.val * 50);
if_inc = 0.5*(I.FD{1}.val - U.FD{1}.val / 50);
uf_ref = U.FD{1}.val - uf_inc;
if_ref = I.FD{1}.val - if_inc;
s11 = uf_ref ./ uf_inc;
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

%%
number = 1;
P_in = 0;
for xn=1:parallel_count
    U = ReadUI( ['port_ut' int2str(number)], 'tmp/', freq ); % time domain/freq domain voltage
    I = ReadUI( ['port_it' int2str(number)], 'tmp/', freq ); % time domain/freq domain current (half time step is corrected)

    P_in = P_in + 0.5*U.FD{1}.val .* conj( I.FD{1}.val );
    number=number+1;
end

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = (0:2:359) - 180;
phiRange = [0 90];
r = 1; % evaluate fields at radius r
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180);

Dlog=10*log10(nf2ff.Dmax);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(Dlog) ' dBi'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

% display phi
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2]);
drawnow


%% calculate 3D pattern
phiRange = 0:3:360;
thetaRange = unique([0:0.5:15 10:3:180]);
disp( 'calculating 3D far field...' );
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180, 'Verbose',2,'Outfile','nf2ff_3D.h5');
figure
plotFF3D(nf2ff);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,1e-3);

%% visualize magnetic fields
% you will find vtk dump files in the simulation folder (tmp/)
% use paraview to visulaize them