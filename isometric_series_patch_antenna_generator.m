% This script is an openEMS geometry generator for isometric series patch
% antennas. This script will have very verbose comments.

% Antenna geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _________________________________________________________________
%|
%| Supports multiple antennas (for MIMO beamforming)
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
%|_________________________________________________________________| |
%\__________________________________________________________________\|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Y+
%|      Z-
%|     /
%|    /
%|   /
%|  /
%| /
%|/_________________ X+

% Script assumes openEMS is already added to path.

% Note: CSX parameter is an openEMS CSX object
% Note: The coordinate system here is defined as where the antenna is
% sideways and the series antenna points in the X+ direction. The feed
% lines are on the left side (X-) and the antenna array is arranged in the
% Y axis. The Z direction is above/below the antenna plane.
% Length is on the X axis, width is on the Y, and thickness is Z.

% Struct definitions (parameter "as"):

% Substrate definition:
% - substrate.thickness = thickness of the substrate (PCB laminate)
% - substrate.length = length of the substrate (PCB length)
% - substrate.width = width of the substrate (PCB length)
% - substrate.cells = number of cells to put inside the substrate
% (thickness)
% - substrate.epsR = substrate dielectric constant
% - substrate.kappa = Defined as: (substrate loss tangent)*(2*pi)*(center
% frequency)*(EPS0 which is vacuum permittivity)*(substrate dielectric
% constant)

% Pad definition (a pad is defined as the large square sections of the patch antenna)
% - pad.length = the length of the pad
% - pad.width = the width of the pad
% - pad.interconnect_length = the length of the connection between pads. 
% the series spacing on a single antenna element between pads.
% - pad.interconnect_width = the width of the connection between pads

% Feed definition (the input)
% - feed.R = the input impedance
% - feed.length = the length of the feed. This is often the same as the pad
% interconnect width
% - feed.width = the width of the feed.

% Antenna parameters
% - series_count = number of series elements per patch element
% - parallel_count = number of patch elements
% - parallel_spacing = the spacing between each array. defined by
% the gap between the center of the elements.
% - f_start = start frequency
% - f_end = end frequency
% - cu_thickness = conductor thickness. allowed to be 0.
% - start_pos = start position in 3D space. enter as a vector. this
% position denotes where the bottom left of the bottom patch element's feed
% line starts. The antenna may extend below this point.
% - unit = unit magnitude relative to a meter

% End struct definitions

function [CSX, port, mesh] = isometric_series_patch_antenna_generator(as)

    physical_constants;

    thickness = as.cu_thickness;
    f_start = as.f_start;
    f_end = as.f_end;
    
    series = as.series_count;
    parallel = as.parallel_count;
    
    feed_length = as.feed.length;
    feed_width = as.feed.width;
    R = as.feed.R;
    
    pad_length = as.pad.length;
    pad_width = as.pad.width;
    interconnect_length = as.pad.interconnect_length;
    interconnect_width = as.pad.interconnect_width;
    parallel_spacing = as.parallel_spacing;
    
    substrate_width = as.substrate.width;
    substrate_length = as.substrate.length;
    substrate_thickness = as.substrate.thickness;
    substrate_cells = as.substrate.cells;
    kappa = as.substrate.kappa;
    epsR = as.substrate.epsR;
    
    translation = as.start_pos;

    SimBox = [substrate_width*2 substrate_length*2 200];
    
    max_res = c0/(f_start+f_end)/as.unit/20;
    CSX = InitCSX();
    
    mesh.x = [-SimBox(1)/2 SimBox(1)/2 -substrate_length/2 substrate_length/2];
    mesh.y = [-SimBox(2)/2 SimBox(2)/2 -substrate_width/2 substrate_width/2];
    mesh.z = [-SimBox(3)/2 linspace(0,substrate_thickness,substrate_cells), SimBox(3)/2];
    mesh.z = SmoothMeshLines(mesh.z, max_res, 1.4);
    
    mesh.x = SmoothMeshLines(mesh.x, max_res, 1.4);
    mesh.y = SmoothMeshLines(mesh.y, max_res, 1.4);
    
    mesh = AddPML(mesh, [8 8 8 8 8 8]);
    CSX = DefineRectGrid(CSX, as.unit, mesh);
    
    CSX = AddMaterial(CSX, 'substrate');
    CSX = SetMaterialProperty(CSX, 'substrate', 'Epsilon', epsR, 'Kappa', kappa);

    start = [-substrate_width/2 -substrate_length/2 0];
    stop = [substrate_width/2 substrate_length/2 substrate_thickness];
    CSX = AddBox(CSX, 'substrate', 1, start, stop);

    CSX = AddMetal(CSX, 'groundplane');
    start = [-substrate_width/2 -substrate_length/2 substrate_thickness];
    stop = [substrate_width/2 substrate_length/2 substrate_thickness];
    CSX = AddBox(CSX, 'groundplane', 10, start, stop);
    
    CSX = AddMetal(CSX, 'pads');
    
    for p = 0:parallel-1
        % Feed geometry
        start = [0 p*parallel_spacing 0] + translation;
        stop = start + [feed_length feed_width thickness];
        CSX = AddBox(CSX, 'pads', 10, start, stop);
        [CSX port{p+1}] = AddLumpedPort(CSX, 999, p+1, R, start, start+[0.01 feed_width thickness], [1 0 0], true);
        
        for s = 0:series-1
            % Interconnect
            start = [feed_length+((s)*pad_length)+(s*interconnect_length) (feed_width/2-interconnect_width/2)+p*parallel_spacing 0]+translation;
            stop = start+[interconnect_length interconnect_width thickness];
            CSX = AddBox(CSX, 'pads', 10, start, stop);
            
            % Pad element
            start = [feed_length+(s*pad_length+(s+1)*interconnect_length) (feed_width/2-pad_width/2)+p*parallel_spacing 0]+translation;
            stop = start+[pad_length pad_width thickness];
            CSX = AddBox(CSX, 'pads', 10, start, stop);
        end
    end
    
end