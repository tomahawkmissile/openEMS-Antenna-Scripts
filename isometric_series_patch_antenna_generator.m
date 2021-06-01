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

% Script assumes openEMS is already added to path.

% Struct definitions:

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
% - pad.spacing = the series spacing on a single antenna element between
% pads
% - pad.interconnect_length = the length of the connection between pads

% Feed definition (the input)
% - feed.R = the input impedance
% - feed.length = the length of the feed. This is often the same as the pad
% interconnect width
% - feed.length = the length of the feed.

% Antenna parameters
% - series_count = number of series elements per patch element
% - parallel_count = number of patch elements
% - fc = center frequency

% End struct definitions

function CSX = series_patch_generate(CSX, antenna_struct)

    

end