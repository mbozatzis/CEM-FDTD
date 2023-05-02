clc
clear

cylinder_options = [3, 1, 3.4, 1.2];
simulation_options = [10, 10, 10, 10*10^9];
boundary = "PML";
PML_options = [8, 2, 10^(-6)];

% To-dos
% CEM: 
% - Matrices with fields
% - Case for PML only at one side
% - How to make videos
% - Start making main (plots comparing boundaries, comparing different
% cylinders, time, PML options, put cylinder in different places, maybe an
% angle diagram for the boundary conditions, show PML plot, show stability and
% make it unstable in purpose, velocity plot (kdx and angle), maybe frequency?)
% - report with theoretical calculations
%
% Scattering:
% - TF/SF Formulation implementation
% - PML at the corners
% - Analytical solution implementation
% - Simulation with different cylinders at different positions and angles
% - Add complex materials

[Ez, Hx, Hy] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary,  PML_options);