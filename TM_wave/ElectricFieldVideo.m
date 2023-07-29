addpath("TM_wave/TFSF/");
addpath("TM_wave/PML/");

clc
clear

% Make a video for point source
cylinder_options = [-3, 1, 3.4, 1.2, 0]; % [x0, R, e_r, sigma_r, y0]
simulation_options = [10, 10, 100, 100, 12, 10*10^9]; % [X0, Y0, Tn, f0]
boundary = "PML"; % Type of boundary condition
boundary_case = "full"; % boundaries only in the left side
PML_options = [16, 2, 10^(-6)]; % [Npml, power, R]

[Ez_p, Hx_p, Hy_p] = CylinderScattering(cylinder_options, simulation_options, ...
    boundary, boundary_case, PML_options);

writerObj = VideoWriter('ElectricFieldPoint.avi');
writerObj.FrameRate = 15; 
open(writerObj);

figure(1);
for k = 1:size(Ez_p,3)
    pcolor(Ez_p(:, :, k));
    shading interp;
    colormap default;
    
    xlabel('X-axis');
    ylabel('Y-axis');
    
    % Capture the frame and add it to the video
    currFrame = getframe(gcf);
    writeVideo(writerObj, currFrame);
end

close(writerObj);


% Make a video for plane wave
cylinder_options_sc = [0, 1, 3.4, 1.2, 0];
simulation_options = [10, 10, 100, 100, 8, 10*10^9]; % [X0, Y0, N_x, N_y Tn, f0]
PML_options = [16, 2, 10^(-6)]; % [Npml, power, R]
TFSF_options = [8, pi/4];

[Ez_sc, ~, ~] = TFSF_formulation(cylinder_options_sc, simulation_options, ...
    PML_options, TFSF_options);

writerObj = VideoWriter('ElectricFieldPlane.avi');
writerObj.FrameRate = 15; 
open(writerObj);

figure(2);
for k = 1:size(Ez_sc,3)
    pcolor(Ez_sc(:, :, k));
    shading interp;
    colormap default;
    
    xlabel('X-axis');
    ylabel('Y-axis');
    
    % Capture the frame and add it to the video
    currFrame = getframe(gcf);
    writeVideo(writerObj, currFrame);
end

close(writerObj);