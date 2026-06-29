%% s19_MeMe_roll_stimuli
% Generate the rolling visual stimuli used by the MeMe network simulation
% (Figure S5D-H). Running this script writes two stimuli to Processed_Data\:
%
%   Roll_60deg_cw.mat  / .mp4    (clockwise roll)
%   Roll_60deg_ccw.mat / .mp4    (counter-clockwise roll)
%
% Pattern: a "sky roll" -- after an exact spherical rotation, pixels whose
% rotated world-z (Zw) is > 0 (the upper / sky hemisphere) are bright (1) and
% the rest (the lower / ground hemisphere) are dark (0). Rotation is at a
% constant angular speed.
%
% The .mat content (pat: arenaheight x arenawidth x nFrames, double 0/1) is the
% stimulus input read by Figures\fig_S5D_E_F_G_H_MeMe_simulation.m.
%
% To change the roll angle, edit total_roll_deg below (it is reflected in the
% file names, e.g. Roll_90deg_cw).
clear all; close all; clc

baseDir = fileparts(fileparts(mfilename('fullpath')));  % resolved from this script's location (portable; no absolute paths)
outDir  = fullfile(baseDir, 'Processed_Data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% --- Settings ---
total_roll_deg = 60;     % total roll angle (deg). File name = Roll_<deg>deg_{cw,ccw}
deg_per_px     = 0.5;    % resolution (deg/pixel)
fps            = 50;     % frame rate
roll_dur       = 2.0;    % sec  roll duration (angular speed = total_roll_deg / roll_dur deg/s)
pre_dur        = 1.0;    % sec  hold before roll (theta = 0)
post_dur       = 1.0;    % sec  hold after roll  (theta = total_roll)

%% --- Resolution / timing ---
arenawidth  = round(360 / deg_per_px);   % 720
arenaheight = round(180 / deg_per_px);   % 360
total_roll  = total_roll_deg * pi/180;   % rad

n_pre   = round(pre_dur  * fps);
n_roll  = round(roll_dur * fps);
n_post  = round(post_dur * fps);
nFrames = n_pre + n_roll + n_post;

% Constant angular speed: pre (hold 0) -> roll (0..total_roll linear) -> post (hold total_roll)
theta_seq = [zeros(1, n_pre), ...
             linspace(0, total_roll, n_roll), ...
             total_roll * ones(1, n_post)];

%% --- Fly-frame pixel direction vectors ---
%   az: azimuth (phi),  el: elevation (theta)
az = (-180 + deg_per_px/2 : deg_per_px : 180 - deg_per_px/2) * pi/180;
el = (-90  + deg_per_px/2 : deg_per_px :  90 - deg_per_px/2) * pi/180;
[AZ, EL] = meshgrid(az, el);
Yf = cos(EL) .* sin(AZ);
Zf = sin(EL);

roll_tag = sprintf('%ddeg', round(total_roll_deg));

%% --- Build the two roll stimuli (cw / ccw) ---
%   sky pattern: pixels with rotated world-z (Zw) > 0 (upper/sky) = 1, rest = 0.
%   sgn = +1 (cw), -1 (ccw).
cases = {
    ['Roll_' roll_tag '_cw'],  +1;
    ['Roll_' roll_tag '_ccw'], -1;
};

for ci = 1:size(cases, 1)
    name = cases{ci, 1};
    sgn  = cases{ci, 2};

    pat = zeros(arenaheight, arenawidth, nFrames);   % double
    for k = 1:nFrames
        th = sgn * theta_seq(k);
        Zw = sin(th) * Yf + cos(th) * Zf;            % exact rotation (no small-angle)
        frame = zeros(arenaheight, arenawidth);
        frame(Zw > 0) = 1;                           % sky: upper bright (1), lower dark (0)
        pat(:, :, k) = frame;
    end

    write_roll_outputs(pat, name, outDir, fps);
end

%% =========================================================
% Local function
%% =========================================================
function write_roll_outputs(pat, name, out_dir, fps)
% pat -> .mp4 (visual check) + .mat (simulation input).
    N = size(pat, 3);
    mp4_path = fullfile(out_dir, [name '.mp4']);
    mat_path = fullfile(out_dir, [name '.mat']);

    try
        vw = VideoWriter(mp4_path, 'MPEG-4');
        vw.Quality = 90;
    catch
        mp4_path = fullfile(out_dir, [name '.avi']);
        vw = VideoWriter(mp4_path, 'Motion JPEG AVI');
    end
    vw.FrameRate = fps;
    open(vw);
    for k = 1:N
        frame_u8 = uint8(round(pat(:, :, k) * 255));
        writeVideo(vw, flipud(frame_u8));   % el grows bottom->top, image rows top->bottom
    end
    close(vw);

    save(mat_path, 'pat', '-v7.3');
    fprintf('Saved: %s\n', mp4_path);
    fprintf('Saved: %s\n', mat_path);
end
