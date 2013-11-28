% open a trr_file file to figure out parameters for subsequent reading, and
% return a trrfile object FS for later use.
% Function to open a trr file and get ready to parse it frame by frame.
% The code is cut-n-pasted fromtrr2matlab.m by Evan Arthur at U Michigan,
% http://www.mathworks.com/matlabcentral/fileexchange/33312-convert-gromacs-v-4-5-trajectory-files-into-matlab-matrix 
% but divided into several chunks to allow frame-by-frame analysis instead
% of reading the whole trajectory at once.
%
% usage:
%
% FS=trrOpen_initRead(trrname); % open trr-file 
% 
%  FR1=trrReadFrame(FS); % read a frame and step forward
%  FR2=trrReadFrame(FS); % read the next frame
% ...
% trrClose(FS);
%
% ML 2013-10-27

function FS = trrOpen_initRead(trajFile,varargin)
% check if input exists
if exist(trajFile, 'file')
    fprintf('\n    Reading trajectory %s... \n', trajFile);
else
    fprintf('\n    Trajectory file not found, exiting program \n');
    return;
end

FS=struct;
FS.trajFile=trajFile;

%% determine single/double precision

% user override
if sum(strcmp(varargin, 'single')) > 0
    fileType = 'single';
end
if sum(strcmp(varargin, 'double')) > 0
    fileType = 'double';
end

% detect single/double
if sum(strcmp(varargin, 'single')) + sum(strcmp(varargin, 'double')) == 0
    
    OPENfile = fopen(trajFile, 'r', 'b');
    fseek(OPENfile, 0 ,'bof');
    precisionTest = fread(OPENfile, [1 9], 'int32');
    fclose(OPENfile);
    
    if precisionTest(9) == 36
        fprintf('single precision detected\n');
        fileType = 'single';
    elseif precisionTest(9) == 72
        fprintf('double precision detected\n');
        fileType = 'double';
    else
        fprintf('no precision dectected, defaulting to single precision\n');
        fileType = 'single';
    end    
end
FS.fileType=fileType;

% open input trr
OPENfile = fopen(trajFile, 'r', 'b');
fseek(OPENfile, 0 ,'bof'); % and point at the beginning

% fix offset for end of file
FS.empty_spacing = fread(OPENfile, [1], '*char');

% take a look at the first frame to get system size etc

if(~feof(OPENfile))
    frame_num = 1;
    data_present_words = {'coordinate ', 'velocity ', 'force '};

    intro_words = fread(OPENfile, [1 51], '*char');
    data_present = fread(OPENfile, [1 3], 'int32');
    num_atoms = fread(OPENfile, 1, 'int32');
    frame_step = fread(OPENfile, [1 2], 'int32');
    frame_time = fread(OPENfile, 1, fileType);
    box_params = fread(OPENfile, [1 10], fileType);
    
    % display statistics in periodic intervals
    fprintf(['\nopening file "%s",\n        ...%s\n' ...
         '    frame number %g contains %g atoms,\n'...
         '    and is located at time step %u\n' ...
         '    frame is set at %g ps from beginning of simulation\n' ...
         '    box dimensions are %g by %g by %g nm\n' ...
         '  %sdata present in this frame\n' ], ...
         trajFile, intro_words(12:25), frame_num, num_atoms, frame_step(1), ...
         frame_time, box_params(2), box_params(6), box_params(10), ... 
         cell2mat(data_present_words(data_present ~= 0)));
else
    fprintf(['No data present in this trr fi../ML04_spherical_protein/flatpatch12by24_runA03.trrle!'])
    FS=struct;
end
fclose(OPENfile);

% now start over to get ready for subsequent reading
% open input trr
FS.OPENfile_read = fopen(trajFile, 'r', 'b');
fseek(FS.OPENfile_read, 0 ,'bof'); % and point at the beginning

% fix offset for end of file
empty_spacing = fread(FS.OPENfile_read, [1], '*char');




