% read data from the current frame of a trr file, and step forward to the
% next frame.
%
% The code is cut-n-pasted fromtrr2matlab.m by Evan Arthur at U Michigan,
% http://www.mathworks.com/matlabcentral/fileexchange/33312-convert-gromacs-v-4-5-trajectory-files-into-matlab-matrix 
% but divided into several chunks to allow frame-by-frame analysis instead
% of reading the whole trajectory at once.
%

function FR=trrReadFrame(FS,displayStats)

if(nargin==1)
    displayStats=false;
end

OPENfile=FS.OPENfile_read;
fileType=FS.fileType;

FR=struct;
FR.trajFile=FS.trajFile;
FR.fileType=FS.fileType;
FR.OPENfile_read=FS.OPENfile_read;

data_present_words = {'coordinate ', 'velocity ', 'force '};

if ~feof(OPENfile)
    
    FR.intro_words = fread(OPENfile, [1 51], '*char');            
    FR.data_present_xvf = fread(OPENfile, [1 3], 'int32');
    FR.num_atoms = fread(OPENfile, 1, 'int32');
    FR.frame_step = fread(OPENfile, [1 2], 'int32');
    FR.frame_time = fread(OPENfile, 1, fileType);
    FR.box_params = fread(OPENfile, [1 10], fileType);
    
    % display statistics in periodic intervals
    if displayStats
    fprintf(['\nopening file "%s",\n        ...%s\n' ...
         '    frame contains %g atoms,\n'...
         '    and is located at time step %u\n' ...
         '    frame is set at %g ps from beginning of simulation\n' ...
         '    box dimensions are %g by %g by %g nm\n' ...
         '  %sdata present in this frame\n' ], ...
         FR.trajFile, FR.intro_words(12:25), FR.num_atoms, FR.frame_step(1), ... 
         FR.frame_time, FR.box_params(2), FR.box_params(6), FR.box_params(10), ... 
         cell2mat(data_present_words(FR.data_present_xvf ~= 0)));
    end
     
    % coordinate data
    if (FR.data_present_xvf(1) ~= 0)        
        rawdata = fread(OPENfile, [1 (FR.num_atoms * 3)], fileType);        
        % structure data as XYZ * atoms * frames 3D matrix
        FR.coord_XYZ = permute(reshape(rawdata,3,FR.num_atoms,[]),[2 1 3]);
        %FR.xyzraw=rawdata;

    end    
    % velocity data
    if (FR.data_present_xvf(2) ~= 0)
        rawdata = fread(OPENfile, [1 (FR.num_atoms * 3)], fileType);
        FR.velocity_XYZ = permute(reshape(rawdata,3,FR.num_atoms,[]),[2 1 3]);
        %FR.vraw=rawdata;
    end    
    % force data
    if (FR.data_present_xvf(3) ~= 0)
        rawdata = fread(OPENfile, [1 (FR.num_atoms * 3)], fileType);
        FR.force_XYZ = permute(reshape(rawdata,3,FR.num_atoms,[]),[2 1 3]);
    end
    
    % fix offset for end of file
    FR.empty_spacing = fread(OPENfile, [1], '*char');
   
else
    disp('end of file rached!')
end