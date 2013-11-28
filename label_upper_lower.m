% script to get the indices of lower and upper CL

clear

addpath /data/Projects/grolab_matlab/grolab/io
addpath ./gromacs_files

%% SCRIPT PARAMETERS:

%nlipids = [120 384 8];      buckling = 0.9;
%nlipids = [120 384 8];      buckling = 0.7;
nlipids = [104 384 24];      buckling = 0.9;
%nlipids = [104 384 24];      buckling = 0.7;
numbeads = [13 13 25];  Nlip = sum(nlipids); lip=3;

Heads = {2,2,1};%{2,2,[2 14]};    %indices for phosphate head groups (PG, PE, CL) NOTE maybe use the GL0 ind = 1 for CL
Tails = {[8, 13], [8,13], [8,13,20,25]}; %indices for innermost tail beads

trr_flatIN = sprintf('traj_lipids_flat_%s_%s_%s.trr', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
FS=trrOpen_initRead(trr_flatIN);

%% calculate indices of tail and head beads across the bilayer
%%
lip = 3;    %indices of CL tails (4 per CL)

CLheads = (sum(nlipids(1:lip-1).*numbeads(1:lip-1)) + Heads{lip}(1) ) : ...  %the 1 is because the head bead is bead 1
    numbeads(lip) : sum(nlipids(1:lip).*numbeads(1:lip)) ;

CLheads = num2cell(CLheads(:));    %cell array, first column bead indices, 2nd column 'u' or 'l'     %mat2cell(CLheads,1,ones(size(CLheads)));

for i=1:nlipids(lip)
    CLtails{i,1} = Tails{lip} + sum(nlipids(1:lip-1).*numbeads(1:lip-1)) + (i-1)*numbeads(lip) ;
end

numCLtailbeads = length(Tails{lip});

%%
% a pdb snapshot from this simulation, to preserve molecule names etc
% need grolab in your path
%M=pdbread('lipids_edc.pdb');
%!rm lipids_aligned.gro


FR=trrReadFrame(FS);
notdone=isfield(FR,'coord_XYZ');
assert(notdone);
%% Do fitting:
upCLheads = {};
upCLtails = {};
lwCLheads = {};
lwCLtails = {};
for nl=1:nlipids(lip)
    
    if FR.coord_XYZ(CLheads{nl},3) > max(FR.coord_XYZ(CLtails{nl},3))
        upCLheads = [upCLheads; CLheads{nl} ];
        upCLtails = [upCLtails; CLtails{nl}];
    elseif FR.coord_XYZ(CLheads{nl},3) < min(FR.coord_XYZ(CLtails{nl},3))
        lwCLheads = [lwCLheads; CLheads{nl}];
        lwCLtails = [lwCLtails; CLtails{nl}];
    end
  
end

assert(length(upCLheads) == length(upCLtails) && length(lwCLheads) == length(lwCLtails) && length(upCLheads) == length(lwCLtails) );

trrClose(FS);
outCL = sprintf('intermediate_files/CL_indices_%s_%s_%s.mat', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));

save(outCL, 'upCLheads', 'upCLtails', 'lwCLheads', 'lwCLtails');
