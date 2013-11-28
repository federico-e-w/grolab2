% script to verify convergence of fourier components
% perform a check every n frames
% (possible flip-flops?)

if(matlabpool('size')==0)
    matlabpool open
end
clear

addpath /data/Projects/grolab_matlab/grolab/io
addpath ./gromacs_files
addpath ./intermediate_files

%% SCRIPT PARAMETERS:

frame_dt=10;
end_frame=20000;

%nlipids = [104 384 24];      buckling = 0.9;
%nlipids = [104 384 24];      buckling = 0.7;
%nlipids = [120 384 8];      buckling = 0.7;
nlipids = [120 384 8];      buckling = 0.9;

numbeads = [13 13 25];

Heads = {2,2,1};    %indices for phosphate head groups (PG, PE, CL) NOTE CL head bead NOT phosphate; phosphates: [2 14]
Tails = {[8, 13], [8,13], [8,13,20,25]}; %indices for innermost tail beads

lip = 3;    % 3 == collect x position of CL

%% filenames:
Nlip = sum(nlipids);

trrIN = sprintf('traj_lipids_%s_%s_%s.trr', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
CL_indices = sprintf('CL_indices_%s_%s_%s.mat', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));

%% calculate indices of tail and head beads across the bilayer
% NOTE this assumes that the order in the lipids in gromacs files is the same as in the
% variables in PARAMETERS section, and that there are no other molecules before those

indH = Heads{1}:numbeads(1):nlipids(1)*numbeads(1);
indT = [Tails{1}(1):numbeads(1):nlipids(1)*numbeads(1) Tails{1}(2):numbeads(1):nlipids(1)*numbeads(1)];

for lip = 2:length(nlipids)
    for i=1:length(Heads{lip})
        indH = [indH, (sum(nlipids(1:lip-1).*numbeads(1:lip-1))+Heads{lip}(i)) : ...
            numbeads(lip) : sum(nlipids(1:lip).*numbeads(1:lip)) ] ; %#ok<AGROW>
    end
    for i=1:length(Tails{lip})
        indT = [indT, (sum(nlipids(1:lip-1).*numbeads(1:lip-1))+Tails{lip}(i)) : ...
            numbeads(lip) : sum(nlipids(1:lip).*numbeads(1:lip)) ] ; %#ok<AGROW>
    end
    
end
%% load files


FS=trrOpen_initRead(trrIN);
load(CL_indices, 'lwCLtails', 'upCLtails', 'lwCLheads', 'upCLheads');

% a pdb snapshot from this simulation, to preserve molecule names etc
% need grolab in your path


FR=trrReadFrame(FS);
notdone=isfield(FR,'coord_XYZ');
%pp_old=[mean(FR.coord_XYZ(:,3)) -5 5];
%% Do fitting:

x_CL_ind = [];
% for nl=1:length(upCLtails)
%     x_CL_ind = [x_CL_ind; upCLtails{nl}(:); lwCLtails{nl}(:) ];
% end

frame=-1;

Fourier_coeff=[];
Fourier_coeff_CL=[];
tails_x_CL=[];
time = [];
while(notdone && frame <= end_frame)
    
    frame=frame+1;
    if mod(frame, frame_dt)~=0
        FR=trrReadFrame(FS);
        notdone=isfield(FR,'coord_XYZ');
        continue
    end
    
    if mod(frame,1000)==0
        fprintf(' %d ', frame);
    end
    
    Lx=FR.box_params(2);
    % int(\delta(X)e^(-i n pi x/L),x)
    F_nx=@(nx,X)(1/Lx*sum(exp(1i*2*pi*X*nx/Lx)));
    
    tails_x = FR.coord_XYZ(indT,1);
    %    tails_x_CL = FR.coord_XYZ(x_CL_ind,1);
    for nl=1:length(upCLheads); %NOTE need to be split in 2 cases if different #CL in leaflets
        tails_x_CL = [tails_x_CL; mean(FR.coord_XYZ(upCLtails{nl},1)); mean(FR.coord_XYZ(lwCLtails{nl},1))];
    end
    
    Fourier_coeff = [Fourier_coeff; F_nx(0,tails_x) F_nx(1,tails_x) F_nx(-1,tails_x) F_nx(2,tails_x) F_nx(-2,tails_x) F_nx(3,tails_x) F_nx(-3,tails_x) F_nx(4,tails_x) F_nx(-4,tails_x) F_nx(6,tails_x) F_nx(-6,tails_x) F_nx(8,tails_x) F_nx(-8,tails_x)];
    Fourier_coeff_CL = [Fourier_coeff_CL; F_nx(0,tails_x_CL) F_nx(1,tails_x_CL) F_nx(-1,tails_x_CL) F_nx(2,tails_x_CL) F_nx(-2,tails_x_CL) F_nx(3,tails_x_CL) F_nx(-3,tails_x_CL) F_nx(4,tails_x_CL) F_nx(-4,tails_x_CL) F_nx(6,tails_x_CL) F_nx(-6,tails_x_CL) F_nx(8,tails_x_CL) F_nx(-8,tails_x_CL)];
    % read next frame if there is one
    FR=trrReadFrame(FS);
    notdone=isfield(FR,'coord_XYZ');
    
    time = [time; frame];
    
end

%fclose(align_fid);
trrClose(FS);
%%

figure(437)
clf
plot(time, sqrt(Fourier_coeff(:,2) .* Fourier_coeff(:,3)), 'b'), hold on
plot(time, sqrt(Fourier_coeff(:,4) .* Fourier_coeff(:,5)), 'r')
plot(time, sqrt(Fourier_coeff(:,6) .* Fourier_coeff(:,7)), 'k')
plot(time, sqrt(Fourier_coeff(:,8) .* Fourier_coeff(:,9)), 'm')
plot(time, sqrt(Fourier_coeff(:,12) .* Fourier_coeff(:,13)), 'g')
xlabel('time [ns]')
ylabel('Fourier coeffcients')
legend('c_1', 'c_2', 'c_3', 'c_4', 'c_8')

figure(438)
clf
plot(time, sqrt(Fourier_coeff_CL(:,2) .* Fourier_coeff_CL(:,3)), 'b'), hold on
plot(time, sqrt(Fourier_coeff_CL(:,4) .* Fourier_coeff_CL(:,5)), 'r')
plot(time, sqrt(Fourier_coeff_CL(:,6) .* Fourier_coeff_CL(:,7)), 'k')
plot(time, sqrt(Fourier_coeff_CL(:,8) .* Fourier_coeff_CL(:,9)), 'm')
plot(time, sqrt(Fourier_coeff_CL(:,12) .* Fourier_coeff_CL(:,13)), 'g')
xlabel('time [ns]')
ylabel('Fourier coeffcients_CL')
legend('c_1', 'c_2', 'c_3', 'c_4', 'c_8')



