% script to try to write a trajectory using a pdb snapshot and the pdb2gro
% file I have written earlier
% ML 2013-10-29 

if(matlabpool('size')==0)
    matlabpool open
end
clear

addpath /data/Projects/grolab_matlab/grolab/io
addpath ./gromacs_files
addpath ./intermediate_files

%% SCRIPT PARAMETERS:

nlipids = [104 384 24];      buckling = 0.9;
%nlipids = [104 384 24];      buckling = 0.7;
%nlipids = [120 384 8];      buckling = 0.7;
%nlipids = [120 384 8];      buckling = 0.9;

numbeads = [13 13 25];

Heads = {2,2,[2 14]};    %indices for phosphate head groups (PG, PE, CL)
Tails = {[8, 13], [8,13], [8,13,20,25]}; %indices for innermost tail beads

lip = 3;    % 3 == collect x position of CL

%% filenames:
Nlip = sum(nlipids);

trrIN = sprintf('traj_lipids_%s_%s_%s.trr', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
CL_indices = sprintf('CL_indices_%s_%s_%s.mat', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
pdbIN = sprintf('lipids_edc_%s_%s_%s.pdb', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));

aligned_gro = sprintf('output_files/lipids_aligned_%s_%s_%s.gro', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
x_positions_mat = sprintf('output_files/x_CLpositions_%s_%s_%s.mat', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));

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
load(CL_indices, 'lwCLtails', 'upCLtails');

% a pdb snapshot from this simulation, to preserve molecule names etc
% need grolab in your path
M=pdbread(pdbIN); 
unix(['rm ' aligned_gro]); % remove previous result (since we append each frame)
align_fid=fopen(aligned_gro,'a');
%%
% figure(1)
% clf
% hold on
% axis equal
% view(-25,5)
% box on
% grid on
%
% datH=plot3(1,1,1,'r.');
% datT2=plot3(1,1,1,'k.');
% datM=plot3(1,1,1,'g.');
% axis([0 21 0 15 -8 8])


%% BUCKLESHAPE
xcomp0=buckling-0.05; xcomp1=buckling+0.05; % x-compression range for lookup table
xcomp_true=buckling;

%construct lookup tables
Lf=linspace(xcomp0,xcomp1,21); % range of compression factors in the table
XF=linspace(0,1,501);    % interpolation points
YF=zeros(length(Lf),length(XF));
clear B

% solve buckling problems numerically
parfor k=1:length(Lf)
    B(k)=buckleshape_half_bvp(Lf(k));
    % move to the interval 0<s<1
    B(k).x=B(k).x+0.5; 
    B(k).y(3,:)=B(k).y(3,:)+Lf(k)/2; % move to the right interval
    disp(int2str([k length(Lf)]))
end
% store results, rescaled to the interval 0 < x < 1
for k=1:length(Lf)
    % all lengths are scaled by 1/Lf to make 0 < x < 1
    x=B(k).y(3,:)/Lf(k);
    y=B(k).y(4,:)/Lf(k);
    
    YF(k,:)  =interp1(x,y,XF,'pchip');
    
%    disp(int2str([k length(Lf)]))
end
clear x y %dydx d2ydx2 kappa kappap xp yp xpp ypp

%% fit functions to individual frames
% buckled shape and derivatives

opt=optimset(@lsqcurvefit);
opt=optimset(opt,'Display','off','MaxFunEvals',1e5,'tolx',1e-9,'tolfun',1e-9);

%x=linspace(-Lx,2*Lx,1000);

FR=trrReadFrame(FS);
notdone=isfield(FR,'coord_XYZ');
pp_old=[mean(FR.coord_XYZ(:,3)) -5 5];
%% Do fitting:

x_CLup = [];
x_CLlw = [];

frame=0;

while(notdone)    
    
    Lx=FR.box_params(2);
    % prelim cos fit:
    fz0=@(p,x)(p(1)+p(3)*(1-cos(2*pi/Lx*(x-p(2)))));
    cosfit_p = lsqcurvefit(fz0, pp_old, FR.coord_XYZ(indT,1),FR.coord_XYZ(indT,3),[],[],opt);
    % buckleshape fit:
    FZ=@(p,x)(p(1)+interp2(XF,Lf,  YF,mod((x-p(2))/Lx,1),p(3),'*linear')*Lx);
    pp=lsqcurvefit(FZ,[cosfit_p(1:2) xcomp_true],FR.coord_XYZ(indT,1)',FR.coord_XYZ(indT,3)',[-inf -inf xcomp0],[inf inf xcomp1],opt);
        
    pp_old=pp;
    % align
    FR.coord_XYZ(:,1)=mod(FR.coord_XYZ(:,1)-pp(2),Lx); %pp(2)=0;
    FR.coord_XYZ(:,3)=FR.coord_XYZ(:,3)-pp(1);%%+20; pp(1)=20;
    
    for lcont=1:nlipids(lip)/2    %NOTE mean for 1 centre of mass tail bead
        x_CLup = [x_CLup; mean( FR.coord_XYZ(upCLtails{lcont},1) )];
        x_CLlw = [x_CLlw; mean( FR.coord_XYZ(lwCLtails{lcont},1) )];
    end
    
%    % plot present frame
%     set(datH,'xdata',FR.coord_XYZ(indH,1),...
%         'ydata',FR.coord_XYZ(indH,2),...
%         'zdata',FR.coord_XYZ(indH,3))
%     set(datT2,'xdata',FR.coord_XYZ(indT,1),...
%         'ydata',FR.coord_XYZ(indT,2),...
%         'zdata',FR.coord_XYZ(indT,3))
%     set(datM,'xdata',FR.coord_XYZ(indH,1),...
%         'ydata',FR.coord_XYZ(indH,2),...
%         'zdata',FZ([0 0 pp(3)],FR.coord_XYZ(indH,1)))
%     drawnow
%     pause(0.01)
    
    
    % write new coordinates to the pdb model, and then to a gro file
    for k=1:size(FR.coord_XYZ,1)
        M.Model.Atom(k).X=10*FR.coord_XYZ(k,1);
        M.Model.Atom(k).Y=10*FR.coord_XYZ(k,2);
        M.Model.Atom(k).Z=10*FR.coord_XYZ(k,3);
    end
    abc=10*FR.box_params([2 6 10]);
    M.Cryst1.a=abc(1);
    M.Cryst1.b=abc(2);
    M.Cryst1.c=abc(3);% may need to increase z box size (or better, make sure bilayer is centered)
    
    pdb2gro(aligned_gro,M,'aligned_test4_martini',FR.frame_time,'a',align_fid);
    % added align_fid so pdb2gro doesnt open and close the file each call
    
    % read next frame if there is one
    FR=trrReadFrame(FS);    
    notdone=isfield(FR,'coord_XYZ');
    
    frame=frame+1;
    if mod(frame,100)==0
        fprintf('%d ', frame)
    end
end

fclose(align_fid);
trrClose(FS);
%% histogram output

x_fin=linspace(-Lx,2*Lx,300);
buckle_fin=FZ(pp,x_fin);

save(x_positions_mat, 'x_CLup', 'x_CLlw', 'Lx', 'buckle_fin', 'x_fin', 'pp')

 
disp('convert to xtc and remove gro file from command line:')
%disp('echo 0 | trjconv -f test4.gro -s ../ML04_spherical_protein/flatpatch12by24_runA03.tpr -o test4.xtc -pbc mol')
disp('echo 0 | trjconv -f lipids_aligned.gro -s PWbilayer_PG120PE384CL8_5buckle0.9.tpr -o lipids_512_8_0.9.xtc -pbc mol')
disp('rm test4.gro')
