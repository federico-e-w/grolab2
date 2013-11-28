% script to verify that CLup and CLlw are correct:
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

check_frames=100;

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


%% BUCKLESHAPE
xcomp0=buckling-0.05; xcomp1=buckling+0.05; % x-compression range for lookup table
xcomp_true=buckling;

%construct lookup tables
Lf=linspace(xcomp0,xcomp1,21); % range of compression factors in the table
XF=linspace(0,1,501);    % interpolation points
YF=zeros(length(Lf),length(XF));
dYF=zeros(length(Lf),length(XF));
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
    dydx=B(k).yp(4,:)./B(k).yp(3,:);
    
    YF(k,:)  =interp1(x,y,XF,'pchip');
    dYF(k,:) =interp1(x,dydx,XF,'pchip');
    
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
    
    frame=frame+1;
    if mod(frame,check_frames)~=0
        FR=trrReadFrame(FS);    
        notdone=isfield(FR,'coord_XYZ');
        continue
    end
    
    if mod(frame,10*check_frames)==0
        disp(frame)
    end
    
    Lx=FR.box_params(2);
    % prelim cos fit:
    fz0=@(p,x)(p(1)+p(3)*(1-cos(2*pi/Lx*(x-p(2)))));
    cosfit_p = lsqcurvefit(fz0, pp_old, FR.coord_XYZ(indT,1),FR.coord_XYZ(indT,3),[],[],opt);
    % buckleshape fit:
    FZ=@(p,x)(p(1)+interp2(XF,Lf,  YF,mod((x-p(2))/Lx,1),p(3),'*linear')*Lx);
    pp=lsqcurvefit(FZ,[cosfit_p(1:2) xcomp_true],FR.coord_XYZ(indT,1)',FR.coord_XYZ(indT,3)',[-inf -inf xcomp0],[inf inf xcomp1],opt);
    
    dFZ =@(p,x)(p(1)+interp2(XF,Lf, dYF,mod((x-p(2))/Lx,1),p(3),'*linear'));
    
    
    pp_old=pp;
    % align
    FR.coord_XYZ(:,1)=mod(FR.coord_XYZ(:,1)-pp(2),Lx); %pp(2)=0;
    FR.coord_XYZ(:,3)=FR.coord_XYZ(:,3)-pp(1);%%+20; pp(1)=20;
    
    %check correct labeling. compare head_z > tail_z if slope < 45 deg. otherwise check x coordinates
    for nl=1:length(upCLheads)
        x_H_up = FR.coord_XYZ(upCLheads{nl},1);
        x_T_up = mean(FR.coord_XYZ(upCLtails{nl},1));
        z_H_up = FR.coord_XYZ(upCLheads{nl},3);
        z_T_up = mean(FR.coord_XYZ(upCLtails{nl},3));
        
        slope_T_up = dFZ([0 0 pp(3)], x_T_up);
        if abs(slope_T_up) >=0 && abs(slope_T_up) <=1
            if z_H_up < z_T_up
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', upCLheads{nl}, frame); 
            end
        elseif slope_T_up > 1
            if x_H_up > x_T_up
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', upCLheads{nl}, frame); 
            end
        elseif slope_T_up < -1
            if x_H_up < x_T_up
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', upCLheads{nl}, frame); 
            end
        else
             %disp('ERROR 84')
             error('Error 84')
        end
     
    end
    for nl=1:length(lwCLheads)
        x_H_lw = FR.coord_XYZ(lwCLheads{nl},1);
        x_T_lw = mean(FR.coord_XYZ(lwCLtails{nl},1));
        z_H_lw = FR.coord_XYZ(lwCLheads{nl},3);
        z_T_lw = mean(FR.coord_XYZ(lwCLtails{nl},3));
        
        slope_T_lw = dFZ([0 0 pp(3)], x_T_lw);
        if abs(slope_T_lw) >=0 && abs(slope_T_lw) <=1
            if z_H_lw > z_T_lw
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', lwCLheads{nl}, frame); 
            end
        elseif slope_T_lw > 1
            if x_H_lw < x_T_lw
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', lwCLheads{nl}, frame); 
            end
        elseif slope_T_lw < -1
            if x_H_lw > x_T_lw
                fprintf('possible flip-flop in molecule with head bead %d, in frame %d \n', lwCLheads{nl}, frame); 
            end
        else
             %disp('ERROR 84')
             error('Error 84')
        end
     
    end
    % read next frame if there is one
    FR=trrReadFrame(FS);    
    notdone=isfield(FR,'coord_XYZ');
    
    
    
end

%fclose(align_fid);
trrClose(FS);

