% script to build bond matrix M, and then write the psf file on order to show bonds in VMD
clear

nlipids = [104 384 24];      buckling = 0.9;
%nlipids = [104 384 24];  buckling=0.7;
%nlipids = [120 384 8];  
numbeads = [13 13 25];

PGbonds = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 3 9; 9 10; 10 11; 11 12; 12 13];
PEbonds = PGbonds;
CLbonds = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 3 9; 9 10; 10 11; 11 12; 12 13; 1 14; 14 15; 15 16; 16 17; 17 18; 18 19; 19 20; 15 21; 21 22; 22 23; 23 24; 24 25 ];

bonds = {PGbonds, PEbonds, CLbonds};
bondslen = [length(PGbonds) length(PEbonds) length(CLbonds)];
BB = zeros(nlipids * bondslen', 2);

for lip = 1:length(nlipids)
    if lip > 1
        lip_ofst = nlipids(1:lip-1) * numbeads(1:lip-1)';
        bnd_ofst = nlipids(1:lip-1) * bondslen(1:lip-1)';
    else
        lip_ofst = 0;
        bnd_ofst = 0;
    end
    
    for ni=1:nlipids(lip)
        BB(bnd_ofst+(ni-1)*bondslen(lip)+1:bnd_ofst+ni*bondslen(lip),:) =  ...
            bonds{lip} + (ni-1)*numbeads(lip) + lip_ofst;
    end
    
end
%% now just call psfwrite:

% first produce the pdb from gro file with trjconv or editconf:
%       editconf -f foo.gro -o foo.pdb
%
% NOTE: I changed the name POPG to GPOP in the gro file (using sed):
%       sed -i 's/POPG/GPOP/g' noWater_PWbilayer_PG120PE384CL8_5buckle0.9_GPOP.gro 
%
%       because pdb only uses 3 letters for residue name 
%       (and thus POPG and POPE were both POP in pdb format)

lip=3; Nlip=sum(nlipids); 

pdbIN = sprintf('lipids_edc_%s_%s_%s.pdb', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));
psfOUT = sprintf('output_files/lipids_edc_%s_%s_%s.psf', num2str(Nlip), num2str(nlipids(lip)), num2str(buckling));

M=pdbread(pdbIN); 
psfwrite(psfOUT, M, BB, 'PG120/PE384/CL8__buckle_0.7');
