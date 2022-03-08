%% Input and Output files
clear;clc;
cd /n02dat01/users/ypwang/AHBA/Github_20220301/Scripts/Step04_Behavior/Reference/
od =  '/n02dat01/users/ypwang/AHBA/Github_20220301/Scripts/Step04_Behavior/Results/';
SubID = textread( '/n02dat01/users/ypwang/AHBA/Github_20220301/Scripts/Step04_Behavior/Reference/ToPCA_211.txt','%s' );

%% Read in FC of ROI 1:17 from FC matrix
cere2cere_allsub = zeros( length( SubID ), 136 ); % 211*136
for i = 1:length( SubID )
    sub = SubID{i};
    brain2brain_path = strcat( '/n02dat01/users/lchai/HCP-FISHER-FC/GSR-cere17-cor17-all/', sub, '.csv');  % The individual hcp fc data, change to your path
    brain2brain_fid = importdata( brain2brain_path );
    brain2brain = brain2brain_fid;
    cere2cere = brain2brain( 115:131, 115:131 ); % 17 * 17
    mask_cere = tril(true(size(cere2cere)),-1);
    cere2cere_noself = cere2cere(mask_cere);
    cere2cere_allsub( i, : ) = transpose(cere2cere_noself);
end
csvwrite( strcat(od, 'FC136_n211_cere.csv' ), cere2cere_allsub);

%% PALM
cd /n02dat01/users/ypwang/AHBA/Github_20220301/Scripts/Step04_Behavior/Results/
wb_command='/n21dat01/kxli/Tools/workbench/bin_rh_linux64/wb_command';

% 211 subjects, 59 behavior
palm -i FC136_n211_cere.csv -d PCAscore_PC1.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC1
palm -i FC136_n211_cere.csv -d PCAscore_PC2.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC2
palm -i FC136_n211_cere.csv -d PCAscore_PC3.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC3
palm -i FC136_n211_cere.csv -d PCAscore_PC4.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC4
palm -i FC136_n211_cere.csv -d PCAscore_PC5.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC5
palm -i FC136_n211_cere.csv -d PCAscore_PC6.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC6
palm -i FC136_n211_cere.csv -d PCAscore_PC7.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC7
palm -i FC136_n211_cere.csv -d PCAscore_PC8.csv -t FCPCA_n211_t1.csv -twotail -zstat -fdr -accel tail -n 10000 -demean -o FC_n211_cere_PC8

%% Read the z score into 17*17 matrix
load_matrix1=zeros(size(cere2cere));
bw=true(size(cere2cere));

for componant = 1:8 
    % z score
    input_path = strcat( od, 'FC_n211_cere_PC', num2str( componant  ),'_dat_ztstat.csv' );
    input = importdata( input_path );
    Z_FC_PC = zeros(size(cere2cere));
    bw=true(size(cere2cere));
    Z_FC_PC(tril(bw, -1))=input ;
    Z_FC_PC =Z_FC_PC + Z_FC_PC';
    Z_FC_PC( Z_FC_PC==0) =1;
    csvwrite( strcat( od,  'FC_n211_cere_PC', num2str( componant  ),'_matrix_z.csv' ), Z_FC_PC );
    
    % permutation p value
    P_FC_PC = zeros(size(cere2cere));
    input_path = strcat( od, 'FC_n211_cere_PC', num2str( componant  ),'_dat_ztstat_uncp.csv' );
    input = importdata( input_path );
    P_FC_PC = zeros(size(cere2cere));
    bw=true(size(cere2cere));
    P_FC_PC(tril(bw, -1))=input ;
    P_FC_PC = P_FC_PC + P_FC_PC';
    P_FC_PC( P_FC_PC==0) =1;
    csvwrite( strcat( od,  'FC_n211_cere_PC', num2str( componant  ),'_matrix_z_uncp.csv' ), P_FC_PC );
    
end
