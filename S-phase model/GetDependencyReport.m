close all
clear
fq = 'S_phase_MCMC_main_treated.m';
[fList,pList] = matlab.codetools.requiredFilesAndProducts(fq);
num_elem = length(fList);
newfList_withextension =  strings(1, num_elem);
newfList_withoutextension = strings(1, num_elem);
newdir = '/Users/chenhuima/MatlabWorkShop/5_FU_CellCycle_2022/MCMC_ModelCalibration_treated_model_1_dtx_GM3_new/Cotreatment_compute';
mkdir(newdir  );
for i = 1: num_elem
    temp = fList{i};
    newStr1 = erase(  temp,'/Users/chenhuima/Matlab_code/BranchingProcess_CellCycle/MCMC_ModelCalibration_treated_model_1_dtx_GM3_new/BranchingProcess_CellCycle/S-phase model/');
    newStr2 = erase(  newStr1,'.m');
    newfList_withoutextension(i) = newStr2;
    newfList_withextension(i) = newStr1 ;
    try
        %copyfile( newStr1, fullfile(newdir, newStr1   ));
    catch e
        continue;
    end
end
for i = 1: num_elem
    fprintf('%s\n', newfList_withoutextension(i));
end

 fprintf('\n');
for i =  1:num_elem
    fprintf('%s \n', newfList_withextension(i));
end
