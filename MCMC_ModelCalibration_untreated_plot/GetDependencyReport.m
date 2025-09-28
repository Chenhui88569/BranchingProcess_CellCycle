close all
clear
fq = 'MCMC_main_untreated.m';
[fList,pList] = matlab.codetools.requiredFilesAndProducts(fq);
num_elem = length(fList);
newfList_withextension =  strings(1, num_elem);
newfList_withoutextension = strings(1, num_elem);
newdir = '/Users/chenhuima/MatlabWorkShop/5_FU_CellCycle_2022/MCMC_ModelCalibration_untreated';
for i = 1: num_elem
    temp = fList{i};
    newStr1 = erase(  temp,'/Users/chenhuima/MatlabWorkShop/5_FU_CellCycle_2022/MCMC_ModelCalibration_untreated/');
    newStr2 = erase(  newStr1,'.m');
    newfList_withoutextension(i) = newStr2;
    newfList_withextension(i) = newStr1 ;
    %copyfile( newStr1, fullfile(newdir, newStr1   ));
end
for i = 1: num_elem
    fprintf('%s ', newfList_withoutextension(i));
end

 fprintf('\n');
for i =  1:num_elem
    fprintf('%s ', newfList_withextension(i));
end
