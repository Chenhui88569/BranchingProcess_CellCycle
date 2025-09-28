function [s, m]  = PdfcdfCal_MeanStd(ChoiceOfDist, a ,b)
if strcmp(ChoiceOfDist,'Gamma' )
    pd = makedist(ChoiceOfDist,'a',a,'b',b);%distribution object
   
elseif  strcmp(ChoiceOfDist,'Weibull' ) 
   pd = makedist(ChoiceOfDist,'A',a,'B',b);%distribution object
   
elseif  strcmp(ChoiceOfDist,'Loglogistic' ) 
   pd = makedist(ChoiceOfDist,'mu',a,'sigma',b);%distribution object
   
elseif  strcmp(ChoiceOfDist,'Exponential' ) 
   pd = makedist(ChoiceOfDist,'mu',a);%distribution object
end
m = mean(pd); % the mean m of the probability distribution pd.
s = std(pd);
end