function [y_pdf,y_cdf]  = PdfCdfCal(ChoiceOfDist,  t_span, a ,b)
if strcmp(ChoiceOfDist,'Gamma' )
    pd = makedist(ChoiceOfDist,'a',a,'b',b);%distribution object
   
elseif  strcmp(ChoiceOfDist,'Weibull' ) 
   pd = makedist(ChoiceOfDist,'A',a,'B',b);%distribution object
elseif  strcmp(ChoiceOfDist,'Exponential' ) 
   pd = makedist(ChoiceOfDist,'mu',a);%distribution object
end
y_cdf = cdf(pd,  t_span);% return vec
y_pdf = pdf(pd, t_span);% return vec\\\
end