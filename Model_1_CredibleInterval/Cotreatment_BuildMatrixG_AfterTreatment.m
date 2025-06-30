function   mat_G = Cotreatment_BuildMatrixG_AfterTreatment(obj,G1_cdf,G2M_cdf,S_cdf,TotalState)

lambda_d_max_G1_G2M = obj.para_unknown(10);
lambda_G2arrest_UR_max = obj.para_unknown(11);
lambda_G2arrest_FR_max = obj.para_unknown(12);
lambda_G2arrest_MS_max = obj.para_unknown(13);
lambda_G1arrest_max = obj.para_unknown(14);
lambda_d_max_S_G2M = obj.para_unknown(15);
lambda_Sarrest_max = obj.para_unknown(16);

tilde_b_G2M =  obj.para_unknown(17);
tilde_b_d_G2M =  obj.para_unknown(18);
n_G2M =  obj.para_unknown(19);
EC50_d_G2M =  obj.para_unknown(20);
EC50_G2M =  obj.para_unknown(21);

Sphase_StartIdx = 21;
lambda_d_max_G2M = obj.para_unknown(Sphase_StartIdx + 8);
lambda_Sarrest_UR_max = obj.para_unknown(Sphase_StartIdx + 9);
lambda_Sarrest_FR_max = obj.para_unknown(Sphase_StartIdx + 10);
lambda_G2Marrest_max = obj.para_unknown(Sphase_StartIdx + 11);
lambda_G1block_max = obj.para_unknown(Sphase_StartIdx + 12);
tilde_b_S =  obj.para_unknown(Sphase_StartIdx +13);
tilde_b_d_S =  obj.para_unknown(Sphase_StartIdx + 14);
n_S =  obj.para_unknown(Sphase_StartIdx + 15);
EC50_d_S =  obj.para_unknown(Sphase_StartIdx + 16);
EC50_S =  obj.para_unknown(Sphase_StartIdx + 17);

Dose_G2M = obj.Dose_1thDrug;
Dose_S = obj.Dose_2thDrug;
%% G2M phase drug
% create G2 arrest
lambda_fun = @(lambda_max) lambda_max * ((obj.Dose_G2M / EC50_G2M)^n_G2M) / (1 + (obj.Dose_G2M/ EC50_G2M)^n_G2M);

lambda_G2Marrest_UR  = lambda_fun(lambda_G2arrest_UR_max);
[~, ~, G2_arrest_UR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_UR  , tilde_b_G2M  ); %A: scale, B , shape
lambda_G2Marrest_FR  = lambda_fun(lambda_G2arrest_FR_max);
[~, ~, G2_arrest_FR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_FR  , tilde_b_G2M  ); %A: scale, B , shape
lambda_G2Marrest_MS = lambda_fun(lambda_G2arrest_MS_max);
[~, ~, G2_arrest_MS ]  =   PdfcdfCal_treated('Weibull',  obj.t_span , lambda_G2Marrest_MS  , tilde_b_G2M  ); %A: scale, B , shape
% Create tilde_G1
lambda_G1arrest = lambda_fun(lambda_G1arrest_max);
lambda_d_max_G1_G2M =  lambda_fun(lambda_d_max_G1_G2M);
[ m_D_G1arrest ,D_G1arrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_G1arrest  , tilde_b_G2M  );
D_G1arrest_pdf( D_G1arrest_pdf == inf) = 0;
[ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_max_G1_G2M  , tilde_b_d_G2M  );
Dd_pdf( Dd_pdf == inf) = 0;
ChoiceOfDist_treated_U = 'Exponential';
[m_U, U_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U, obj.t_span , 1 ,[] );
U_pdf( U_pdf == inf) = 0;
tilde_G1_G2M  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_G1arrest, G1_pdf , U_pdf, Dd_pdf, D_G1arrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize   )  ;

% Create tilde_S
lambda_Sarrest = lambda_fun(lambda_Sarrest_max);
lambda_d_max_S_G2M =  lambda_fun(lambda_d_max_S_G2M);
[ m_D_Sarrest ,D_Sarrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_Sarrest  , tilde_b_G2M  );
D_Sarrest_pdf( D_Sarrest_pdf == inf) = 0;
[ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span ,  lambda_d_max_S_G2M  , tilde_b_d_G2M  );
Dd_pdf( Dd_pdf == inf) = 0;


tilde_S_G2M = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_Sarrest, S_pdf , U_pdf, Dd_pdf, D_Sarrest_pdf ,  obj.num_t, obj.t_span , obj.StepSize);

%% S phase drug
% create S arrest
lambda_fun = @(lambda_max) lambda_max * ((obj.Dose_S  / EC50_S)^n_S) / (1 + (obj.Dose_S  / EC50_S)^n_S);

% create G1 block
lambda_G1block  = lambda_fun(lambda_G1block_max );
[~, ~, G1block ]  =   PdfcdfCal_treated('Weibull',  obj.t_span_AfterTreatment ,  lambda_G1block  , tilde_b_S  ); %A: scale, B , shape

% create S arrest
lambda_Sarrest_UR  = lambda_fun(lambda_Sarrest_UR_max);
[~, ~, S_arrest_UR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span_AfterTreatment , lambda_Sarrest_UR  , tilde_b_S  ); %A: scale, B , shape
lambda_Sarrest_FR  = lambda_fun(lambda_Sarrest_FR_max);
[~, ~, S_arrest_FR ]  =   PdfcdfCal_treated('Weibull',  obj.t_span_AfterTreatment , lambda_Sarrest_FR  , tilde_b_S  ); %A: scale, B , shape
% Create tilde_G2M
lambda_G2Marrest = lambda_fun(lambda_G2Marrest_max);
lambda_d_G2M =  lambda_fun(lambda_d_max_G2M);
[ m_D_G2Marrest ,D_G2Marrest_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span_AfterTreatment ,   lambda_G2Marrest   , tilde_b_S  );
D_G2Marrest_pdf(D_G2Marrest_pdf  == inf) = 0;
[ m_Dd , Dd_pdf , ~ ]  =   PdfcdfCal_treated('Weibull',  obj.t_span_AfterTreatment ,  lambda_d_G2M  , tilde_b_d_S  );
Dd_pdf( Dd_pdf == inf) = 0;
ChoiceOfDist_treated_U = 'Exponential';
[m_U, U_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U, obj.t_span, 1 ,[] );
U_pdf( U_pdf == inf) = 0;
tilde_G2M_S  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D_G2Marrest, G1_pdf , U_pdf, Dd_pdf, D_G2Marrest_pdf ,  obj.num_t, obj.t_span, obj.StepSize   )  ;

mat_G = zeros(TotalState, TotalState , obj.num_t);

mat_G(1,1,:) = G1_cdf;
mat_G(2,2,:) = G1block;
mat_G(3,3,:) = S_cdf;
mat_G(4,4,:) =S_arrest_UR ;
mat_G(5,5,:) =S_arrest_FR ;
mat_G(6,6,:) =G2M_cdf;
mat_G(7,7,:) =G2_arrest_UR  ;
mat_G(8,8,:) =G2_arrest_FR ;
mat_G(9,9,:) =G2_arrest_MS ;
mat_G(10,10,:) = G1_cdf;
mat_G(11,11,:) = G1block;
mat_G(12,12,:) = tilde_G1_G2M;
mat_G(13,13,:) = S_cdf;
mat_G(14,14,:) =S_arrest_UR ;
mat_G(15,15,:) =S_arrest_FR ;
mat_G(16,16,:) =G2M_cdf;
mat_G(17,17,:) =G2_arrest_UR  ;
mat_G(18,18,:) =G2_arrest_FR ;
mat_G(19,19,:) =G2_arrest_MS ;

for i = 1: obj.num_RroundsOfTreatedcells-1
    startnum = 20 +( i-1 )*19;
    endnum =  20 + i*19 -1;
    mat_G(startnum   : endnum, startnum :  endnum,:    ) = mat_G( 1 :19, 1:19,:);
end

True_num_GenOfHealthycells = obj.num_GenOfHealthycells -1 ;
if True_num_GenOfHealthycells > 0
    for i = 1:True_num_GenOfHealthycells
        startnum = 19 *num_RroundsOfTreatedcells +1 + (i-1)*3;
        mat_G(startnum , startnum,: ) = G1_cdf;
        mat_G(startnum+1, startnum+1,: ) = S_cdf;
        mat_G(startnum+2, startnum+2,: ) = G2M_cdf;
    end
end
num_t_AfterTreatment = length(obj.t_span_AfterTreatment );
mat_G = mat_G(:,:,1:num_t_AfterTreatment );
[~, ~, G_apop_cdf]  =   PdfcdfCal_treated('Exponential',  t_span_AfterTreatment  ,5 ,  [] );
mat_G(end,  end,:   )  =    G_apop_cdf; % time for eliminating the apoptotic cell
end