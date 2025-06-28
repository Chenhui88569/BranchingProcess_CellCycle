function tilde_cdf  = fun_treated_PDFWatingTime_plot( )


para_unknown = [
    3
    20
    3
    50
    2

 ];
%%  ChoiceOfDist
ChoiceOfDist_untreated = 'Gamma';   ChoiceOfDist_treated_Dd = 'Weibull' ;
ChoiceOfDist_treated_D = 'Weibull';% 'Loglogistic';
ChoiceOfDist_treated_U = 'Exponential';
%% Time steps
StepSize= 0.1;
MaxT =  48; %unit in h
t_span = 0: StepSize: MaxT;
if t_span(end) ~=  MaxT
    t_span  = [t_span MaxT];
end
num_t = length(t_span);

Model_Index = 1;
%% scale idx to time zero
switch  Model_Index
    case 1
        idx = find(t_span == 250);  % Find where the trajectory converges
        Para_treated = [10.1563816036989;11.0674534469914;10.1470681258926;1.79695482015536];
        
    case 3
        idx = find(t_span == 258);  % Find where the trajectory converges
        Para_treated = [26.5129408241056;15.7185318119013;15.4319050483568;2.36165108138079;0.747086296925364;0.508772132761863;0.435713998475407];
end

a = Para_treated(1:3);
b =  Para_treated(4);
b_scale = 1/b;

%% G1
[G1_pdf, G1_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  t_span , a(1) ,b_scale);
%% S
[S_pdf, S_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  t_span , a(2) ,b_scale);
%% G2
[G2M_pdf, G2M_cdf]  =  PdfCdfCal(ChoiceOfDist_untreated,  t_span , a(3) ,b_scale);
phase_pdf  = G2M_pdf;
%%
[m_Dd, Dd_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_Dd,  t_span , para_unknown(2), para_unknown(1)  );
%[m_Dd, Dd_pdf, Dd_cdf]  =   PdfcdfCal_treated(ChoiceOfDist_treated_Dd,  t_span , para_unknown(2) , para_unknown(1));
Dd_pdf(Dd_pdf == inf) = 0;
%     % distribution of downtime along n_i -> n_j, log-logistic distribution
%     [m_D_G1, D_G1_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_D,  t_span , para_unknown(4)*DecresingFactor_D, para_unknown(3)*DecresingFactor_D );
%     D_G1_pdf( D_G1_pdf == inf) =  0;
%
%     [m_D_S, D_S_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_D,  t_span , para_unknown(6)*DecresingFactor_D, para_unknown(5)*DecresingFactor_D );
%     D_S_pdf( D_S_pdf == inf) =  0;

[m_D, D_pdf, ~]  =   PdfcdfCal_treated(ChoiceOfDist_treated_D,  t_span , para_unknown(4), para_unknown(3) );
D_pdf( D_pdf == inf) =  0;
%% distribution of uptime along n_i -> n_d &  n_i -> n_j

[m_U, U_pdf, U_cdf]  =   PdfcdfCal_treated(ChoiceOfDist_treated_U,  t_span , para_unknown(5),[] );
U_pdf( U_pdf == inf) = 0;


mathscr_D = 1/m_D * myfun_ImproperIntegral_Cal(t_span,D_pdf);
P_mathscr_D = myfun_ImproperIntegral_Cal(t_span,mathscr_D );
mathscr_Dd  = 1/m_Dd *  myfun_ImproperIntegral_Cal(t_span,Dd_pdf);
Pd_mathscr_Dd = myfun_ImproperIntegral_Cal(t_span,mathscr_Dd   );
tilde_p_i1mod3_i =  zeros(1,num_t  );
tilde_p_d_i =  m_U/(m_U+m_Dd );
%%
ImproperIntegral_U = myfun_ImproperIntegral_Cal(t_span,U_pdf );
for i =  2:num_t %tilde_p_i1_mod_3(i)
    temp_col = zeros(1,i);
    for j = 0:i-1   % integral 0 to s
        temp_col(j+1) = D_pdf(j+1) *ImproperIntegral_U(i-j);
    end
    tilde_p_i1mod3_i(i) = 1- trapz( t_span(1:i),    temp_col  );
end

% PDF of transition from n_j to n_{j+1} for three states
f_itoi1 = zeros(1,num_t);
tilde_cdf  = zeros(1,num_t);
for i =  2:num_t
    first_term_itoi1 = phase_pdf(i) * (1- tilde_p_i1mod3_i(i) );
    second_term_itoi1 = zeros(1,i);
    for j = 0:i-1
        second_term_itoi1(j+1)  = phase_pdf(j+1)*(1-tilde_p_d_i) * Pd_mathscr_Dd(i-j) * tilde_p_i1mod3_i(j+1) * mathscr_D(i-j);
    end
    % trapz(X,Y) integrates Y with respect to the coordinates or scalar spacing specified by X.
    second_term_itoi1_trapz = trapz(0:StepSize: t_span(i) ,second_term_itoi1);
    f_itoi1(i) =  first_term_itoi1  + second_term_itoi1_trapz ;
    tilde_cdf(i)   = trapz(0:StepSize: t_span(i) ,f_itoi1(1:i));
    %integral() is strictly for integrating a function (which might possibly be multi-valued), and is never for use for calculating the integral at a restricted list of points.
    % trapz() is for doing a trapezoid numeric integration for the given points. It returns a scalar.
end
tilde_cdf  =  tilde_cdf./ trapz(0:StepSize: t_span(end) ,f_itoi1);

figure
plot( t_span ,  f_itoi1./ trapz(0:StepSize: t_span(end) ,f_itoi1) );
hold on
plot(t_span , phase_pdf);
legend

figure
plot(t_span , tilde_cdf );

end