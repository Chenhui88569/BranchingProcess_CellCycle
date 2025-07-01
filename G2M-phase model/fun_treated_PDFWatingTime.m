function tilde_cdf  = fun_treated_PDFWatingTime( m_U, m_Dd , m_D, phase_pdf, U_pdf, Dd_pdf, D_pdf , num_t,t_span,  StepSize)

    mathscr_D = 1/m_D * myfun_ImproperIntegral_Cal(t_span,D_pdf);
    P_mathscr_D = myfun_ImproperIntegral_Cal(t_span,mathscr_D );
    mathscr_Dd  = 1/m_Dd *  myfun_ImproperIntegral_Cal(t_span,Dd_pdf); %pdf
    mathscr_Dd_cdf =  pdf_to_cdf(mathscr_Dd , StepSize);  
    Pd_mathscr_Dd =  1- mathscr_Dd_cdf;  %myfun_ImproperIntegral_Cal(t_span,mathscr_Dd   ); 
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
    for i =  2:num_t
        first_term_itoi1 = phase_pdf(i) * (1- tilde_p_i1mod3_i(i) );
        second_term_itoi1 = zeros(1,i);
        for j = 0:i-1
            second_term_itoi1(j+1)  = phase_pdf(j+1)*(1-tilde_p_d_i) * Pd_mathscr_Dd(i-j) * tilde_p_i1mod3_i(j+1) * mathscr_D(i-j);
        end
        % trapz(X,Y) integrates Y with respect to the coordinates or scalar spacing specified by X.
        second_term_itoi1_trapz = trapz(0:StepSize: t_span(i) ,second_term_itoi1);
        f_itoi1(i) =  first_term_itoi1  + second_term_itoi1_trapz ;
        %integral() is strictly for integrating a function (which might possibly be multi-valued), and is never for use for calculating the integral at a restricted list of points.
        % trapz() is for doing a trapezoid numeric integration for the given points. It returns a scalar.
    end
    tilde_cdf  = pdf_to_cdf(f_itoi1, StepSize);
   % tilde_cdf  =  tilde_cdf./ trapz(0:StepSize: t_span(end) ,f_itoi1);
    
end