function [mat_Gm_tilde_mixture,I_mat_G_tilde_mixture]  = myfun_Gm_tilde_mixture_modified(mat_G, mat_Gm, mat_m, I_mat_G ,...
    a, b_scale,para_unknown, t_span, StepSize, States_convergence, TreatmentEffectOnPhases, leftout_CyclingG2M)

%% \tildeG matrix stores of waiting time prolonged by the checkpoint activation
num_t = length(t_span);
d_G1 = para_unknown(11) ;  d_S = para_unknown(12) ;  d_G2M = para_unknown(13);
m12 = 1;
m1d = 0;
m23 = 1;
m2d = 0;
m31 = 2;
m3d = 0;

%para_unknown = para_unknown(:);
%mixture_para = num2cell( [10,30 , 0.5, 0.6, 1, 3,0.5,0.5  ]) ;
mixture_para = num2cell(  para_unknown( 18 : 29) );
[factor_G2M_1, factor_G2M_2, p2_G2M,  p3_G2M, factor_G1_1, factor_S_1, p1_G1, p1_S, factor_G1_2, factor_S_2, p2_G1, p2_S ] = mixture_para{:} ;
factor_rest = [factor_G1_1, factor_G1_2
    factor_S_1,  factor_S_2  ] ;
p_rest = [p1_G1, p2_G1; 
      p1_S, p2_S];
%%
AbsolutePhase = mod(States_convergence,3);
mat_Gm_tilde_mixture = mat_Gm;
I_mat_G_tilde_mixture =  I_mat_G ;
m_col = [ m12*(1-d_G1)
    m23*(1-d_S)
    ];
d_col = [ d_G1 ;  d_S];
switch TreatmentEffectOnPhases
    case 0 % radiation
        
    case 2
    case 3 %G2M
        if  ~isempty(find(AbsolutePhase == 0, 1)) % Cycling G2M phase starts from G2M
            idx_CyclingG2M = find(AbsolutePhase == 0);
            %if length(idx_CyclingG2M) == 1
                for  i = 1:length( States_convergence)
                    CurrentState =  States_convergence(i);
                    AbsolutePhase = mod(CurrentState,3);  
                    if i <  idx_CyclingG2M(1) %not affected by treatment
                        continue
                    elseif i < leftout_CyclingG2M(end) + 3
                        continue
                    %elseif ismember(i,  leftout_CyclingG2M)
                        %continue
                    elseif   ismember(i, idx_CyclingG2M) %&&  ~ismember(i,  leftout_CyclingG2M)
                        %% G2M
                        alpha = [a(3)  a(3)*factor_G2M_1  a(3)*factor_G2M_2 ]; %] a(3)*factor_G2M_3   ];
                        beta = [b_scale   b_scale  b_scale  b_scale];
                        p1_G2M = 0.3;
                        pesos = [p1_G2M, p2_G2M, p3_G2M ];% p4_G2M ];
                        pesos =   pesos./sum( pesos);
                        tilde_G2M_pdf =  0;
                        for conta = 2 :2
                            %plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
                            tilde_G2M_pdf =  tilde_G2M_pdf + pesos(conta)*gampdf(t_span,alpha(conta),beta(conta));
                        end
                        tilde_G2M_cdf = zeros(1,num_t);
                        for t = 2:num_t
                            tilde_G2M_cdf(t)   = trapz(0:StepSize: t_span(t) ,  tilde_G2M_pdf(1:t));
                        end
                       
                        mat_G( CurrentState, CurrentState, :   ) =    tilde_G2M_cdf;
                        mat_Gm_tilde_mixture( CurrentState, CurrentState+1 ,: ) =     tilde_G2M_cdf.* m31*( 1-d_G2M);
                        I_mat_G_tilde_mixture(CurrentState, CurrentState,: ) = 1 -  tilde_G2M_cdf; 
                        mat_Gm_tilde_mixture(CurrentState, end,:) =   tilde_G2M_cdf * d_G2M; %mat_m(i,j) 

                    else
                        alpha = [a( AbsolutePhase )   a( AbsolutePhase  )*factor_rest( AbsolutePhase,1)   a( AbsolutePhase  )*factor_rest( AbsolutePhase,2)];
                        beta = [b_scale   b_scale   b_scale ];
                        pesos = [0.5  p_rest(AbsolutePhase,1)   p_rest(AbsolutePhase,2)];
                        pesos = pesos./sum(pesos);
                        %pesos = [4/5, 1/10, 1/10];
                        tilde_pdf =  0;
                        for conta = 2: 2
                            %plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
                            tilde_pdf =  tilde_pdf + pesos(conta)*gampdf(t_span,alpha(conta),beta(conta));
                        end
                        tilde_cdf = zeros(1,num_t);
                        for t = 2:num_t
                            tilde_cdf(t)   = trapz(0:StepSize: t_span(t) ,  tilde_pdf(1:t));
                        end
                         mat_G( CurrentState, CurrentState, :   ) = tilde_cdf ;
                         mat_Gm_tilde_mixture( CurrentState, CurrentState+1 ,: ) = tilde_cdf .*m_col( AbsolutePhase );   %mat_m( CurrentState, CurrentState+1 );
                         I_mat_G_tilde_mixture(CurrentState, CurrentState,: ) = 1 - tilde_cdf ; 
                         mat_Gm_tilde_mixture(CurrentState, end,:) =  tilde_cdf *d_col(AbsolutePhase) ; %mat_m(i,j) 
                    
                    end
                end
            %end
        else
            mat_Gm_tilde_mixture =mat_Gm;
            I_mat_G_tilde_mixture = I_mat_G;
            
        end
        
        
end
end