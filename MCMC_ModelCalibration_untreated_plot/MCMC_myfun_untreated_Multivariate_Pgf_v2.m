function [SteadyState,  ClinicalTime, CellFractions , t_span] = MCMC_myfun_untreated_Multivariate_Pgf_v2(theta )

% newdir  = 'Figs_untreated_v2';
% if ~isfolder(newdir)
%     mkdir(newdir  );
% end

%% Initialization
ChoiceOfDist = 'Gamma'; %
MaxT =  48*8; %unit in h
StepSize= 0.1;
a =  theta(1:3);
b =  theta(4);
b_scale = 1/b;
N_0_1 =     1.0845e+12 * 0.3847 ;  %1.0845e+12 *   0.4126; % initial number of cells
N_0_2 =  0 ;                         %1.0845e+12 * 0.3847 ;
N_0_3 =  0  ;                      %1.0845e+12*     0.2028;
N_0  = [N_0_1; N_0_2 ; N_0_3];
t_span = 0: StepSize: MaxT;
if t_span(end) ~=  MaxT
    t_span  = [t_span MaxT];
end
%% G1
[G1_pdf, G1_cdf]  =  PdfCdfCal(ChoiceOfDist,  t_span , a(1) ,b_scale);
%% S
[S_pdf, S_cdf]  =  PdfCdfCal(ChoiceOfDist,  t_span , a(2) ,b_scale);
%% G2
[G2_pdf, G2M_cdf]  =  PdfCdfCal(ChoiceOfDist,  t_span , a(3) ,b_scale);

H_G1_cdf =  G1_cdf;
H_S_cdf =  S_cdf;
H_G2M_cdf =  G2M_cdf;
%% 
m12 = 1 ;
m23 = 1;
m31 = 2;

num_t = length(t_span); 
IncreaseGen_flag = 0;
if  IncreaseGen_flag 
    num_gen = 45;
else
    num_gen = 30;
end

%% matrix m
temp = diag( repmat([m12, m23, m31],1,num_gen )  ,1);
%mat_m= temp;
mat_m= [temp(1: 3*num_gen, 1: 3*num_gen),   zeros(3*num_gen,1 )    ; 
    zeros(1, 3*num_gen+1 )  ];


%% matrix G
mat_G = zeros(3*num_gen+1, 3*num_gen+1 , num_t);
mat_G(1,1,:) = H_G1_cdf;
mat_G(2,2,:) = H_S_cdf;
mat_G(3,3,:) = H_G2M_cdf;
%single dose
for k = 2:num_gen
    mat_G( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = G1_cdf;
    mat_G( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = S_cdf; 
    mat_G( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = G2M_cdf;      
end
%mat_G(end,  end,:   )  =   G1_cdf;      

%% matrix I-G
I_mat_G = zeros(3*num_gen+1 , 3*num_gen+1, num_t);
I_mat_G(1,1,:) = 1-H_G1_cdf;
I_mat_G(2,2,:) = 1-H_S_cdf;
I_mat_G(3,3,:) = 1-H_G2M_cdf;
for k = 2:num_gen
        I_mat_G( 1+ 3*(k-1),  1+ 3*(k-1) ,:   ) = 1-  G1_cdf;
        I_mat_G( 2+ 3*(k-1),  2+ 3*(k-1),:   )  = 1-  S_cdf; 
        I_mat_G( 3+ 3*(k-1),  3+ 3*(k-1),:   )  = 1-  G2M_cdf; 
end
%I_mat_G(end,end,:) = 1- G1_cdf;
%% matrix Gm
mat_Gm = zeros(3*num_gen+1,3*num_gen+1, num_t);
for i = 1: size(mat_m,1)
        if  i < 3*num_gen-1
            mat_Gm(i,i+1 ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t  ) * mat_m(i,i+1); %mat_m(i,j) scalar value
        %mat_Gm(i, end ,: ) = reshape( mat_G( i,  i,:   ) , 1, num_t  ) * mat_m(i,end); %mat_m(i,j) scalar value
        end
end
%% 
% stopping rule
ctr = 0; % the number of iterates
while ctr<=  3*num_gen
    if ctr  == 0
        % k=0
        next_mat_M = I_mat_G;
    elseif ctr == 1    % k = 1
        curr_mat_M = next_mat_M;
        next_T_k_iterate = mat_Gm;
        next_mat_M_term_k  =   zeros(3*num_gen+1,3*num_gen+1, num_t);
        for i = 1: size(mat_m,1) 
            for  j  = 1: size(mat_m,1)
                 temp_mat_M2_conv_full  =  conv( reshape(I_mat_G(j,j,:),1,num_t),   TakeDerivative(reshape(  next_T_k_iterate(i,j,:),1,num_t) ),'full')*StepSize;
                % temp_mat_M2_conv_full =  conv( reshape( next_T_k_iterate(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full')*StepSize;
                 temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                next_mat_M_term_k(i,j,:) =   temp_mat_M2_conv_same ;
            end
        end
        next_mat_M = curr_mat_M+next_mat_M_term_k;
    else  %k >= 2
        curr_mat_M = next_mat_M;
        curr_T_k_iterate = next_T_k_iterate;
        iterate_T_k =  zeros(3*num_gen+1,3*num_gen+1 , num_t);
        %Calculate current iterate T^k
        for i = 1: size(mat_m,1)
            %for  j  = 1: size(mat_m,1)
               % InterVar = 0 ;
               % for m = 1: size(mat_m,1)
               if i+ctr <  size(mat_m,1) % curr_T_k_iterate(i,i+ctr-1,:),1,num_t) is non-zero, computed in the pevious interation.
                   temp_mat_M1_conv_full = conv( reshape( curr_T_k_iterate(i,i+ctr-1,:),1,num_t),   TakeDerivative(reshape(mat_Gm(i+ctr-1,i+ctr,:),1,num_t)),'full')*StepSize;
                   temp_mat_M1_conv_same =  temp_mat_M1_conv_full(1:num_t);
                   %temp_mat_M1_conv = Conv_Cal( reshape( curr_T_k_iterate(i,m,:),1,num_t),  reshape(mat_Gm(m,j,:),1,num_t),t_span  );
                   %InterVar  = InterVar  + temp_mat_M1_conv_same;
                    %end
                   iterate_T_k(i,i+ctr,:) = temp_mat_M1_conv_same;       
                  
               end
               
           % end
        end
        next_mat_M_term_k  =   zeros(3*num_gen+1,3*num_gen +1, num_t);
        for i = 1: size(mat_m,1)
            %for  j  = 1: size(mat_m,1)
                if i+ctr <  size(mat_m,1) % curr_T_k_iterate(i,i+ctr-1,:),1,num_t) is non-zero, computed in the pevious interation.
                    temp_mat_M2_conv_full =  conv( reshape(I_mat_G(i+ctr,i+ctr,:),1,num_t),   TakeDerivative(reshape(  iterate_T_k(i,i+ctr,:),1,num_t) ),'full')*StepSize;
                    % temp_mat_M2_conv_full =  conv( reshape(  iterate_T_k(i,j,:),1,num_t),  reshape(I_mat_G(j,j,:),1,num_t) ,'full')*StepSize;
                    temp_mat_M2_conv_same =   temp_mat_M2_conv_full(1:num_t);
                    next_mat_M_term_k(i,i+ctr,:) =  temp_mat_M2_conv_same ;
                end
            %end
        end
        next_T_k_iterate  = iterate_T_k ; %update iterate
        next_mat_M  = curr_mat_M + next_mat_M_term_k;

    end  
    ctr = ctr +1;
end

 
%% 
mat_M =  next_mat_M;
max_idx =  1;

N_col_G1_AcrossGen = zeros( num_gen, num_t -  max_idx+1 );
N_col_S_AcrossGen = zeros( num_gen, num_t - max_idx+1 );
N_col_G2M_AcrossGen = zeros( num_gen, num_t- max_idx+1 );

for k = 1:num_gen
     
    N_col_G1_AcrossGen(k,:) =  N_0(1)* reshape(mat_M(1, 3*(k-1)+1 ,:), 1, num_t) +  N_0(2)*reshape(mat_M(2,3*(k-1)+1,:), 1,num_t) + N_0(3)*reshape(mat_M(3,3*(k-1)+1,:), 1, num_t) ;
    N_col_S_AcrossGen(k,:) = N_0(1)* reshape(mat_M(1,3*(k-1)+2,:), 1, num_t) +  N_0(2)*reshape(mat_M(2,3*(k-1)+2,:), 1, num_t) + N_0(3)*reshape(mat_M(3,3*(k-1)+2,:), 1, num_t) ;
    N_col_G2M_AcrossGen(k,: )= N_0(1)* reshape(mat_M(1,3*(k-1)+3,:), 1, num_t) +  N_0(2)*reshape(mat_M(2,3*(k-1)+3,:), 1,num_t) + N_0(3)*reshape(mat_M(3,3*(k-1)+3,:), 1, num_t) ;

end
N_col_DiffPhases =  [ sum( N_col_G1_AcrossGen,1)   ;
    sum( N_col_S_AcrossGen,1);
    sum( N_col_G2M_AcrossGen,1)];
N_total = sum( N_col_DiffPhases,1  ) ; 

%% 
Threshold = 10e-2;
CellFractions = N_col_DiffPhases./repmat(N_total,3,1);

idx_5h = find(t_span == 250 );

t_idx = idx_5h + 1;

FlagSteadyState = false;
while ~FlagSteadyState && t_idx <= length(t_span)
    curr_t = t_span(t_idx);
    % find the cell fractions 5h before t_idx
    idx_MaxTminus5  = find(t_span>   curr_t -5 ,1); 
    CellFractions_5h = CellFractions(:,idx_MaxTminus5  );
    
    if  sum( abs(CellFractions_5h - CellFractions(:,t_idx)) <= Threshold  ) == 3
        SteadyState = CellFractions(:,t_idx);
        FlagSteadyState = true;
    else
        FlagSteadyState = false;
        SteadyState = [];
        %     set(0,'currentfigure' ,f_FalsePlot )
        %     plot(CellFractions','LineWidth',1.5)
        %     set(gca,'FontSize',14);
        %     legend(["G1" ,"S" ,"G2M"])
        %     xlabel('Time(Hours)')
        %     ylabel('Cell fraction')
        %     title( 'Cell fractions' )
    end
    t_idx  = t_idx +1;
end
ClinicalTime = t_idx-1;  %idx



%t_span_new = 0:StepSize:t_span( num_t - max_idx +1 );
% f = figure;
% set(f, 'Position', get(0, 'Screensize'));
% subplot(1,3,1)
% semilogy(t_span_new,   N_col_G1_AcrossGen , 'LineWidth',1.5);
% %plot(t_span,   N_col_G1_AcrossGen );
% num_col = string(1:1:num_gen);
% legend(append("n", num_col),'Location','best');
% xlabel('Time(Hours)')
% ylabel('Cell number')
% set(gca,'FontSize',14);
% grid on
% %ylim([0  max( max(N_col_G1_AcrossGen))])
% title('G1')
% 
% subplot(1,3,2)
% semilogy(t_span_new,  N_col_S_AcrossGen, 'LineWidth',1.5);
% %plot(t_span,   N_col_S_AcrossGen );
% num_col = string(1:1:num_gen);
% legend(append("n", num_col),'Location','best');
% xlabel('Time(Hours')
% ylabel('Cell number')
% set(gca,'FontSize',14);
% grid on
% title('S')
% 
% subplot(1,3,3)
% semilogy(t_span_new, N_col_G2M_AcrossGen, 'LineWidth',1.5 );
% %plot(t_span,   N_col_G2M_AcrossGen );
% num_col = string(1:1:num_gen);
% legend(append("n", num_col),'Location','best');
% xlabel('Time(Hours)')
% ylabel('Cell number')
% set(gca,'FontSize',14);
% grid on
% title('G2/M')
% sgtitle(strcat(num2str(num_gen)," generations of offspring after the devision of 0th generation" ) ,'FontSize', 14 ) 
% exportgraphics(f, fullfile(newdir , strcat("CellsInDifferentGens_treated_gen" , num2str(num_gen), '.tif' )  ), 'Resolution', 300 )



    function output = TakeDerivative(T)
        T_temp = diff(T )./StepSize;
        output  =   [ T_temp  T_temp(end)]; 
           %output  = T;
    end
end
