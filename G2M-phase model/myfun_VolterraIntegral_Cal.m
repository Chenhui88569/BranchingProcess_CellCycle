function Output_vec = myfun_VolterraIntegral_Cal(FirstTerm, InitValue,  LS_SecondTerm, t_span, StepSize)
num_t = length(t_span);
Output_vec = zeros(1,num_t);
Output_vec(1) = InitValue ;
n_try =20;
for i =  2:num_t
      Output_vec(i) = Output_vec(i-1); %The initial estimate for the iteration
      k_vec =  KernelGen(i,1:i).*Output_vec(1:i);
      for j = 1:n_try
          Output_vec(i)  =  FirstTerm(i) + StepSize*( sum( k_vec(2:i-1) ) +  (k_vec(1) +  k_vec(i))/2  );
          if  Output_vec(i) <0
              Output_vec(i) = 0;
          end
          k_vec(i) =   KernelGen(i,i).*Output_vec(i);
      end
end

    function f_k =  KernelGen(x_i, x_1toi)
        f_k = zeros(size(x_1toi));
        for t = 1: length(x_1toi)
            f_k(t) = LS_SecondTerm(x_i - x_1toi(t) +1 );
       end
    end
end


