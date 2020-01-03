
%Adesoji Bello
%MAP Posterior estimate
%Absolute Value estimate
%Mean Square estimate
syms tf Tf U y(tf)
%------PDF for the estimated value (tf) | t>= 0

f_tf = (1/100) * exp (-tf/100);

%---------PDF for Observation (Tf | tf)----------
fTf_tf2 =400*tf^2/Tf^3 ;                        %0<= tf <= 0.05*Tf
fTf_tf1 =((Tf-10*tf)*tf*40)/Tf^3 ;              %0.05*Tf<= tf <= 0.10*Tf

%-----------PDF for Tf at both boundary conditions
fTf = int(fTf_tf1*f_tf,tf,[0.05*Tf 0.1*Tf]) +int(fTf_tf2*f_tf,tf,[0 0.05*Tf]);

%PDF for Radio Usage

f_U1 = U/25 ;    %0<U<5
f_U2 = (10-U)/25; %5<U<10;
%--------A Posterior PDF f(tf|Tf)--------------

ftf_Tf1 = (fTf_tf1 * f_tf) /fTf;
ftf_Tf2 = (fTf_tf2 * f_tf) /fTf;
Tf_sub = 1000;
ftf_Tf1k = subs(ftf_Tf1,Tf,Tf_sub);
ftf_Tf2k = subs(ftf_Tf2,Tf,Tf_sub);

tfv = linspace(0,100,100);
y = piecewise(tf> 0 & tf <=0.05*Tf_sub, ftf_Tf2k, tf > 0.05*Tf_sub & tf <0.1*Tf_sub, ftf_Tf1k);
y_val = subs(y,tf,tfv);
y_val = vpa(y_val);


%--------------Figure 1------------------------------------------
% figure(1)
% plot(tfv,y_val)
% title('Plot of A Posterior Density Function f(tf | Tf) at Tf = 10000 ')
% xlabel('tf')
% %ylabel('f(tf|Tf)')
%grid on


%-------------Using the MEAN SQUARE Estimate condition------------
Tf_val = [1000,5000,10000];
est_ms= zeros(1,3);
for i=1:1:3
 est_ms(i) = int(tf*subs(ftf_Tf2,Tf,Tf_val(i)),tf,[0 0.05*Tf_val(i)])+ int(tf*subs(ftf_Tf1,Tf,Tf_val(i)),tf,[0.05*Tf_val(i) 0.10*Tf_val(i)]);            

end
disp('The Mean Square Estimate for the 3 values of Tf are:');
vpa(est_ms)

%-------------Using the Absolute Value Estimate condition------------
Tf_val = [1000,5000,10000];
tf_ABS = zeros(1,3);
for i =1:1:3
    ftf_Tf1k = subs(ftf_Tf1,Tf,Tf_val(i));
    ftf_Tf2k = subs(ftf_Tf2,Tf,Tf_val(i));
    y = piecewise(tf>= 0 & tf <=0.05*Tf_val(i), ftf_Tf2k, tf > 0.05*Tf_val(i) & tf <=0.1*Tf_val(i), ftf_Tf1k);
    tf_aval = linspace(0,Tf_val(i)*0.1,Tf_val(i)*0.1);
    y_v = subs(y,tf,tf_aval);
    tf_abs = cumsum(vpa(y_v))';
    
    for j=1:1:length(tf_abs)
        if tf_abs(j) >=0.5
             kval = j-1;            %Lowers the step size
             break;
         else
            continue;
        end 
    end
    tf_ABS(i) = kval;
end

disp('The Absolute Value Estimates for the 3 values are:');
tf_ABS

%----------FINDING the MAP Estimates ------------------------------------------------
Tf_valm = [1000,5000,10000];                %Substituting Values for Tf
tf_MAP = zeros(1,3);
for i =1:1:3
    ftf_Tf1k = subs(ftf_Tf1,Tf,Tf_valm(i));
    ftf_Tf2k = subs(ftf_Tf2,Tf,Tf_valm(i));
    y = piecewise(tf>= 0 & tf <=0.05*Tf_valm(i), ftf_Tf2k, tf > 0.05*Tf_valm(i) & tf <0.1*Tf_valm(i), ftf_Tf1k);
    tfv1 = linspace(0,Tf_valm(i)*0.1,1000);
    y_v = vpa(subs(y,tf,tfv1));
    [m,k] = max(y_v);
    tf_MAP(i) = tfv1(k);
end
disp('The MAP Estimate for the 3 values are:');
vpa(tf_MAP)


%--------Generate 10000 values of tf & Tf(Mean Time to Failure)----------------
N = 10000;                %Number of samples
tf_val = exprnd(100,1,N);
U = randi([1 10],1,N);

Tf_val10k = zeros(1,N);

for i = 1:1:N
    Tf_val10k(i) = 100*tf_val(i)/U(i);
end

% -------------------Mean & Variance Estimation Error Computation for MSE--------------
start = tic;
est_ms10k = zeros(1,N);
est_int = int(tf*ftf_Tf2,tf,[0 0.05*Tf])+ int(tf*ftf_Tf1,tf,[0.05*Tf 0.10*Tf]);
for i = 1:1:N
    
 est_ms10k(i) = vpa(subs(est_int,Tf,Tf_val10k(i)));
end

fprintf('The Mean error estimate for MSE with N = %d Tf is:',N);
Mean_MSE = (1/N) * sum(tf_val-est_ms10k)
var = 0;
for i =1:1:N
 var = var + ((tf_val(i)-est_ms10k(i)) - Mean_MSE).^2;

end
fprintf('The Variance error estimate for MSE with N = %d Tf is:',N);
var_MSE = (1/(N-1)) * var
%Cost Function
fprintf('The average absolute value cost & MSE cost for MSE with N = %d Tf is:',N);
J_MSE = (1/N-1) * sum(tf_val-est_ms10k)^2
J_MSEabs = (1/N-1) * sum(abs(tf_val-est_ms10k))

toc(start);


% -------------------Mean & Variance Computation for Absolute Value--------------
start = tic;
est_abs10k = zeros(1,N);
 for i =1:1:N
     ftf_Tf1k = subs(ftf_Tf1,Tf,Tf_val10k(i));
    ftf_Tf2k = subs(ftf_Tf2,Tf,Tf_val10k(i));
    y = piecewise(tf>= 0 & tf <=0.05*Tf_val10k(i), ftf_Tf2k, tf > 0.05*Tf_val10k(i) & tf <=0.1*Tf_val10k(i), ftf_Tf1k);
    tf_aval = linspace(0,Tf_val10k(i)*0.1,round(Tf_val10k(i)*0.1));
    y_v = subs(y,tf,tf_aval);
    tf_abs = cumsum(vpa(y_v))';
  
    for j=1:1:length(tf_abs)
       if tf_abs(j) >=0.5
            kval = j;
            break;
        else
            continue;
        end  
    end
   est_abs10k(i) = kval;
end

fprintf('The Mean error estimate for Absolute Value with N = %d Tf is:',N);
Mean_Abs = (1/N) * sum(tf_val-est_abs10k)
var = 0;
for i =1:1:N
 var = var + ((tf_val(i)-est_abs10k(i)) - Mean_Abs)^2;

end
fprintf('The Variance error estimate for Absolute Value with N = %d Tf is:',N);
var_Abs = (1/(N-1)) * var
%cost Function
fprintf('The average absolute value cost & MSE cost for Absolute Value with N = %d Tf is:',N);
J_ABS = (1/N-1) * sum(abs(tf_val-est_abs10k))
J_ABSab = (1/N-1) * sum(abs(tf_val-est_abs10k))



%----------FINDING the MAP Estimates ------------------------------------------------

Using Tf_val10k
est_map10k = zeros(1,N);
for i =1:1:N
    ftf_Tf1k = subs(ftf_Tf1,Tf,Tf_val10k(i));
    ftf_Tf2k = subs(ftf_Tf2,Tf,Tf_val10k(i));
    y = piecewise(tf> 0 & tf <=0.05*Tf_val10k(i), ftf_Tf2k, tf > 0.05*Tf_val10k(i) & tf <0.1*Tf_val10k(i), ftf_Tf1k);
    tfv1 = linspace(0,Tf_val10k(i)*0.1,100);
    y_v = subs(y,tf,tfv1);
    y_v1 = vpa(y_v);
    [val, k] = max(y_v1);
    est_map10k(i) = tfv1(k);
end
fprintf('The Mean error estimate for MAP with N = %d Tf is:',N);
Mean_MAP = (1/N) * sum(tf_val-est_map10k)
var = 0;
for i =1:1:N
 var = var + ((tf_val(i)-est_map10k(i)) - Mean_MAP).^2;

end
fprintf('The Variance error estimate for MAP with N = %d Tf is:',N);
var_MAP = (1/(N-1)) * var
%Average Absolute Value Cost & MSE Cost
fprintf('The average absolute value cost & MSE cost for MAP with N = %d Tf is:',N);
J_MAP = (1/N-1) * sum(tf_val-est_map10k)^2
J_MAPabs = (1/N-1) * sum(abs(tf_val-est_map10k))

toc(start);
