clc
clear all
%calculations of h
% DT = 250; %C
% L = 7 ; %m
% h = 1.52*(DT)^(1/3); % W/m^2.K
%به این دلیل که جریان در تنور متلاطم است از رابطه ی بالا استفاده میکنیم
k = 0.223; % W/m.K
ro = 720; % Kg/m^3
Cp = 3600; % J/Kg.K
h = 9.57; % W/m^2.K
Bi = h*0.0025/k;
velocity= 0.1667; % m/s
delta_t = 0.0150; % s
alpha = k/(ro*Cp);
fo = alpha*delta_t/0.0025;
q = 1350; %W
V = 0.00025; %m^3
q_dot = q/V;
Q = q/(ro*Cp*V);
q_w = 20000; % W/m^2
st = 5.67*10^-8; % W/m^2.K^4
%T_inf calculations
T_inf_x = @(x) ((250/7)*x + 50);
x = linspace(0,7,2801);
T_inf = zeros(1,2801);
for i = 1:2801
    T_inf(i) = T_inf_x(x(i));
end
T_inf2 = zeros(1 , 2800);
for i = 1:2800
T_inf2(i) =T_inf(2801-i);
end
T_inf3 = [T_inf , T_inf2 , T_inf ];
T = zeros(3,201,5601);
T0 = ones(3,201)*25;
T(:,:,1) = T0;
% P=1;
% claculations for temperature  of each node in each time (khaste 1)
for P = 1:5600
    T_inf4 = T_inf3(1+P:201+P);
    coef_T = zeros(603,603);
    csc_T = zeros(603,1);
    T_first = T(:,:,P);
    %left_down
    csc_T(1,1) = T_first(1,1) + q_dot*delta_t/(ro*Cp) + 2*h*delta_t/(ro*0.0025*Cp)*T_inf4(1) + 2*q_w*delta_t/(ro*Cp*0.0025);
    coef_T(1,1) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
    coef_T(1,2) = -2*fo;
    coef_T(1,202) = -2*fo;
    %center points in down plate
    for i = 2:200
        csc_T(i,1) = T_first(1,i) + q_dot*delta_t/(ro*Cp) + 2*q_w*delta_t/(ro*Cp*0.0025);
        coef_T(i,i) = 1 + 4*fo ;
        coef_T(i,i-1) = -fo;
        coef_T(i,i+1) = -fo;
        coef_T(i,i+201) = -2*fo;
    end
    %right_down
    csc_T(201,1) = T_first(1,201) + q_dot*delta_t/(ro*Cp) + 2*q_w*delta_t/(ro*Cp*0.0025) + 2*h*delta_t/(ro*Cp*0.0025)*(T_inf4(201));
    coef_T(201,201) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
    coef_T(201,200) = -2*fo;
    coef_T(201,402) = -2*fo;
    % mid_left
    csc_T(202,1) = T_first(2,1) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(1) + q_dot*delta_t/(ro*Cp);
    coef_T(202,202) = 4*fo + 2*h*delta_t/(ro*Cp*0.0025) + 1;
    coef_T(202,1) = -fo;
    coef_T(202,403) = -fo;
    coef_T(202,203) = -2*fo;
    %center points in middle plate
    for i = 203:401
       for j = 2:200
           csc_T(i,1) = T_first(2,j) + q_dot*delta_t/(ro*Cp);
           coef_T(i,i) = 4*fo+1;
           coef_T(i,i+1) = -fo;
           coef_T(i,i-1) = -fo;
           coef_T(i,i+201) = -fo;
           coef_T(i,i-201) = -fo;
       end
   end
   %mid_right
   csc_T(402,1) = q_dot*delta_t/(ro*Cp) + T_first(2,201) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(201);
   coef_T(402,402) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
   coef_T(402,201) = -fo;
   coef_T(402,603) = -fo;
   coef_T(402,401) = -2*fo;
   %up_left
   csc_T(403,1) = 4*h*delta_t/(ro*Cp*0.0025)*T_inf4(1) + q_dot*delta_t/(ro*Cp) + T_first(3,1);
   coef_T(403,403) = 4*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
   coef_T(403,404) = -2*fo;
   coef_T(403,202) = -2*fo;
   %center points in up plate
   for i = 404:602
       for j = 2:200
           csc_T(i,1) = q_dot*delta_t/(ro*Cp) + T_first(3,j) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(j);
           coef_T(i,i) = 2*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
           coef_T(i,i+1) = -fo;
           coef_T(i,i-1) = -fo;
           coef_T(i,i-201) = -2*fo;
      end
  end
  %up_right
  csc_T(603,1) = 4*h*delta_t/(ro*Cp*0.0025)*T_inf4(201) + q_dot*delta_t/(ro*Cp) + T_first(3,201);
  coef_T(603,603) = 4*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
  coef_T(603,602) = -2*fo;
  coef_T(603,402) = -2*fo;
  temperatures = inv(coef_T)*csc_T;
  for i = 1:3
    for j = 1:201
        T(i,j,P+1) = temperatures((((i-1)*201)+j),1);
    end
  end
  disp(P)
end

disp(T(:,:,5600))
% plotting (khaste 3)
g = zeros(5600,201);
for i = 2:5600
    g(1,:) = 25;
    g(i,:) = T(2,:,i);
end
mesh(g)
% calculating velocity when we have radiation (khaste 2)
% velocity_prime = 0.112; % m/s
% delta_t = 0.0223;%s
% alpha = k/(ro*Cp);
% fo = alpha*delta_t/0.0025;
% for P = 1:5600
%     T_inf4 = T_inf3(1+P:201+P);
%     coef_T = zeros(603,603);
%     csc_T = zeros(603,1);
%     T_first = T(:,:,P);
%     %left_down
%     csc_T(1,1) = T_first(1,1) + q_dot*delta_t/(ro*Cp) + 2*h*delta_t/(ro*0.0025*Cp)*T_inf4(1) + 2*q_w*delta_t/(ro*Cp*0.0025);
%     coef_T(1,1) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
%     coef_T(1,2) = -2*fo;
%     coef_T(1,202) = -2*fo;
%     %center points in down plate
%     for i = 2:200
%         csc_T(i,1) = T_first(1,i) + q_dot*delta_t/(ro*Cp) + 2*q_w*delta_t/(ro*Cp*0.0025);
%         coef_T(i,i) = 1 + 4*fo ;
%         coef_T(i,i-1) = -fo;
%         coef_T(i,i+1) = -fo;
%         coef_T(i,i+201) = -2*fo;
%     end
%     %right_down
%     csc_T(201,1) = T_first(1,201) + q_dot*delta_t/(ro*Cp) + 2*q_w*delta_t/(ro*Cp*0.0025) + 2*h*delta_t/(ro*Cp*0.0025)*(T_inf4(201));
%     coef_T(201,201) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
%     coef_T(201,200) = -2*fo;
%     coef_T(201,402) = -2*fo;
%     % mid_left
%     csc_T(202,1) = T_first(2,1) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(1) + q_dot*delta_t/(ro*Cp);
%     coef_T(202,202) = 4*fo + 2*h*delta_t/(ro*Cp*0.0025) + 1;
%     coef_T(202,1) = -fo;
%     coef_T(202,403) = -fo;
%     coef_T(202,203) = -2*fo;
%     %center points in middle plate
%     for i = 203:401
%        for j = 2:200
%            csc_T(i,1) = T_first(2,j) + q_dot*delta_t/(ro*Cp);
%            coef_T(i,i) = 4*fo+1;
%            coef_T(i,i+1) = -fo;
%            coef_T(i,i-1) = -fo;
%            coef_T(i,i+201) = -fo;
%            coef_T(i,i-201) = -fo;
%        end
%    end
%    %mid_right
%    csc_T(402,1) = q_dot*delta_t/(ro*Cp) + T_first(2,201) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(201);
%    coef_T(402,402) = 4*fo + 1 + 2*h*delta_t/(ro*Cp*0.0025);
%    coef_T(402,201) = -fo;
%    coef_T(402,603) = -fo;
%    coef_T(402,401) = -2*fo;
%    %up_left
%    csc_T(403,1) = 4*h*delta_t/(ro*Cp*0.0025)*T_inf4(1) + q_dot*delta_t/(ro*Cp) + T_first(3,1) + 2*delta_t*st/(ro*Cp*0.0025)*((T_first(3,1)+273)^4 - (T_inf4(1)+273)^4);
%    coef_T(403,403) = 4*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
%    coef_T(403,404) = -2*fo;
%    coef_T(403,202) = -2*fo;
%    %center points in up plate
%    for i = 404:602
%        for j = 2:200
%            csc_T(i,1) = q_dot*delta_t/(ro*Cp) + T_first(3,j) + 2*h*delta_t/(ro*Cp*0.0025)*T_inf4(j) +  4*delta_t*st/(ro*Cp*0.0025)*((T_first(3,j)+273)^4 -(T_inf4(j)+273)^4);
%            coef_T(i,i) = 2*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
%            coef_T(i,i+1) = -fo;
%            coef_T(i,i-1) = -fo;
%            coef_T(i,i-201) = -2*fo;
%       end
%   end
%   %up_right
%   csc_T(603,1) = 4*h*delta_t/(ro*Cp*0.0025)*T_inf4(201) + q_dot*delta_t/(ro*Cp) + T_first(3,201) +  2*delta_t*st/(ro*Cp*0.0025)*((T_first(3,201)+273)^4 - (T_inf4(201)+273)^4);
%   coef_T(603,603) = 4*h*delta_t/(ro*Cp*0.0025) + 1 + 4*fo;
%   coef_T(603,602) = -2*fo;
%   coef_T(603,402) = -2*fo;
%   temperature = inv(coef_T)*csc_T;
%   for i = 1:3
%     for j = 1:201
%         T(i,j,P+1) = temperature((((i-1)*201)+j),1);
%     end
%   end
%   disp(P)
% end
% disp(T(:,:,5600))
% بنابراین سرعت در بخش الف برابر 16.67 دلتا تی برابر0.015 شد
% سرعت در بخش ب برابر 11.2 و دلتا تی برابر 0.0223 شد
% تابش را فقط از سطح بالایی نان و با فرض این که توزیع دمای سقف نان با توزیع
% دمای هوای اطراف یکی است انجام دادیم
% velocity at part a = 16.7 cm/s and delta_t = 0.015
%velocity at part b =  11.2  cm/s  and delta_t = 0.0223
% سرعت را بصورت دستی عوض کردیم و انقدر تغییر دادیم تا در نهایت در لحظه ی
% اخر دمای تمام نقاط به 200 درجه ی سانتی گراد رسید
% برای بدست آوردن توزیع دمای خواسته ی 2 لطفا کد های مربوط به این بخش را از
% حالت کامنت درآورید

   

