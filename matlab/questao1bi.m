clc, clear, close all

L1 = Link('revolute', 'd', .290, 'alpha', -pi/2, 'qlim', 11/12*[-pi pi]);

L2 = Revolute('a', .270, 'offset', -pi/2, 'qlim',11/18*[-pi pi]);

L3 = Revolute('a', .070, 'alpha', -pi/2, 'qlim', [-(11/18)*pi (7/18)*pi]);

L4 = Revolute('d', .302, 'alpha', pi/2, 'qlim', 8/9*[-pi pi]);

L5 = Revolute('alpha', -pi/2, 'qlim', 2/3*[-pi pi]);

L6 = Revolute('d', 0.072, 'offset', pi, 'qlim', 20/9*[-pi pi]);

robo1 = SerialLink([L1 L2 L3 L4 L5 L6], 'name', 'IRB120')

q = [0 0 0 0 -pi/2  0];

qvlim = [25*pi/18 25*pi/18 25*pi/18 16*pi/9 16*pi/9 7*pi/3];
%____________________________


%_________________________________--
syms t;
wn = pi/10;

pds(t) = [0.050*sin(wn*t)+0.428 0.02 0.050*cos(wn*t)+0.669];
%pds(t) = [0.050*sin(wn*t)+0.428 0.2 0.050*cos(wn*t)+0.669];
%pds(t) = [0.020*(sin(wn*t)+sin(4*wn*t))+0.428 0.02 0.020*(cos(wn*t)+cos(4*wn*t))+0.669];

pddots = diff(pds);

figure('position',[280 70 560 420])
robo1.plot(q)
view([60,30])

pd = [0.478 0.020 0.669];
%Rd = SO3();
%Rd= Rd.R;
%Td = SE3(Rd,pd);
%rpyd = rotm2eul(Rd);

lambda = 1;
epsilon = 2e-2;
e = inf(6,1);
cont0 = tic;
t0 = tic;

contf = toc(cont0);

i=0;

%rpyd = [0,0,0]

while (norm(e) > epsilon) % Critério de parada
    T = robo1.fkine(q);
    Jc = robo1.jacob0(q, 'rpy'); % Jacobiana analítica
    J = Jc(1:3,:)
    p = transl(T);
    p_til = pd - p;
   
    R = SO3(T);
    R = R.R; % Extrai rotação do efetuador
    rpy = rotm2eul(R);
   
    %rpy_til = rpyd-rpy;
   
    e = [p_til']; % Vetor de erro
       
    u = pinv(J)*lambda*e; % Lei de controle

    tf= toc(t0);
    t0= tic;
   
    for k = 1:6
        if u(k) > qvlim(k)
            u(k) = qvlim(k);
        elseif u(k) < -qvlim(k)
            u(k) = -qvlim(k);
        end
       
       
        if q(k)+u(k)*tf < robo1.qlim(k,1)
            u(k) = (robo1.qlim(k,1) - q(k))/tf;
        elseif q(k)+u(k)*tf > robo1.qlim(k,2)
            u(k) = (robo1.qlim(k,2) - q(k))/tf;
        end
    end
   
    tf_robo =  @(t,qi) [u(1);u(2);u(3);u(4);u(5);u(6)];
    [~, qi] = ode45(tf_robo, 0:tf:tf, zeros(1,6));
   
    control_sig(:,i+1) = 180*u/pi;
    err(i+1)= norm(e);
    traj(:,i+1) = p;
    vrpy(:,i+1) = 180*rpy/pi;
    erro(:,i+1) = e;
   
    for k = 1:6
        qplot(k,i+1) = 180*q(k)/pi;
        qmax(k,i+1) = 180*robo1.qlim(k,2)/pi;
        qmin(k,i+1) = 180*robo1.qlim(k,1)/pi;
        qvlimmax(k,i+1) = 180*qvlim(k)/pi;
        qvlimmin(k,i+1) = -180*qvlim(k)/pi;
    end
   
    time(i+1) = i;
   
    q = q+qi(end,:);
           
    robo1.plot(q);
   
     
     
    i=i+1;
end

for j = 1:1:36
    T = robo1.fkine(q);
    Jc = robo1.jacob0(q, 'rpy'); % Jacobiana analítica
    J = Jc(1:3,:)
   
    p = transl(T);    
    pd = double(pds(j));
    p_til = pd - p;

    R = SO3(T);
    R = R.R; 
    rpy = rotm2eul(R);
    
    e = [p_til']; % Vetor de erro
       
    pddot = [double(pddots(j))];
   
    lambda= 1;
   
    u = pinv(J)*(pddot'+lambda*e); % Lei de controle

    contf = toc(cont0);
    cont0= tic;
   
    tf= toc(t0);
    t0= tic;
   
    for k = 1:6
        if u(k) > qvlim(k)
            u(k) = qvlim(k);
        elseif u(k) < -qvlim(k)
            u(k) = -qvlim(k);
        end
       
       
        if q(k)+u(k)*tf < robo1.qlim(k,1)
            u(k) = (robo1.qlim(k,1) - q(k))/tf;
        elseif q(k)+u(k)*tf > robo1.qlim(k,2)
            u(k) = (robo1.qlim(k,2) - q(k))/tf;
        end
    end
   
    tf_robo =  @(t,qi) [u(1);u(2);u(3);u(4);u(5);u(6)];
    [~, qi] = ode45(tf_robo, 0:tf:tf, zeros(1,6));
   
    control_sig(:,i+1) = 180*u/pi;
    err(i+1)= norm(e);
    traj(:,i+1) = p;
    vrpy(:,i+1) = 180*rpy/pi;
    erro(:,i+1) = e;
    
    for k = 1:6
        qplot(k,i+1) = 180*q(k)/pi;
        qmax(k,i+1) = 180*robo1.qlim(k,2)/pi;
        qmin(k,i+1) = 180*robo1.qlim(k,1)/pi;
        qvlimmax(k,i+1) = 180*qvlim(k)/pi;
        qvlimmin(k,i+1) = -180*qvlim(k)/pi;        
    end
   
    time(i+1) = i;
   
    q = q+qi(end,:);
           
    robo1.plot(q);
    
    i=i+1;
end   
hold off

%% Plot sinal de controle e norma do erro

hold on
figure(2)
sgtitle('Sinais de controle(velocidade das juntas)')
subplot(2,3,1)
hold on
grid on
title('Junta 1')
plot(time,qvlimmax(1,:),'r')
plot(time,qvlimmin(1,:),'g')
plot(time,control_sig(1,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_1(graus/s)')
legend('qvlim_[1max]', 'qvlim_[1min]','u_1','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(1) 1.1*(180/pi)*qvlim(1)])

subplot(2,3,2)
hold on
grid on
title('Junta 2')
plot(time,qvlimmax(2,:),'r')
plot(time,qvlimmin(2,:),'g')
plot(time,control_sig(2,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_2(graus/s)')
legend('qvlim_[2max]', 'qvlim_[2min]','u_2','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(2) 1.1*(180/pi)*qvlim(2)])

subplot(2,3,3)
hold on
grid on
title('Junta 3')
plot(time,qvlimmax(3,:),'r')
plot(time,qvlimmin(3,:),'g')
plot(time,control_sig(3,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_3(graus/s)')
legend('qvlim_[3max]', 'qvlim_[3min]','u_3','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(3) 1.1*(180/pi)*qvlim(3)])

subplot(2,3,4)
hold on
grid on
title('Junta 4')
plot(time,qvlimmax(4,:),'r')
plot(time,qvlimmin(4,:),'g')
plot(time,control_sig(4,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_4(graus/s)')
legend('qvlim_[4max]', 'qvlim_[4min]','u_4','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(4) 1.1*(180/pi)*qvlim(4)])

subplot(2,3,5)
hold on
grid on
title('Junta 5')
plot(time,qvlimmax(5,:),'r')
plot(time,qvlimmin(5,:),'g')
plot(time,control_sig(5,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_5(graus/s)')
legend('qvlim_[5max]', 'qvlim_[5min]','u_5','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(5) 1.1*(180/pi)*qvlim(5)])

subplot(2,3,6)
hold on
grid on
title('Junta 6')
plot(time,qvlimmax(6,:),'r')
plot(time,qvlimmin(6,:),'g')
plot(time,control_sig(6,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_6(graus/s)')
legend('qvlim_[6max]', 'qvlim_[6min]','u_6','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(6) 1.1*(180/pi)*qvlim(6)])

hold on

figure(3)
sgtitle('Deslocamento das juntas')
subplot(2,3,1)
hold on
grid on
title('Junta 1')
plot(time,qmax(1,:),'r')
plot(time,qmin(1,:),'g')
plot(time,qplot(1,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_1(graus)')
legend('q_[1max]', 'q_[1min]','q_1','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(1,1) 1.1*(180/pi)*robo1.qlim(1,2)])

subplot(2,3,2)
hold on
grid on
title('Junta 2')
plot(time,qmax(2,:),'r')
plot(time,qmin(2,:),'g')
plot(time,qplot(2,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_2(graus)')
legend('q_[2max]', 'q_[2min]','q_2','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(2,1) 1.1*(180/pi)*robo1.qlim(2,2)])

subplot(2,3,3)
hold on
grid on
title('Junta 3')
plot(time,qmax(3,:),'r')
plot(time,qmin(3,:),'g')
plot(time,qplot(3,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_3(graus)')
legend('qlim_[3max]', 'qlim_[3min]','q_3','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(3,1) 1.1*(180/pi)*robo1.qlim(3,2)])

subplot(2,3,4)
hold on
grid on
title('Junta 4')
plot(time,qmax(4,:),'r')
plot(time,qmin(4,:),'g')
plot(time,qplot(4,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_4(graus)')
legend('qlim_[4max]', 'qlim_[4min]','q_4','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(4,1) 1.1*(180/pi)*robo1.qlim(4,2)])

subplot(2,3,5)
hold on
grid on
title('Junta 5')
plot(time,qmax(5,:),'r')
plot(time,qmin(5,:),'g')
plot(time,qplot(5,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_5(graus)')
legend('qlim_[5max]', 'qlim_[5min]','q_5','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(5,1) 1.1*(180/pi)*robo1.qlim(5,2)])

subplot(2,3,6)
hold on
grid on
title('Junta 6')
plot(time,qmax(6,:),'r')
plot(time,qmin(6,:),'g')
plot(time,qplot(6,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_6(graus)')
legend('qlim_[6max]', 'qlim_[6min]','q_6','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(6,1) 1.1*(180/pi)*robo1.qlim(6,2)])

hold on
figure (4)
sgtitle('Trajetória do efetuador')
subplot(2,2,1)
hold on
grid on
plot3(traj(1,:),traj(2,:),traj(3,:))
view(3)
hold off
legend('Caminho percorrido(m)','Location','Best');

subplot(2,2,2)
hold on
grid on
plot(time,vrpy(1,:))
hold off
xlabel('Iterações')
ylabel('Ângulo(Roll)')
legend('Roll','Location','Best');

subplot(2,2,3)
hold on
grid on
plot(time,vrpy(2,:))
hold off
xlabel('Iterações')
ylabel('Ângulo(Pitch)')
legend('Pitch','Location','Best');

subplot(2,2,4)
hold on
grid on
plot(time,vrpy(3,:))
hold off
xlabel('Iterações')
ylabel('Ângulo(Yaw)')
legend('Yaw','Location','Best');

hold on

figure(5)
sgtitle('Erros de posição')
subplot(2,2,1)
hold on
grid on
plot(time,erro(1,:))
hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em x','Location','Best');


subplot(2,2,2)
hold on
grid on
plot(time,erro(2,:))
hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em y','Location','Best');


subplot(2,2,3)
hold on
grid on
plot(time,erro(3,:))
hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em z','Location','Best');


hold on

figure (6)
sgtitle('Norma do erro')
plot(time,err)