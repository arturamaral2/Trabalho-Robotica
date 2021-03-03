clc, clear, close all

L1 = Link('prismatic', 'alpha', pi/2, 'qlim', [0 .5]);

L2 = Revolute('d', .290, 'alpha', -pi/2, 'qlim', 11/12*[-pi pi]);

L3 = Revolute('a', .270, 'offset', -pi/2, 'qlim',11/18*[-pi pi]);

L4 = Revolute('a', .070, 'alpha', -pi/2, 'qlim', [-(11/18)*pi (7/18)*pi]);

L5 = Revolute('d', .302, 'alpha', pi/2, 'qlim', 8/9*[-pi pi]);

L6 = Revolute('alpha', -pi/2, 'qlim', 2/3*[-pi pi]);

L7 = Revolute('d', 0.072, 'offset', pi, 'qlim', 20/9*[-pi pi]);

robo1 = SerialLink([L1 L2 L3 L4 L5 L6 L7], 'name', 'IRB120')

robo1.base = trotx(-pi/2);

q = [0 0 0 0 0 -pi/2  0];

qvlim = [.500 25*pi/18 25*pi/18 25*pi/18 16*pi/9 16*pi/9 7*pi/3];
%____________________________
sim=remApi('remoteApi');
sim.simxFinish(-1); 
clientID=sim.simxStart('127.0.0.1',19999,true,true,5000,5);
%_________________________________



figure('position',[280 70 560 420])
robo1.plot(q)
view([60,30])

pd = [0.380 0.580 0.600];%3
Rd = SO3();
Rd= Rd.R;
Td = SE3(Rd,pd);
rpyd = rotm2eul(Rd);

lambda = 0.8;
epsilon = 2e-2;
e = inf(6,1);  

t0 = tic;



w=0;


rpyd = [0,0,0]

while (norm(e) > epsilon) % Critério de parada
    T = robo1.fkine(q);
    J = robo1.jacob0(q, 'rpy'); % Jacobiana analítica

    p = transl(T);
    p_til = pd - p;
   
    R = SO3(T);
    R = R.R; % Extrai rotação do efetuador
    rpy = rotm2eul(R);
   
    rpy_til = rpyd-rpy;
   
    e = [p_til';rpy_til']; % Vetor de erro
       
    u = pinv(J)*lambda*e; % Lei de controle
    tf= toc(t0);
    t0= tic;
   
    for k = 1:7
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
   
    tf_robo =  @(t,qi) [u(1);u(2);u(3);u(4);u(5);u(6);u(7)];
    [~, qi] = ode45(tf_robo, 0:tf:tf, zeros(1,7));%4
   
    control_sig(:,w+1) = 180*u/pi;
    control_sig(1,w+1) = u(1);
    err(w+1)= norm(e);
    traj(:,w+1) = p;
    vrpy(:,w+1) = 180*rpy/pi;
    erro(:,w+1) = e;

    qplot(1,w+1) = q(1);
    qmax(1,w+1) = robo1.qlim(1,2);
    qmin(1,w+1) = robo1.qlim(1,1);
    qvlimmax(1,w+1) = qvlim(1);
    qvlimmin(1,w+1) = qvlim(1);  
   
    for k = 2:7
        qplot(k,w+1) = 180*q(k)/pi;
        qmax(k,w+1) = 180*robo1.qlim(k,2)/pi;
        qmin(k,w+1) = 180*robo1.qlim(k,1)/pi;
        qvlimmax(k,w+1) = 180*qvlim(k)/pi;
        qvlimmin(k,w+1) = -180*qvlim(k)/pi;
    end
   
    time(w+1) = w;   
    q = q+qi(end,:);           
    robo1.plot(q);
   
if (clientID>-1)
    disp('Connected to remote API server');
h = [0,0,0,0,0,0,0];
    [r,h(2)]=sim.simxGetObjectHandle(clientID,'joint_1',sim.simx_opmode_blocking);
    [r,h(3)]=sim.simxGetObjectHandle(clientID,'joint_2',sim.simx_opmode_blocking);
    [r,h(4)]=sim.simxGetObjectHandle(clientID,'joint_3',sim.simx_opmode_blocking);
    [r,h(5)]=sim.simxGetObjectHandle(clientID,'joint_4',sim.simx_opmode_blocking);
    [r,h(6)]=sim.simxGetObjectHandle(clientID,'joint_5',sim.simx_opmode_blocking);
    [r,h(7)]=sim.simxGetObjectHandle(clientID,'joint_6',sim.simx_opmode_blocking);
    [r,h(1)]=sim.simxGetObjectHandle(clientID,'joint_7',sim.simx_opmode_blocking);%2

    for i=1:7
    sim.simxSetJointTargetPosition(clientID,h(i),q(i),sim.simx_opmode_streaming);
    end      
else
    disp('Failed connecting to remote API server');  
end
    sim.delete(); % call the destructor!
   
    disp('Program ended');                 
    w=w+1;

end
           
hold off

%% Plot sinal de controle e norma do erro
hold on
figure(2)
sgtitle('Sinais de controle(velocidade das juntas)')
subplot(2,4,1)
hold on
grid on
title('Trilho')
plot(time,qvlimmax(1,:),'r')
plot(time,qvlimmin(1,:),'g')
plot(time,control_sig(1,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_1(m/s)')
legend('qvlim_[1max]', 'qvlim_[1min]','u_1','Location','Best');
axis ([0 length(time)-1 -qvlim(1) qvlim(1)])

subplot(2,4,2)
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

subplot(2,4,3)
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

subplot(2,4,4)
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

subplot(2,4,5)
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

subplot(2,4,6)
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

subplot(2,4,7)
hold on
grid on
title('Junta 7')
plot(time,qvlimmax(7,:),'r')
plot(time,qvlimmin(7,:),'g')
plot(time,control_sig(7,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_7(graus/s)')
legend('qvlim_[7max]', 'qvlim_[7min]','u_7','Location','Best');
axis ([0 length(time)-1 -1.1*(180/pi)*qvlim(7) 1.1*(180/pi)*qvlim(7)])

hold on

figure(3)
sgtitle('Deslocamento das juntas')
subplot(2,4,1)
hold on
grid on
title('Trilho')
plot(time,qmax(1,:),'r')
plot(time,qmin(1,:),'g')
plot(time,qplot(1,:),'b')
hold off
xlabel('Iterações')
ylabel('qlim_1(m)')
legend('q_[1max]', 'q_[1min]','q_1','Location','Best');
axis ([0 length(time)-1 robo1.qlim(1,1) robo1.qlim(1,2)])

subplot(2,4,2)
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

subplot(2,4,3)
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

subplot(2,4,4)
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

subplot(2,4,5)
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

subplot(2,4,6)
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

subplot(2,4,7)
hold on
grid on
title('Junta 7')
plot(time,qmax(7,:),'r')
plot(time,qmin(7,:),'g')
plot(time,qplot(7,:),'b')
hold off
xlabel('Iterações')
ylabel('qvlim_7(graus)')
legend('qlim_[7max]', 'qlim_[7min]','q_7','Location','Best');
axis ([0 length(time)-1 1.1*(180/pi)*robo1.qlim(7,1) 1.1*(180/pi)*robo1.qlim(7,2)])

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
sgtitle('Erros de posição e orientação')
subplot(2,3,1)
hold on
grid on
plot(time,erro(1,:))
hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em x','Location','Best');


subplot(2,3,2)
hold on
grid on
plot(time,erro(2,:))

hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em y','Location','Best');


subplot(2,3,3)
hold on
grid on
plot(time,erro(3,:))
hold off
xlabel('Iterações')
ylabel('Erro(m)')
legend('Erro em z','Location','Best');

subplot(2,3,4)
hold on
grid on
plot(time,erro(4,:))
hold off
xlabel('Iterações')
ylabel('Erro(°)')
legend('Erro em roll','Location','Best');

subplot(2,3,5)
hold on
grid on
plot(time,erro(5,:))
hold off
xlabel('Iterações')
ylabel('Erro(°)')
legend('Erro em pitch','Location','Best');

subplot(2,3,6)
hold on
grid on

plot(time,erro(6,:))
hold off
xlabel('Iterações')
ylabel('Erro(°)')
legend('Erro em yaw','Location','Best');

hold on

figure (6)
sgtitle('Norma do erro')
plot(time,err)