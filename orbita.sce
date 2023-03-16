function RK=rka(state, tau, erro);
// calcula um passo da integração da órbita por 
// Runge-Kutta-4 "adaptive step", redimensionando, 
// o intervalo de tempo, tau;
S1=.9; S2=2;
  ////  atualização de state e tau  ////
 h=tau/2;hh=h/2;
   ////  os dois passos pequenos ////
 F1=fun(state);
 F2=fun(state+hh*F1);
 F3=fun(state+hh*F2);
 F4=fun(state+h*F3);
 F=(F1+2*F2+2*F3+F4)/6;
 stmed=state+h*F;
 F1=fun(stmed);
 F2=fun(stmed+hh*F1);
 F3=fun(stmed+hh*F2);
 F4=fun(stmed+h*F3);
 F=(F1+2*F2+2*F3+F4)/6;
 stfin=stmed+h*F;
    ////  passo único grande  ////
 F1=fun(state);
 F2=fun(state+h*F1);
 F3=fun(state+h*F2);
 F4=fun(state+tau*F3);
 F=(F1+2*F2+2*F3+F4)/6;
 stg=state+tau*F;
    ////  computa o erro corrente////
 dif2=(stg(1)-stfin(1))^2+(stg(2)-stfin(2))^2;
 ds2=(stfin(1)-state(1))^2+(stfin(2)-state(2))^2;
 errc=sqrt(dif2/(ds2));
   ////  decide novo tau ////
 if(errc > erro);    
    tau=S1*tau*(erro/(errc+%eps))^(.2);
 elseif(errc < erro/S2)
    tau=S1*tau*(erro/(errc+%eps))^(.2);
    tau=min(tau, S2*tau);
 end
RK=[stfin, tau, errc];



// Programa: orbita.sce; Calcula a órbita de um cometa 
// usando "adaptive Runge-Kutta", redimensionando tau 
// após cada passo; dá opção de continuação do cálculo.  
// Traça também a órbita exata;  
// Ex: [x0=1000, v0=.002, dt=5, erro=1.e-6, N=200] 

par=input('deh [distancia, vel.tang., dt, erro, N],\n');
x0=par(1); v0=par(2); dt=par(3);
erro=par(4); N=par(5);
// órbita exata:
L=x0*v0; E0=L^2/(2*x0^2)-1/x0; //mom. angular e energia
a=1/(2*E0); // semi-eixo maior 
ec2=1+2*E0*L^2; ec=sqrt(ec2); // ecentricidade
pi=%pi;
teta=-pi:pi/101:pi; I=ones(1,length(teta));
re=a*(1-ec2)*I./(1+ec*cos(teta)); // raio exato
xe=re.*cos(teta); ye=re.*sin(teta);
// cálculo numérico por RK4:
estado=[x0, 0, 0, v0];  // a posição e a vel. inicial 
x=[];y=[];vx=[];vy=[];t=[];tau=[]; erroc=[];
x(1)=x0;y(1)=0; vx(1)=0; vy(1)=v0; n=1; 
tau(1)=dt; t(1)=0;
function F=fun(R);
    r3=(sqrt(R(1)^2 + R(2)^2))^3;
    F=[R(3), R(4), -R(1)/r3, -R(2)/r3];
endfunction

// introduz rka:
while(N >= 1)
  x=[x,zeros(1,N)];
  y=[y,zeros(1,N)];
  vx=[vx,zeros(1,N)];
  vy=[vy,zeros(1,N)];
  tau=[tau,zeros(1,N)];
  t=[t,zeros(1,N)]; 
  erroc=[erroc,zeros(1,N)];
for it=n+1:n+N;
 R=rka(estado, dt, erro);
 estado=R(1:4); 
 x(it)=R(1); y(it)=R(2);
 vx(it)=R(3); vy(it)=R(4); 
 tau(it)=R(5);
 if(it>1); t(it)=t(it-1)+tau(it); end
 erroc(it)=R(6); // erro calculado
 dt=R(5);
end;
xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
xset('window',0); clf; isoview(xmin,xmax,ymin,ymax);
plot(x,y,'k',0,0,'ro'); 
legend('orbita calculada')
n=n+N;
N=input('quantos passos mais? 0:encerra, \n');
end
// Energia
Ec=(vx.^2 + vy.^2)/2; raio=sqrt(x.^2+y.^2); 
V=-ones(raio)./raio; E=Ec+V; 
xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
xset('window',0); clf; isoview(xmin,xmax,ymin,ymax);
plot(xe,ye,x,y,0,0,'r.');
legend('exata','rka')
xset('window',1); clf;
plot(t,Ec,t,V,t,E);
legend('En. cinetica', 'En. potencial', 'En. total')
xset('window',2); clf;
plot(t,tau);
legend('tau(t)') 
xset('window',3); clf;
plot(t,erroc);
legend('erro computado',4) 
