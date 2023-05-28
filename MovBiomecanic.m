clear 
close all
clc

%% Ingreso el archivo y lo leo
archivo = uigetfile(' .csv');
archivo = archivo(1:end-4);
M = csv2strc(archivo);

%% Obtención de marcadores, medidas antropometricas, tiempo y los Frames
Markers = M.Sets.Subject.Raw;
Antrop  = M.Anthropometry.Subject;
NFrames = M.Frames;
Time = M.Time;

%% Matriz 3D Marcadores

Nomb_Marcadores = fieldnames(Markers);
[NMarkers,~]= size(Nomb_Marcadores);
Marcadores = zeros(NFrames,NMarkers,3);

for i=1:NMarkers
    NMarcador = Nomb_Marcadores{i,1};
    Marcadores(:,i,1)= Markers.(sprintf('%s', NMarcador))(:,1);
    Marcadores(:,i,2)= Markers.(sprintf('%s', NMarcador))(:,2);
    Marcadores(:,i,3)= Markers.(sprintf('%s', NMarcador))(:,3);
end

%% Varias Antroprometricas
% Recordar que: A.A5,A.A6 y A.A9, A.A10 No son necesarias
A.A0= Antrop.TotalHeight.Value/1000;
A.A1= Antrop.TotalBodyMass.Value;
A.A2= Antrop.ASISBreadth.Value/1000;
A.A3= Antrop.RightThighLength.Value/1000;
A.A4= Antrop.LeftThighLength.Value/1000;
A.A7=Antrop.RightCalfLength.Value/1000;
A.A8=Antrop.LeftCalfLength.Value/1000;
A.A11=Antrop.RightKneeDiameter.Value/1000;
A.A12=Antrop.LeftKneeDiameter.Value/1000;
A.A13=Antrop.RightFootLength.Value/1000;
A.A14=Antrop.LeftFootLength.Value/1000;
A.A15=Antrop.RightMalleolusHeight.Value/1000;
A.A16=Antrop.LeftMalleolusHeight.Value/1000;
A.A17=Antrop.RightMalleolusWidth.Value/1000;
A.A18=Antrop.LeftMalleolusWidth.Value/1000;
A.A19=Antrop.RightFootBreadth.Value/1000;
A.A20=Antrop.LeftFootBreadth.Value/1000;


%% Centros

for  v=1: NFrames
    P1=Markers.R_Met(v,:);
    P2=Markers.R_Heel(v,:);
    P3=Markers.R_Mall(v,:);
    P4=Markers.R_Bar2(v,:);
    P5=Markers.R_Knee1(v,:);
    P6=Markers.R_Bar1(v,:);
    P7=Markers.R_Asis(v,:);
    P8=Markers.L_Met(v,:);
    P9=Markers.L_Heel(v,:);
    P10=Markers.L_Mall(v,:);
    P11=Markers.L_Bar2(v,:);
    P12=Markers.L_Knee1(v,:);
    P13=Markers.L_Bar1(v,:);
    P14=Markers.L_Asis(v,:);
    P15=Markers.Sacrum(v,:);
    
    
    
    %Cadera 
    
    vPelvis=(P14-P7)/norm(P14-P7);
    wPelvis=cross(P7-P15,P14-P15)/norm(cross(P7-P15,P14-P15));
    uPelvis=cross(vPelvis,wPelvis);
    
    PRHip= P15+0.598*A.A2*uPelvis-0.344*A.A2*vPelvis-0.29*A.A2*wPelvis;
    PLHip= P15+0.598*A.A2*uPelvis+0.344*A.A2*vPelvis-0.29*A.A2*wPelvis;
    
    Centers.R_Hip(v,:)=PRHip;
    Centers.L_Hip(v,:)=PLHip;
    
    %Rodilla Derecha 
    
    vRodilla=(P3-P5)/norm(P3-P5);
    uRodilla=cross(P4-P5,P3-P5)/norm(cross(P4-P5,P3-P5));
    wRodilla=cross(uRodilla,vRodilla);
   
    PRKnee = P5 + 0.5*A.A11*wRodilla;
    
    Centers.R_Knee(v,:)=PRKnee;
    
    %Rodilla Izquierda
   
    vRodillaI=(P10-P12)/norm(P10-P12);
    uRodillaI=cross(P10-P12,P11-P12)/norm(cross(P10-P12,P11-P12));
    wRodillaI=cross(uRodillaI,vRodillaI);
   
    PLKnee = P12 - 0.5*A.A12*wRodillaI;
    
    Centers.L_Knee(v,:)=PLKnee;
    
    %Pie Derecho
    
    uPie=(P1-P2)/norm(P1-P2);
    wPie=cross(P1-P3,P2-P3)/norm(cross(P1-P3,P2-P3));
    vPie=cross(wPie,uPie);
    
    PRAnkle=P3+0.16*A.A13*uPie+0.392*A.A15*vPie+0.478*A.A17*wPie;
    PRToe=P3+0.742*A.A13*uPie+1.074*A.A15*vPie-0.187*A.A19*wPie;
    
    Centers.R_Ankle(v,:)=PRAnkle;
    Centers.R_Toe(v,:)=PRToe;
    
    %Pie Izquierdo
    
    uPieI=(P8-P9)/norm(P8-P9);
    wPieI=cross(P8-P10,P9-P10)/norm(cross(P8-P10,P9-P10));
    vPieI=cross(wPieI,uPieI);
    
    PLAnkle=P10+0.16*A.A14*uPieI+0.392*A.A16*vPieI-0.478*A.A18*wPieI;
    PLToe=P10+0.742*A.A14*uPieI+1.074*A.A16*vPieI+0.187*A.A20*wPieI;
    
    Centers.L_Ankle(v,:)=PLAnkle;
    Centers.L_Toe(v,:)=PLToe;
    
    
    PRHeel=P2;
    PLHeel=P9;
    
    %% Centros de Masa
    
    %1
    
    PRThighCG = PRHip+0.39*(PRKnee-PRHip);
    MidP= ((P7-P14)/2)+P14;
    PPelvisCG=((P15-MidP)/2)+MidP;
    
    Segmentos.Pelvis.CG(v,:)=PPelvisCG;
    Segmentos.R_Thigh.CG(v,:)=PRThighCG;
    
    %2
    
    PLThighCGL= PLHip+0.39*(PLKnee-PLHip);
    MidPL= ((P14-P7)/2)+P7;
    PPelvisCGL=((P15-MidPL)/2)+MidPL;
    
    Segmentos.Pelvis.CG(v,:)=PPelvisCGL;
    Segmentos.L_Thigh.CG(v,:)=PLThighCGL;
    
    %3
    
    PRCalfCG=PRKnee+0.42*(PRAnkle-PRKnee);
    Segmentos.R_Knee.CG(v,:)=PRCalfCG;
    
    %4
    
    PLCalfCG=PLKnee+0.42*(PLAnkle-PLKnee);
    Segmentos.L_Knee.CG(v,:)=PLCalfCG;
    
    %5
    
    PRFootCG=PRHeel+0.44*(PRToe-PRHeel);
    Segmentos.R_Heel.CG(v,:)=PRFootCG;
    
    %6
    
    PLFootCG=PLHeel+0.44*(PLToe-PLHeel);
    Segmentos.L_Heel.CG(v,:)=PLFootCG;
    
    %% Segmentos i-j-k
    
    %Pelvis
    iPelvis=wPelvis;
    jPelvis=uPelvis;
    kPelvis=vPelvis;
    
    Segmentos.Pelvis.i(v,:)=iPelvis;
    Segmentos.Pelvis.j(v,:)=jPelvis;
    Segmentos.Pelvis.k(v,:)=kPelvis;
    
    %Muslo Derecho
    i1=(PRHip-PRKnee)/norm(PRHip-PRKnee);
    j1= cross(P6-PRHip,PRKnee-PRHip)/norm(cross(P6-PRHip,PRKnee-PRHip));
    k1= cross(i1,j1);
    
    Segmentos.R_Thigh.i(v,:)= i1;
    Segmentos.R_Thigh.j(v,:)= j1;
    Segmentos.R_Thigh.k(v,:)= k1;
    
    %Muslo Izquierdo
    i2=(PLHip-PLKnee)/norm(PLHip-PLKnee);
    j2= cross(PLKnee-PLHip,P13-PLHip)/norm(cross(PLKnee-PLHip,P13-PLHip));
    k2= cross(i2,j2);
    
    Segmentos.L_Thigh.i(v,:)= i2;
    Segmentos.L_Thigh.j(v,:)= j2;
    Segmentos.L_Thigh.k(v,:)= k2;
    
    %Pierna Derecha
    i3=(PRKnee-PRAnkle)/norm(PRKnee-PRAnkle);
    j3=cross(P4-PRKnee,PRAnkle-PRKnee)/norm(cross(P4-PRKnee,PRAnkle-PRKnee));
    k3=cross(i3,j3);
    
    Segmentos.R_Ankle.i(v,:)= i3;
    Segmentos.R_Ankle.j(v,:)= j3;
    Segmentos.R_Ankle.k(v,:)= k3;
    
    %Pierna Izquierda
    i4=(PLKnee-PLAnkle)/norm(PLKnee-PLAnkle);
    j4=cross(PLAnkle-PLKnee,P11-PLKnee)/norm(cross(PLAnkle-PLKnee,P11-PLKnee));
    k4=cross(i4,j4);
    
    Segmentos.L_Ankle.i(v,:)= i4;
    Segmentos.L_Ankle.j(v,:)= j4;
    Segmentos.L_Ankle.k(v,:)= k4;
    
    %Pie derecho
    i5=(P2-PRToe)/norm(P2-PRToe);
    k5=cross(PRAnkle-P2,PRToe-P2)/norm(cross(PRAnkle-P2,PRToe-P2));
    j5= cross(k5,i5);
    
    Segmentos.R_Toe.i(v,:)= i5;
    Segmentos.R_Toe.j(v,:)= j5;
    Segmentos.R_Toe.k(v,:)= k5;
    
    %Pie Izquierdo 
    i6=(P9-PLToe)/norm(P9-PLToe);
    k6=cross(PLAnkle-P9,PLToe-P9)/norm(cross(PLAnkle-P9,PLToe-P9));
    j6= cross(k6,i6);
    
    Segmentos.L_Toe.i(v,:)= i6;
    Segmentos.L_Toe.j(v,:)= j6;
    Segmentos.L_Toe.k(v,:)= k6;
    
    
    %% Sistemas de Coordenadas
    %Angulos Articulares.
    %Cadera Derecha
    iRHip=cross(kPelvis,i1)/norm(cross(kPelvis,i1));
    
    aRHip=acosd(dot(iRHip,jPelvis))*(dot(iRHip,iPelvis)/abs(dot(iRHip,iPelvis)));
    bRHip= asind(dot(kPelvis,i1));
    yRHip=-acosd(dot(iRHip,j1))*(dot(iRHip,k1)/abs(dot(iRHip,k1)));
    
    Angulos.R_Hip.i(v,:)=iRHip;
    Angulos.R_Hip.a(v,1)=aRHip;
    Angulos.R_Hip.b(v,1)=bRHip;
    Angulos.R_Hip.y(v,1)=yRHip;
    
    %Cadera Izquierda
    iLHip=cross(kPelvis,i2)/norm(cross(kPelvis,i2));
    
    aLHip=acosd(dot(iLHip,jPelvis))*(dot(iLHip,iPelvis)/abs(dot(iLHip,iPelvis)));
    bLHip= -asind(dot(kPelvis,i2));
    yLHip=acosd(dot(iLHip,j2))*(dot(iLHip,k2)/abs(dot(iLHip,k2)));
  
    Angulos.L_Hip.i(v,:)=iLHip;
    Angulos.L_Hip.a(v,1)=aLHip;
    Angulos.L_Hip.b(v,1)=bLHip;
    Angulos.L_Hip.y(v,1)=yLHip;
    
    %Rodilla Derecha
    iRKnee=cross(k1,i3)/norm(cross(k1,i3));
   
    aRKnee=-acosd(dot(iRKnee,j1))*(dot(iRKnee,i1)/abs(dot(iRKnee,i1)));
    bRKnee= asind(dot(k1,i3));
    yRKnee=-acosd(dot(iRKnee,j3))*(dot(iRKnee,k3)/abs(dot(iRKnee,k3)));
    
    Angulos.R_knee.i(v,:)=iRKnee;
    Angulos.R_knee.a(v,1)=aRKnee;
    Angulos.R_knee.b(v,1)=bRKnee;
    Angulos.R_knee.y(v,1)=yRKnee;
    
    %Rodilla Izquierda
    iLKnee=cross(k2,i4)/norm(cross(k2,i4));
   
    aLKnee=-acosd(dot(iLKnee,j2))*(dot(iLKnee,i2)/abs(dot(iLKnee,i2)));
    bLKnee= -asind(dot(k2,i4));
    yLKnee=acosd(dot(iLKnee,j4))*(dot(iLKnee,k4)/abs(dot(iLKnee,k4)));
    
    Angulos.L_knee.i(v,:)=iLKnee;
    Angulos.L_knee.a(v,1)=aLKnee;
    Angulos.L_knee.b(v,1)=bLKnee;
    Angulos.L_knee.y(v,1)=yLKnee;
    
    %Tobillo Derecho
    iRAnkle=cross(k3,i5)/norm(cross(k3,i5));
    
    aRAnkle=-asind(dot(iRAnkle,j3));
    bRAnkle=asind(dot(k3,i5));
    yRAnkle=asind(dot(iRAnkle,k5));
    
    Angulos.R_Ankle.i(v,:)=iRAnkle;
    Angulos.R_Ankle.a(v,1)=aRAnkle;
    Angulos.R_Ankle.b(v,1)=bRAnkle;
    Angulos.R_Ankle.y(v,1)=yRAnkle;
    
    %Tobillo Izquierdo
    iLAnkle=cross(k4,i6)/norm(cross(k4,i6));
    
    aLAnkle=-asind(dot(iLAnkle,j4));
    bLAnkle=-asind(dot(k4,i6));
    yLAnkle=-asind(dot(iLAnkle,k6));
    
    Angulos.L_Ankle.i(v,:)=iLAnkle;
    Angulos.L_Ankle.a(v,1)=aLAnkle;
    Angulos.L_Ankle.b(v,1)=bLAnkle;
    Angulos.L_Ankle.y(v,1)=yLAnkle; 
 end

xmin = min(min(Marcadores(:,:,1)));
xmax = max(max(Marcadores(:,:,1)));
ymin = min(min(Marcadores(:,:,2)));
ymax = max(max(Marcadores(:,:,2)));
% zmin = min(min(Marcadores(:,:,3)));
zmax = max(max(Marcadores(:,:,3))); 
 
Ejes = [xmin-0.2 xmax+0.2 ymin-0.5 ymax+0.1 0 zmax+0.2];

%% Matriz 3D Centros
Center{1,1}='R_Toe';
Center{2,1}='R_Ankle';
Center{3,1}='R_Knee';
Center{4,1}='R_Hip';
Center{5,1}='L_Hip';
Center{6,1}='L_Knee';
Center{7,1}='L_Ankle';
Center{8,1}='L_Toe';
NombresCenters = length(Center);
Centros = zeros(NFrames,NombresCenters,3);

for i=1:NombresCenters
   
    NCentro = Center{i,1};
    Centros(:,i,1)= Centers.(sprintf('%s', NCentro))(:,1);
    Centros(:,i,2)= Centers.(sprintf('%s', NCentro))(:,2);
    Centros(:,i,3)= Centers.(sprintf('%s', NCentro))(:,3);
end
 
 %% GRAFICO DE LOS CENTROS ARTICULARES CON LOS EJES
tam_seg = 0.15;

%  
% for f=1:30:NFrames   
%     plot3(Centros(f,:,1), Centros(f,:,2), Centros(f,:,3),'*-k');
%     frameref([0 0 1], [0 1 0], [1 0 0], [0 0 0], 'k', 0.20, 'G');
%    
%      
%     xlabel('EJE X');
%     ylabel('EJE Y');
%     zlabel('EJE Z');
%     axis equal;
%     axis(Ejes)
%     grid on;
%     drawnow;
%     hold off
% end
figure
for f=1:30:NFrames
    subplot(1,2,1)
    plot3(Centros(f,:,1), Centros(f,:,2), Centros(f,:,3),'*-k');
     xlabel('EJE X');
     ylabel('EJE Y');
     zlabel('EJE Z');
     axis equal;
     axis(Ejes)
     grid on;
     drawnow;
     hold on
     hold off
     
    subplot(1,2,2)
    plot3(Centros(f,:,1), Centros(f,:,2), Centros(f,:,3),'k');
    hold on
    
    frameref([0 0 1], [0 1 0], [1 0 0], [0 0 0], 'k', 0.20, 'G');
    frameref(Segmentos.Pelvis.i(f,:),Segmentos.Pelvis.j(f,:),Segmentos.Pelvis.k(f,:),Segmentos.Pelvis.CG(f,:),'b',tam_seg,'L');
    
    frameref(Segmentos.R_Thigh.i(f,:),Segmentos.R_Thigh.j(f,:),Segmentos.R_Thigh.k(f,:),Segmentos.R_Thigh.CG(f,:),'g',tam_seg,'L');
    frameref(Segmentos.L_Thigh.i(f,:),Segmentos.L_Thigh.j(f,:),Segmentos.L_Thigh.k(f,:),Segmentos.L_Thigh.CG(f,:),'r',tam_seg,'L');
   
    frameref(Segmentos.R_Ankle.i(f,:),Segmentos.R_Ankle.j(f,:),Segmentos.R_Ankle.k(f,:),Segmentos.R_Knee.CG(f,:),'g',tam_seg,'L');
    frameref(Segmentos.L_Ankle.i(f,:),Segmentos.L_Ankle.j(f,:),Segmentos.L_Ankle.k(f,:),Segmentos.L_Knee.CG(f,:),'r',tam_seg,'L');
   
    frameref(Segmentos.R_Toe.i(f,:),Segmentos.R_Toe.j(f,:),Segmentos.R_Toe.k(f,:),Segmentos.R_Heel.CG(f,:),'g',tam_seg,'L');
    frameref(Segmentos.L_Toe.i(f,:),Segmentos.L_Toe.j(f,:),Segmentos.L_Toe.k(f,:),Segmentos.L_Heel.CG(f,:),'r',tam_seg,'L');
    
    xlabel('EJE X');
    ylabel('EJE Y');
    zlabel('EJE Z');
    axis equal;
    axis(Ejes)
    grid on;
    drawnow;
    hold off
      % Guardar el plot actual como una imagen PNG
    filename = sprintf('plot_frame_%d.png', f);
    saveas(gcf, filename, 'png');
    close(gcf);
end




ts = 100/NFrames;%Transformacion para que en las graficas apreciera el ciclo de marcha entre %0 y 100
T=ts:ts:NFrames*ts; %Arreglo para el avance, arranca en el frame 1 que se convierte en el 0, y va hasta el 100

%%%%%%%%%%%%

figure (10)

subplot(3,3,1)
plot (T,Angulos.R_Hip.a,'g')
hold on
grid on
title ('Sagital')
ylabel([{'Cadera'};{'Ext(-)    Flx(+)'}]);
plot (T,Angulos.L_Hip.a,'r')


subplot(3,3,2);
plot (T,Angulos.R_Hip.b,'g')
hold on
grid on
title ('Frontal')
ylabel('Aduc(-)    Abduc(+)');
plot (T,Angulos.L_Hip.b,'r')

    
subplot(3,3,3);
plot (T,Angulos.R_Hip.y,'g')
hold on
grid on
title ('Transversal')
ylabel('R Ext(-)    R Int(+)');
plot (T,Angulos.L_Hip.y,'r')



subplot(3,3,4)
plot (T,Angulos.R_knee.a,'g')
hold on
grid on
title ('Sagital')
ylabel([{'Rodilla'};{'Ext(-)    Flx(+)'}]);
plot (T,Angulos.L_knee.a,'r')


subplot(3,3,5);
plot (T,Angulos.R_knee.b,'g')
hold on
grid on
title ('Frontal')
ylabel('Varo(-)    Valgo(+');
plot (T,Angulos.L_Hip.b,'r')

    
subplot(3,3,6);
plot (T,Angulos.R_knee.y,'g')
hold on
grid on
title ('Transversal')
ylabel('R Ext(-)    R Int(+)');
plot (T,Angulos.L_knee.y,'r')



subplot(3,3,7);
plot (T,Angulos.R_Ankle.a,'g')
hold on
grid on
title ('Sagital')
ylabel([{'Tobillo'};{'Planti(-)    Dorsi(+)'}]);
plot (T,Angulos.L_Ankle.a,'r')

subplot(3,3,8);
plot (T,Angulos.R_Ankle.y,'g')
hold on
grid on
title ('Frontal')
ylabel('Ever(-)    Inv(+)');
plot (T,Angulos.L_Ankle.y,'r')

subplot(3,3,9);
plot (T,Angulos.R_Ankle.b,'g')
hold on
grid on
title ('Transversal')
ylabel('P Int(-)    P Ext(+)');
plot (T,Angulos.L_Ankle.y,'r')