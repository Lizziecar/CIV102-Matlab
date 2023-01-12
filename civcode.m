%%% Intial Parameters for both designs and cases:
n = 1250;                  % Number of locations to evaluate bridge failure 
L = 1250;                  % Length of bridge 
  
x = linspace(0, L, n);     % Define x coordinate 

locationA = 1;
locationB = 1060;

%Material properties 
SigT = 30; 
SigC = 6; 
E    = 4000; 
TauU = 4; 
TauG = 2; 
mu   = 0.2; 

%%% Design Zero: Train
SFD_TrainLoad = zeros(1, n);      % Initialize SFD(x) 
P = 400/6;
TLocation = 0;

%diaphragms
a = zeros(1,length(x));
a(1:550) = 550;
a(551:1060) = 510;
a(1061:1250) = 190;

%Cross Section
%% 2. Define cross-sections
% Design 0
dx = n/L;
xc = [1 dx*L]; % x locations of changes in cross section
bVector = [100 1.27 1.27 10 10 80];
hVector = [1.27 72.46 72.46 1.27 1.27 1.27];
yVector = [74.365 37.5 37.5 73.1 73.1 0.635]; % distances between local centroid of subcomponents of sections from bottom
 
I = zeros(1, length(x));
yBar = zeros(1, length(x));
[I, yBar] = SectionProperties( I, yBar, xc(1), xc(2), bVector, hVector, yVector);
 
cutoffAreaCent = [ 77.46*1.27, 1.27*yBar(1), 1.27*yBar(1)];
yCentCutoff = [ yBar(1)-0.635, yBar(1)/2, yBar(1)/2];
 
cutoffAreaGlue = [127];
yGlue = [75-yBar(1)-0.635]; % centroid of the cutoff area
 
Qglue = zeros(1, length(x));
Qmax = zeros(1, length(x));
[Qmax, Qglue] = calculateQ( Qmax, Qglue, xc(1), xc(2), cutoffAreaCent, cutoffAreaGlue, yCentCutoff, yGlue);


%Assuming I and Q
%I = zeros(1,length(x));
%I(:) = 415.684*10^3;
%y = zeros(1,length(x));
%y(:) = 41.701;
b = zeros(1,length(x));
b(:) = 2*1.27;
bGlue = zeros(1,length(x));
bGlue(:) = 10;
%Qcent = zeros(1,length(x));
%Qcent(:) = 6248.5;
%Qglue = zeros(1,length(x));
%Qglue(:) = 4148.2;
t = zeros(1,length(x));
t(:) = 1.27;
height = zeros(1,length(x));
height(:) = 75;

%% Train Load:
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176 + 164, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176 + 164 + 176, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176 + 164 + 176 + 164, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176 + 164 + 176 + 164 + 176, P, x, SFD_TrainLoad, locationA, locationB);
[SFD_TrainLoad, BMD_TrainLoad] = ApplyPL(TLocation + 52 + 176 + 164 + 176 + 164 + 176 + 52, P, x, SFD_TrainLoad, locationA, locationB);

%Find applied stresses
%Shear Stress demand
ShearStress = VStress(I, b, Qmax, SFD_TrainLoad);
%ShearStressGlue

%Shear Stress Capacity
ShearBuckle = VBuck(I, t, b, height, E, mu, yBar, a, Qmax);

%Bending Stress
%Bending Stress Demand
BendingStressT =  MStressT( I,yBar, height, BMD_TrainLoad);
BendingStressC =  MStressC( I,yBar, height, BMD_TrainLoad);

%Shear Stress Capacity
BendingBuck1 = MBuck( I, t, (100-80)/2, E, mu, 2);
BendingBuck2 = MBuck( I, t, 80-2*1.27, E, mu, 1 );
BendingBuck3 = MBuck( I, t, 75-41.701-1.27, E, mu, 3 );
BendingBuck4 = MBuck( I, t, (100-80)/2, E, mu, 2 );
BendingBuck5 = MBuck( I, t, 75-41.701-1.27, E, mu, 3 );
BendingBuck6 = MBuck( I, t, 41.701-1.27, E, mu, 3 );
BendingBuck7 = MBuck( I, t, 80-2*1.27, E, mu, 1 );
BendingBuck8 = MBuck( I, t, 41.701-1.27, E, mu, 1 );


FOSBridge = FOSCalc(ShearStress, BendingStressT, BendingStressC, SigT, SigC, BendingBuck1, BendingBuck2, BendingBuck3, BendingBuck4, BendingBuck5, BendingBuck6, BendingBuck7, BendingBuck8, TauU, ShearBuckle);

%plot(x,BMD_TrainLoad)
%set(gca,'YDir','reverse

VisualizeTL(x,SFD_TrainLoad, BMD_TrainLoad, ShearStress, TauU, ShearBuckle, BendingStressT, BendingStressC, SigT, SigC, BendingBuck1, BendingBuck2, BendingBuck3, BendingBuck4, BendingBuck5, BendingBuck6, BendingBuck7, BendingBuck8, FOSBridge)



%Design 0 Point Loads
SFD_PL = zeros(1, n);      % Initialize SFD(x) 
P = 1; %P = 1 * whatever value, so I can set p to 1

%diaphragms
a = zeros(1,length(x));
a(1:550) = 550;
a(551:1060) = 510;
a(1061:1250) = 190;

%Cross Section
%% 2. Define cross-sections
% Design 0
dx = n/L;
xc = [1 dx*L]; % x locations of changes in cross section
bVector = [100 1.27 1.27 10 10 80];
hVector = [1.27 72.46 72.46 1.27 1.27 1.27];
yVector = [74.365 37.5 37.5 73.1 73.1 0.635]; % distances between local centroid of subcomponents of sections from bottom
 
I = zeros(1, length(x));
yBar = zeros(1, length(x));
[I, yBar] = SectionProperties( I, yBar, xc(1), xc(2), bVector, hVector, yVector);
 
cutoffAreaCent = [ 77.46*1.27, 1.27*yBar(1), 1.27*yBar(1)];
yCentCutoff = [ yBar(1)-0.635, yBar(1)/2, yBar(1)/2];
 
cutoffAreaGlue = [127];
yGlue = [75-yBar(1)-0.635]; % centroid of the cutoff area
 
Qglue = zeros(1, length(x));
Qmax = zeros(1, length(x));
[Qmax, Qglue] = calculateQ( Qmax, Qglue, xc(1), xc(2), cutoffAreaCent, cutoffAreaGlue, yCentCutoff, yGlue);

%Assuming I and Q
%I = zeros(1,length(x));
%I(:) = 415.684*10^3;
%y = zeros(1,length(x));
%y(:) = 41.701;
b = zeros(1,length(x));
b(:) = 2*1.27;
bGlue = zeros(1,length(x));
bGlue(:) = 10;
%Qcent = zeros(1,length(x));
%Qcent(:) = 6248.5;
%Qglue = zeros(1,length(x));
%Qglue(:) = 4148.2;
t = zeros(1,length(x));
t(:) = 1.27;
height = zeros(1,length(x));
height(:) = 75;

%% 1. Point Loading Analysis (SFD, BMD) 
[SFD_PL, BMD_PL] = ApplyPL(550, P, x, SFD_PL, locationA, locationB); % Construct SFD, BMD  
[SFD_PL, BMD_PL] = ApplyPL(L, P, x, SFD_PL, locationA, locationB);   % Construct SFD, BMD


%% Calculate Failure Moments and Shear Forces 
V_Mat = Vfail(I, b, Qmax, TauU);
V_Buck = VfailBuck( I, t, b, height, E, mu, yBar, a, TauU, Qmax); 
V_Glue = Vfail(I, bGlue, Qglue, TauU);
  
M_MatT = MfailMatT(I, yBar, height, SigT, BMD_PL); 
M_MatC = MfailMatC(I, yBar, height, SigC, BMD_PL);
M_Buck1 = MfailBuck(I, t, (100-80)/2, E, mu, BMD_PL, 2, yBar, height, SigC); 
M_Buck2 = MfailBuck(I, t, 80-2*1.27, E, mu, BMD_PL, 1, yBar, height, SigC); 
M_Buck3 = MfailBuck(I, t, 75-41.701-1.27, E, mu, BMD_PL, 3, yBar, height, SigC); 
M_Buck4 = MfailBuck(I, t, (100-80)/2, E, mu, BMD_PL, 2, yBar, height, SigC);
M_Buck5 = MfailBuck(I, t, 75-41.701-1.27, E, mu, BMD_PL, 1, yBar, height, SigC);
M_Buck6 = MfailBuck(I, t, 41.701-1.27, E, mu, BMD_PL, 3, yBar, height, SigC);
M_Buck7 = MfailBuck(I, t, 80-2*1.27, E, mu, BMD_PL, 1, yBar, height, SigC);
M_Buck8 = MfailBuck(I, t, 41.701-1.27, E, mu, BMD_PL, 3, yBar, height, SigC);

%%Calculate Failure Load 
Pf = FailLoad(P, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1); 
  
%% Visualization 
VisualizePL(x, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3,M_Buck4,M_Buck5,M_Buck6, M_Buck7, M_Buck8, Pf); 
  
%% 5. Curvature, Slope, Deflections 
%Defls = Deflections(x, BMD_PL, I, E);
%disp(Defls)






%% Train failure Functions:

function [ ShearFail ] = VStress(I, b, Qcent, SFD ) 
% Calculates shear forces at every value of x that would cause a matboard shear failure 
% Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property) 
% Output: V_fail a 1-D array of length n 
    
    ShearFail = ((SFD .* Qcent) ./ (I .* b));
end

function [ ShearBuckle ] = VBuck( I, t, b, height, E, mu, y, a, Qcent )  
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web 
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property) 
% Output: V_Buck a 1-D array of length n
   
    %pi
    pi = 355/113;

    for i = 1:length(I)      
        ShearBuckle(i) = ((5*E*pi^2)/(12*(1-mu^2))) * ((t(i)/height(i))^2 + (t(i)/a(i))^2);
    end

end

function [ BendingTensionStress ] = MStressT( I,y, height, BMD )  %Tension
% Calculates bending moments at every value of x that would cause a matboard tension failure 
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
% Output: M_MatT a 1-D array of length n 
  
for i = 1 : length(I)  

        ybot = y; 
        ytop = height(i) - y; %I need somesort of height vector

        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom 
            BendingTensionStress(i) = (BMD(i) * ybot(i)) / I(i);
        elseif BMD < 0 % If the moment is negative, the tension failure will be at the top 
            BendingTensionStress(i) = (BMD(i) * ytop(i)) / I(i);
        end 
end 

end 

function [ BendingCompressionStress ] = MStressC( I,y, height, BMD )  %Compression
% Calculates bending moments at every value of x that would cause a matboard Compression failure 
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
% Output: M_MatT a 1-D array of length n 
  
for i = 1 : length(I)  

        ybot = y; 
        ytop = height(i) - y; %I need somesort of height vector

        if BMD(i) > 0 % If the moment is positive, the compression failure will be at the top 
            BendingCompressionStress(i) = (BMD(i) * ytop(i)) / I(i);
        elseif BMD(i) > 0 % If the moment is negative, the compression failure will be at the bottom 
            BendingCompressionStress(i) = (BMD(i) * ybot(i)) / I(i);
        end 
end 

    
end 

function [ MBuckStress ] = MBuck( I, t, bBuck, E, mu, caseNum )  
% Calculates bending moments at every value of x that would cause a buckling failure 
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array) 
% Output: M_MatBuck a 1-D array of length n 

    if caseNum == 1
        k = 4;
    end
    if caseNum == 2
        k = 0.425;
    end
    if caseNum == 3
        k = 6;
    end

    for i = 1:length(I)
        MBuckStress(i) = ((k*E*pi^2)/(12*(1-mu^2))) * (t(i)/bBuck)^2;
    end
    

end

function [ FOS ] = FOSCalc( ShearFail, BendingTensionStress, BendingCompressionStress, TensionCapacity, CompressionCapacity, CompressionBucklingCapacity1, CompressionBucklingCapacity2, CompressionBucklingCapacity3, CompressionBucklingCapacity4, CompressionBucklingCapacity5, CompressionBucklingCapacity6, CompressionBucklingCapacity7, CompressionBucklingCapacity8, TauCapacity, TauBuckling )
    ShearFOS = min(TauCapacity ./ abs(ShearFail));
    ShearBucklingFOS = min(TauBuckling ./ abs(ShearFail));
    TensionFOS = (TensionCapacity ./ abs(BendingTensionStress));
    CompressionFOS = min(CompressionCapacity ./ abs(BendingCompressionStress));
    CompressionBucklingFOS1 = min(CompressionBucklingCapacity1 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS2 = min(CompressionBucklingCapacity2 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS3 = min(CompressionBucklingCapacity3 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS4 = min(CompressionBucklingCapacity4 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS5 = min(CompressionBucklingCapacity5 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS6 = min(CompressionBucklingCapacity6 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS7 = min(CompressionBucklingCapacity7 ./ abs(BendingCompressionStress) );
    CompressionBucklingFOS8 = min(CompressionBucklingCapacity8 ./ abs(BendingCompressionStress) );


    FOSes = [ShearFOS ShearBucklingFOS TensionFOS CompressionFOS CompressionBucklingFOS1 CompressionBucklingFOS2 CompressionBucklingFOS3 CompressionBucklingFOS4 CompressionBucklingFOS5 CompressionBucklingFOS6 CompressionBucklingFOS7 CompressionBucklingFOS8];
    %disp(FOSes)
    FOS = min(FOSes);
    disp(FOS)
end

%Cross Section Functions
function [I, yBar] = SectionProperties( I, yBar, xStart, xEnd, bVector, hVector, yVector ) % Calculates important sectional properties. Including but not limited to ybar, I, Q, etc.
   % Input: Geometric Inputs. Format will depend on user
   % Output: Sectional Properties at every value of x. Each property is a 1-D array of length n
 
   % finding global centroid
   for i = xStart : xEnd
       sumAY = 0;
       sumA = 0;
 
       for j = 1:length(bVector)
           sumAY = sumAY + bVector(j) * hVector(j) * yVector(j);
           sumA = sumA + bVector(j) * hVector(j);
       end
 
       yBar(i) = sumAY/sumA;
 
       %finding I
       for j = 1:length(bVector)
           I(i) = I(i) + (bVector(j)*hVector(j)^3)/12 + bVector(j)*hVector(j)*(yVector(j) - yBar(i))^2;
       end
 
   end
 
end
 
 
function [Qmax, Qglue] = calculateQ( Qmax, Qglue, xStart, xEnd, cutoffAreaCent, cutoffAreaGlue, yCentCutoff, yGlueCutoff)
   
   for i = xStart : xEnd
       %finding Q at centroid
       for j = 1:length(cutoffAreaCent)
           Qmax(i) = Qmax(i) + cutoffAreaCent(j)*yCentCutoff(j);
       end
 
       %finding Q at glue
       for j = 1:length(cutoffAreaGlue)
           Qglue(i) = Qglue(i) + cutoffAreaGlue(j)*yGlueCutoff(j);
       end
 
   end
 
end

function [] = VisualizeTL(x, SFD, BMD, ShearStressDemand, TauU, ShearBuckleCapacity, TensionDemand, CompressionDemand, SigT, SigC, BendingBuck1, BendingBuck2, BendingBuck3, BendingBuck4, BendingBuck5, BendingBuck6, BendingBuck7, BendingBuck8, FOS)
   ShearTau = zeros(1,length(x));
   ShearTau(:) = TauU;
   TensionSig = zeros(1,length(x));
   TensionSig(:) = SigT;
   
   CompressionSig(:) = SigC;

 
%Plots all outputs of design process
   figure;
 
   hold on;
   %shear force diagram
   subplot(2, 3, 1)
   plot(x, SFD, "k")
   title1 = strcat("SFD from Train ", "FOS of bridge = ", int2str(FOS), " N");
   title(title1)
  
   %shear stress vs mat shear fail
   subplot(2, 3, 2)
   plot(x, ShearStressDemand, "k")
   hold on
   %plot(x, ShearGlue, "r")
   hold on
   plot(x, ShearTau, "g")
   hold on
   legend("Shear Demand", "Matboard Shear Failure", "Buckle Failure",'Location','northwest','NumColumns',2)
   title("Shear Demand vs Material Shear Failure")
   hold off
 
   subplot(2, 3, 3)
   plot(x, ShearStressDemand, "k")
   hold on
   plot(x, ShearBuckleCapacity, "g")
   hold on
   legend("Shear Demand", "Buckling Failure", 'Location','northwest','NumColumns',2)
   title("ShearDemand vs Shear Buckling Failure")
   hold off
  
   %BMD Diagrams
   subplot(2, 3, 4)
   plot(x, BMD, "k")
   hold on
   title("BMD of Train Load")
   set(gca, 'YDir', 'reverse')
   hold off
 
   subplot(2, 3, 5)
   plot(x, TensionDemand, "k")
   hold on
   plot(x, TensionSig, "r")
   hold on
   legend("Tension Demand", "Matboard Tension Failure", 'Location','northwest','NumColumns',2)
   title("Tension Demand vs Matboard Tension Capacity")
   hold off
 
   subplot(2, 3, 6)
   plot(x, CompressionDemand, "-k")
   hold on
   plot(x, CompressionSig, "-r")
   hold on
   plot(x, BendingBuck1, "-r")
   hold on
   plot(x, BendingBuck2, "-b")
   hold on
   plot(x, BendingBuck3, "-g")
   hold on
   plot(x, BendingBuck4, "-y")
   hold on
   plot(x, BendingBuck5, "-p")
   hold on
   plot(x, BendingBuck6, "-c")
   hold on
   plot(x, BendingBuck7, "-m")
   hold on
   plot(x, BendingBuck8, "-o")
   hold on
   legend("Bending Demand", "Compression Capacity", 'Location','northwest','NumColumns',2)
   title("Compression Demand vs Moment Buckling Failures")
   hold off
 
 
end


%%Point Load Functions:
function [ SFD, BMD ] = ApplyPL( xP, P, x, SFD, locationA, locationB ) %(Location, Force, x coordinate, )
% Constructs SFD and BMD from application of 1 Point Load. Assumes fixed location of supports 
% Input: location and magnitude of point load. The previous SFD can be entered as input to  
%  construct SFD of multiple point loads 
% Output: SFD, BMD both 1-D arrays of length n 
    
    %find index of x value @ 1000
    indexSuppB = length(SFD);
    for i = 1:length(x)
        if indexSuppB == length(SFD)
            if x(i) >= locationB
                indexSuppB = i;
            end
        end
    end
    
    SupportB = ((xP * P) / (locationB)); %Using Sum of Moments = 0
    SupportA = (P - SupportB); %Using sum of vertical force = 0
    
    %find index of x value = xp
    indexXP = length(SFD);
    for i = 1:length(x)
        if indexXP == length(SFD)
            if x(i) >= xP
                indexXP = i;
            end
        end
    end
    
    SFDLocal = zeros(1, length(SFD));
    if indexXP < indexSuppB
        for i = 1:(indexXP) %SFD before new load
            SFDLocal(i) = SupportA;
        end
        
        for i = indexXP+1:indexSuppB %SFD after new load
            SFDLocal(i) = (SupportA - P);
        end
        
        for i = 1:length(SFD) %add the new SFD with the old ones
            SFD(i) = SFD(i) + SFDLocal(i);
        end
    
    else
        for i = 1:(indexSuppB) %SFD before new load
            SFDLocal(i) = SupportA;
        end
        
        for i = indexSuppB+1:(indexXP) %SFD after new load
            SFDLocal(i) = (SupportA + SupportB);
        end

        for i = indexXP+1:length(SFD)
           SFDLocal(i) = 0;
        end

        for i = 1:length(SFD) %add the new SFD with the old ones
            SFD(i) = SFD(i) + SFDLocal(i);
        end
    end
    
    %form BMD:
    dx = x(2) - x(1);
    BMD = cumsum((SFD * dx)/1000);
end

function [ V_fail ] = Vfail(I, b, Qcent, TauU ) 
% Calculates shear forces at every value of x that would cause a matboard shear failure 
% Input: Sectional Properties (list of 1-D arrays), TauU (scalar material property) 
% Output: V_fail a 1-D array of length n 
 
    V_fail = TauU .* I .* b ./ Qcent;  
end
 
function [ V_Buck ] = VfailBuck( I, t, b, height, E, mu, y, a, TauU, Qcent )  
% Calculates shear forces at every value of x that would cause a shear buckling failure in the web 
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property) 
% Output: V_Buck a 1-D array of length n
   
    %pi
    pi = 355/113;

    for i = 1:length(I)      
        sigma(i) = ((5*E*pi^2)/(12*(1-mu^2))) * ((t(i)/height(i))^2 + (t(i)/a(i))^2);
    end
    
    V_Buck = Vfail(I, b, Qcent, sigma);

    if min(sigma) < TauU
        disp(min(sigma));
    end
end

function [ M_MatT ] = MfailMatT( I,y, height, SigT, BMD )  
% Calculates bending moments at every value of x that would cause a matboard tension failure 
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
% Output: M_MatT a 1-D array of length n 
  
for i = 1 : length(I)  

        ybot = y; 
        ytop = height(i) - y; %I need somesort of height vector

        if BMD(i) > 0 % If the moment is positive, the tension failure will be at the bottom 
            M_MatT(i) = (SigT * I(i)) / ybot(i); 
        elseif BMD(i) < 0 % If the moment is negative, the tension failure will be at the top 
            M_MatT(i) = (SigT * I(i)) / ytop(i); 
        end 
    end 
end 
 
function [ M_MatC ] = MfailMatC( I,y, height, SigC, BMD )  
% Calculates bending moments at every value of x that would cause a matboard tension failure 
% Input: Sectional Properties (list of 1-D arrays), SigT (material property), BMD (1-D array) 
% Output: M_MatT a 1-D array of length n 

ybot = y; 
ytop = height - y; %I need somesort of height vector
  
for i = 1 : length(I)    
        if BMD(i) > 0 % If the moment is positive, the compression failure will be at the top
            M_MatC(i) = (SigC * I(i)) / ytop(i); 
        elseif BMD(i) < 0 % If the moment is negative, the compression failure will be at the bottom
            M_MatC(i) = (SigC * I(i)) / ybot(i); 
        end 
    end 
end 
 
function [ M_Buck ] = MfailBuck( I, t, bBuck, E, mu, BMD, caseNum, y, height, SigC )  
% Calculates bending moments at every value of x that would cause a buckling failure 
% Input: Sectional Properties (list of 1-D arrays), E, mu (material property), BMD (1-D array) 
% Output: M_MatBuck a 1-D array of length n 

    if caseNum == 1
        k = 4;
    end
    if caseNum == 2
        k = 0.425;
    end
    if caseNum == 3
        k = 6;
    end

    for i = 1:length(I)
        sigma(i) = ((k*E*pi^2)/(12*(1-mu^2))) * (t(i)/bBuck)^2;
    end
    
    M_Buck = MfailMatC(I,y, height, min(sigma), BMD);

end

function [ Pf ] = FailLoad( P, SFD, BMD, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_buck )  
% Calculates the magnitude of the load P that will cause one of the failure mechanisms to occur 
% Input: SFD, BMD under the currently applied points loads (P) (each 1-D array of length n) 
%  {V_Mat, V_Glue, ... M_MatT, M_MatC, ... } (each 1-D array of length n) 
% Output: Failure Load value Pf 

%Flexural Stress
CompressionFailure = M_MatC ./ abs(BMD);
TensionFailure = M_MatT ./ abs(BMD);
BendingBucklingFailure = M_buck ./ abs(BMD);

%Shear Stress
ShearFailureMat = V_Mat ./ abs(SFD);
ShearFailureGlue = V_Buck ./ abs(SFD);
ShearBuckliingFailure = V_Buck ./ abs(SFD);

%Find the lowest Pf
Pforces = [min(CompressionFailure) min(TensionFailure) min(BendingBucklingFailure) min(ShearFailureMat) min(ShearFailureGlue) min(ShearBuckliingFailure)];
Pf = min(Pforces);

for i = 1:6
   if Pforces(i) == Pf
       disp(i)
   end
end


end
     
function [] = VisualizePL(x, SFD_PL, BMD_PL, V_Mat, V_Glue, V_Buck, M_MatT, M_MatC, M_Buck1, M_Buck2, M_Buck3, M_Buck4, M_Buck5, M_Buck6, M_Buck7, M_Buck8, Pf) 
   %Plots all outputs of design process
   figure;
 
   hold on;
   %shear force diagram
   subplot(2, 3, 1)
   plot(x, SFD_PL, "k")
   title1 = strcat("SFD from P = ", int2str(floor(Pf)), "N Pfail = ", int2str(Pf), " N");
   title(title1)
  
   %shear force vs mat shear fail
   subplot(2, 3, 2)
   plot(x, SFD_PL, "k")
   hold on
   plot(x, V_Mat, "r")
   hold on
   plot(x, V_Glue, "g")
   hold on
   legend("Shear Force", "Matboard Shear Failure", "Glue Shear Failure",'Location','northwest','NumColumns',2)
   title("SFD vs Material Shear Failure")
   hold off
 
   subplot(2, 3, 3)
   plot(x, SFD_PL, "k")
   hold on
   plot(x, V_Buck, "r")
   hold on
   legend("Shear Force", "Buckling Failure", 'Location','northwest','NumColumns',2)
   title("SFD vs Shear Buckling Failure")
   hold off
  
   %BMD Diagrams
   subplot(2, 3, 4)
   plot(x, BMD_PL, "k")
   hold on
   title2 = strcat("BMD from P = ", int2str(floor(Pf)), "N Pfail = ", int2str(Pf), " N");
   title(title2)
   set(gca, 'YDir', 'reverse')
   hold off
 
   subplot(2, 3, 5)
   plot(x, BMD_PL, "k")
   hold on
   plot(x, M_MatT, "r")
   hold on
   plot(x, M_MatC, "b")
   hold on
   legend("Bending Moment", "Matboard Tension Failure", "Matboard Compression Failure", 'Location','northwest','NumColumns',2)
   title("BFD vs Material Moment Failures")
   hold off
 
   subplot(2, 3, 6)
   plot(x, BMD_PL, "-k")
   hold on
   plot(x, M_Buck1, "-r")
   hold on
   plot(x, M_Buck2, "-b")
   hold on
   plot(x, M_Buck3, "-g")
   hold on
   plot(x, M_Buck4, "-y")
   hold on
   plot(x, M_Buck5, "-p")
   hold on
   plot(x, M_Buck6, "-c")
   hold on
   plot(x, M_Buck7, "-m")
   hold on
   plot(x, M_Buck8, "-o")
   hold on
   %legend("Bending Moment", "Mid Flange Buckling", "Side Flange Buckling", "Web Compression Buckling", 'Location','northwest','NumColumns',2)
   title("BMD vs Moment Buckling Failures")
   hold off
 
 
end

function [ Defls ] = Deflections( x, BMD, I, E )  
% Calculates deflections 
% Input: I(1-D arrays), E (material property), BMD (1-D array) 
% Output: Deflection for every value of x (1-D array) or for the midspan only   
    
% %phi = BMD ./ (I .* E);
% %dx = x(2) - x(1);
% %TotalArea = zeros(1,length(x));
% 
% %for i = 1:length(x)
% %    yPHI = phi(1:i);
% %    TotalArea(i) = trapz(yPHI);
% %end
% 
% xes = (cumsum(phi .* dx) .* x) ./ TotalArea;
% disp(xes)
% Defls = ((cumsum(phi .* dx) .* x) ./ TotalArea) .* cumsum(phi .* dx);

P = 200;

SFD_PL = zeros(1, length(x));      % Initialize SFD(x) 

locationA = 1;
locationB = 1060;

[SFD_200, BMD_200] = ApplyPL(550, P, x, SFD_PL, locationA, locationB); % Construct SFD, BMD  
[SFD_200, BMD_200] = ApplyPL(1250, P, x, SFD_PL, locationA, locationB);   % Construct SFD, BMD

phi = (BMD_200*1000) ./ (I * E);
dx = x(2) - x(1);
xBar1 = x - x(530);
delta_BA = sum(phi .* dx .* xBar1);

xBar2 = x - x(1000);
delta_DA = sum(phi .* dx .* xBar2) ;

Defls = ((x(1000)/x(500)) * abs(delta_DA)) - abs(delta_BA);


    
end
