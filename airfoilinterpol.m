function [p,p2,Cl,Cd]=airfoilinterpol(Re,PN,e)
 Rtext=num2str(round(Re,0));
% Rtext=num2str(10000);
if e<=8
    Airfoil  = 'E471.dat';
else
    Airfoil  = 'CH10.dat'; 
end                                                        
AoA1   = '0'; 
AoA2   ='10';
AoA2stall =  '-10';
AoAstep ='0.5';
AoAstep2 ='0.5';
saveFlnmAoA = 'CH10alpha.txt'; 
saveFlnmAoA2 = 'CH10alpha2.txt'; 
if (exist(saveFlnmAoA,'file'))
    delete(saveFlnmAoA);
end
if (exist(saveFlnmAoA2,'file'))
    delete(saveFlnmAoA2);
end

fid = fopen('xfoil_input1.txt','w');
fprintf(fid,['Load ' Airfoil '\n']);
fprintf(fid,'GDES\n');
fprintf(fid,'CADD\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'CADD\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'CADD\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'exec\n');
fprintf(fid,'\n');
fprintf(fid,'PPAR\n');
fprintf(fid,'n\n');
fprintf(fid,'300\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'oper\n');
fprintf(fid,'v\n');
fprintf(fid,[Rtext '\n']);
fprintf(fid,'iter\n');
fprintf(fid,'110\n');
fprintf(fid, 'vpar\n');
fprintf(fid, 'n\n');
fprintf(fid, ['1' '\n']);
fprintf(fid,'\n');
fprintf(fid,'a 0 \n');
fprintf(fid,'pacc\n');
fprintf(fid,[saveFlnmAoA '\n']);
fprintf(fid,'\n');
if (~PN)
fprintf(fid,'aseq\n');
fprintf(fid,[AoA1 '\n']);
fprintf(fid,[AoA2 '\n']);
fprintf(fid,[AoAstep '\n']);
fprintf(fid,'\n');
% fprintf(fid,'pacc\n');
% fprintf(fid,'pacc\n');
% fprintf(fid,[saveFlnmAoA2 '\n']);
else
fprintf(fid,'aseq\n');
fprintf(fid,[AoA1 '\n']);
fprintf(fid,[AoA2stall '\n']);
fprintf(fid,[AoAstep2 '\n']);
fprintf(fid,'\n');
end
fprintf(fid,'quit\n');
% fprintf(fid,['aseq' AoA1 AoA2 AoAstep '\n']);
% fprintf(fid,['aseq' AoA1 AoA2 AoAstep '\n']);
fclose(fid);
cmd = 'xfoil.exe < xfoil_input1.txt';
[status,result] = system(cmd);



saveFlnmAoA = 'CH10alpha.txt';
fidsaveFlnmAoA = fopen(saveFlnmAoA);
dataBuffer = textscan(fidsaveFlnmAoA,'%f %f %f %f %f %f %f ','HeaderLines',12,...                  % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
fclose(fidsaveFlnmAoA);
delete(saveFlnmAoA);
alpha = dataBuffer{1,1}(:,1);                                               
Cl= dataBuffer{1,1}(:,2);                                              
Cd = dataBuffer{1,1}(:,3); 
% saveFlnmAoA2 = 'CH10alpha2.txt';
% fidsaveFlnmAoA2 = fopen(saveFlnmAoA2);
% dataBuffer2 = textscan(fidsaveFlnmAoA2,'%f %f %f %f %f %f %f ','HeaderLines',12,...                  % Ready data from file
%                             'CollectOutput',1,...
%                             'Delimiter','');
% fclose(fidsaveFlnmAoA2);
% delete(saveFlnmAoA2);
% alpha2 = dataBuffer2{1,1}(:,1);                                               
% Cl2= dataBuffer2{1,1}(:,2);                                              
% Cd2 = dataBuffer2{1,1}(:,3);
% p3 = polyfit(alpha2,Cl2,5);
% p4 = polyfit(alpha2,Cd2,5);
p = polyfit(alpha,Cl,5);
p2 = polyfit(alpha,Cd,5);
% x1 = linspace(0,15);
% y1 = polyval(p,x1);
% x2 = linspace(0,15);
% y2 = polyval(p2,x2);
% fclose(fidsaveFlnmAoA);
% fclose(fidsaveFlnmAoA2);
% delete(saveFlnmAoA);
% delete(saveFlnmAoA2);
% figure
% plot(alpha,Cl,'o')
% hold on
% plot(x1,y1)
% hold off
% figure
% plot(alpha,Cd,'*')
% hold on
% plot(x2,y2)
end

