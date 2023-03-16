clear all
CH10 = 'example.txt';
fidCH10 = fopen(CH10);
dataBuffer = textscan(fidCH10,'%f %f %f %f %f %f %f ','HeaderLines',12,...                  % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
alpha = dataBuffer{1,1}(:,1);                                               
Cl= dataBuffer{1,1}(:,2);                                              
Cd = dataBuffer{1,1}(:,3); 
p = polyfit(alpha,Cl,3);
p2 = polyfit(alpha,Cd,5);
x1 = linspace(0,15);
y1 = polyval(p,x1);
x2 = linspace(0,15);
y2 = polyval(p2,x2);
figure
plot(alpha,Cl,'o')
hold on
plot(x1,y1)
hold off
figure
plot(alpha,Cd,'*')
hold on
plot(x2,y2)
