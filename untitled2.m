syms alpha alphast Clstall Cdstall
Cl=sin(2*alpha)+(Clstall-sin(2*alphast))*((sin(alphast))/(cos(alphast)^2))*((cos(alpha)^2)/(sin(alpha)));
Cd=2*sin(alpha)^2+(Cdstall-2*sin(alphast)^2)*(cos(alpha)/cos(alphast));