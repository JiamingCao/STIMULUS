function [x,y,z] = polar_to_cartesian(r,theta1,phi1)
    
x = r.*sin(theta1).*cos(phi1);
y = r.*sin(theta1).*sin(phi1);
z = r.*cos(theta1);
