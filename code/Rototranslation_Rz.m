function [rotation_matrix] = Rototranslation_Rz(x,y,z, rotZ)
import casadi.*

rotation_matrix = SX.eye(4);
rotation_matrix(1:2, 1:2) = [ cos(rotZ), -sin(rotZ) ; ...
    sin(rotZ), cos(rotZ)] ;
rotation_matrix(1:3,4) = [x, y, z]';