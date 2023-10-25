function [point_in0] = Rototranslate(R1to0, point_in1)

point_in0 = R1to0 * [point_in1;1];
point_in0 = point_in0(1:3);

end