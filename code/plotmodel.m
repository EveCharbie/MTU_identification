function plotmodel(muscle_origins, muscle_insertions, joint_centers)
% PLOTMODEL Plots a 3D representation of muscle origins, insertions, and joint centers.
%
% USAGE:
%   plotmodel(muscle_origins, muscle_insertions, joint_centers)
%
% INPUTS:
%   muscle_origins       - A 3xN matrix where each column represents the 3D coordinates of a muscle's origin.
%   muscle_insertions    - A 3xN matrix where each column represents the 3D coordinates of a muscle's insertion.
%   joint_centers        - A 3xN matrix where each column represents the 3D coordinates of a joint center.
%
% DESCRIPTION:
%   This function creates a 3D plot where muscle origins are represented as red dots,
%   muscle insertions as blue dots, and joint centers as black dots. It also draws lines
%   between joint centers and between muscle origins and their respective insertions.
%
% NOTE:
%   It is assumed that the matrices have at least 4 columns. If not, the function might error out.
%
% EXAMPLE:
%   plotmodel(rand(3,4), rand(3,4), rand(3,4));

Onumeric = full(muscle_origins);
Inumeric = full(muscle_insertions); 
Mnumeric = full(joint_centers);
%%
figure("Name","model", "Color",[1 1 1])
plot3(Onumeric(1,:),Onumeric(2,:),Onumeric(3,:),'or')
hold on 
plot3(Inumeric(1,:),Inumeric(2,:),Inumeric(3,:),'ob')
plot3(Mnumeric(1,:),Mnumeric(2,:),Mnumeric(3,:),'ok')

for i = 1:3
    line([Mnumeric(1,i),Mnumeric(1,i+1)], [Mnumeric(2,i),Mnumeric(2,i+1)], [Mnumeric(3,i),Mnumeric(3,i+1)],'Color','black')
    line([Onumeric(1,i),Inumeric(1,i)], [Onumeric(2,i),Inumeric(2,i)], [Onumeric(3,i),Inumeric(3,i)],'Color','red')
end
legend([{'muscle origin'},{'muscle insersion'},{'joint center'}])
xlim([-1 1])
ylim([-1 1])
end
