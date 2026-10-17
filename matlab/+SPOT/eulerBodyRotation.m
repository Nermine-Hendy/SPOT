function R=eulerBodyRotation(angles)
%EULERBODYROTATION Active body-to-GCRF Rz(yaw)*Ry(pitch)*Rx(roll), degrees.
a=angles(1); b=angles(2); c=angles(3);
Rx=[1 0 0;0 cosd(a) -sind(a);0 sind(a) cosd(a)];
Ry=[cosd(b) 0 sind(b);0 1 0;-sind(b) 0 cosd(b)];
Rz=[cosd(c) -sind(c) 0;sind(c) cosd(c) 0;0 0 1];
R=Rz*Ry*Rx;
end
