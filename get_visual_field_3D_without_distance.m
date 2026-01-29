function [r,XvRot] = get_visual_field_3D_without_distance(v,xyz,index)
%v_copy=[theta,phi]
v_copy=zeros(numel(v(:,1)),2);
v_copy(:,1)=atan2(v(:, 3),sqrt(v(:,1).^2 +v(:,2).^2));
v_copy(:,2)=atan2(v(:, 2), v(:, 1));
Xv=xyz-xyz(index,:);

%XvRot:rotate and translate each position around the focal individual index
XvRot=zeros(size(Xv));
XvRot(:,1)=Xv(:,1)*cos(-v_copy(index,2))-Xv(:,2)*sin(-v_copy(index,2));
XvRot(:,2)=Xv(:,1)*sin(-v_copy(index,2))+Xv(:,2)*cos(-v_copy(index,2));
XvRot(:,3)=Xv(:,3);

Xv=XvRot;
XvRot(:,1)=Xv(:,1)*cos(-v_copy(index,1))-Xv(:,3)*sin(-v_copy(index,1));
XvRot(:,2)=Xv(:,2);
XvRot(:,3)=Xv(:,1)*sin(-v_copy(index,1))+Xv(:,3)*cos(-v_copy(index,1));

% r=[r,theta,phi]
r=cart2sph(XvRot);

% Xv=xyz-xyz(index,:);
% for i = 1 : 7
%     rotated_vector(i,:) = rotateVectorToReference(Xv(i,:), v(index,:));
%     rotated_vector1(i,:) = rotateVectorToReference1(Xv(i,:), v(index,:));
%     rotated_vector2(i,:) = rotateVectorToReference2(Xv(i,:), v(index,:));
% 
% end
end

function r=cart2sph(x)
r=zeros(numel(x(:,1)),3);
r(:,1)=sqrt(x(:,1).^2 +x(:,2).^2+x(:,3).^2);
r(:,2)=atan2(x(:, 3),sqrt(x(:,1).^2 +x(:,2).^2));
r(:,3)=atan2(x(:, 2), x(:, 1));
end

function rotated_vector = rotateVectorToReference(vector, reference)
    % Normalize the input vectors
    vector_norm = vector / norm(vector);
    reference_norm = reference / norm(reference);

    % Calculate the axis of rotation using cross product
    axis = cross(vector_norm, reference_norm);

    % Calculate the angle of rotation using dot product
    angle = acos(dot(vector_norm, reference_norm));

    % Create the rotation matrix
    c = cos(angle);
    s = sin(angle);
    rot_matrix = [c + axis(1)^2*(1-c) axis(1)*axis(2)*(1-c) - axis(3)*s axis(1)*axis(3)*(1-c) + axis(2)*s;
                  axis(2)*axis(1)*(1-c) + axis(3)*s c + axis(2)^2*(1-c) axis(2)*axis(3)*(1-c) - axis(1)*s;
                  axis(3)*axis(1)*(1-c) - axis(2)*s axis(3)*axis(2)*(1-c) + axis(1)*s c + axis(3)^2*(1-c)];

    % Rotate the vector to the reference direction
    rotated_vector = rot_matrix * vector';

    % Transpose the result back to a column vector
    rotated_vector = rotated_vector';

end

function rotated_vector = rotateVectorToReference1(vector, reference)
    % Normalize the input vectors
    vector_norm = vector / norm(vector);
    reference_norm = reference / norm(reference);

    % Calculate the angle between the input vector and the reference vector
    angle = atan2(norm(cross(vector_norm, reference_norm)), dot(vector_norm, reference_norm));

    % Create the rotation matrix around the axis perpendicular to both vectors
    axis = cross(vector_norm, reference_norm);
    R = vrrotvec2mat([axis, angle]);

    % Apply the rotation matrix to the original vector
    rotated_vector = (R * vector')';

end

function rotated_vector = rotateVectorToReference2(vector, reference)
    % Normalize the input vectors
    vector_norm = vector / norm(vector);
    reference_norm = reference / norm(reference);

    % Calculate the rotation matrix to align the x-axis with the reference vector
    cos_theta = dot([1, 0, 0], reference_norm);
    axis = cross([1, 0, 0], reference_norm);
    sin_theta = norm(axis);
    kross = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];
    R = eye(3) + kross + kross^2 * (1 - cos_theta) / sin_theta^2;

    % Apply the rotation matrix to the original vector
    rotated_vector = (R * vector')';

end