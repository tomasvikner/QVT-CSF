function [Tangent_V, V2, V3] = calc_tangent(branchList)

d = 2; % dist. behind/ahead of current pt for tangent plane calc (d=2->5pts)
Tangent_V = zeros(0,3);
for n = 1:max(branchList(:,4))
    branchActual = branchList(branchList(:,4)==n,:);
    dir_temp = zeros(size(branchActual,1),3);
    for i = 1:size(branchActual,1)
        % Extract normal to cross-section
        if i < d+1 %if near 1st endpoint
            dir = (branchActual(i+d,1:3) - branchActual(i,1:3));
        elseif i >= size(branchActual,1)-d %if near 2nd endpoint 
            dir = (branchActual(i,1:3) - branchActual(i-d,1:3));
        else % calculate tangent from d points ahead/behind curr point
            dir = (branchActual(i+d,1:3) - branchActual(i-d,1:3));
        end
        dir_temp(i,:) = dir/norm(dir); %tangent vector with magnitude of 1
    end
    Tangent_V = [Tangent_V;dir_temp]; %#ok<*AGROW> %add all tangents to large list
end

% This will find a normalized vector perpendicular to the tangent vector
[~,idx_max] = max(abs(Tangent_V),[],2); %get max unit along rows
idx_max(idx_max==2) = 1; %flatten to 2D
max_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_max);
temp = zeros(size(Tangent_V));
temp(max_pts) = 1; %binary matrix of location of max unit vectors
[~,idx_shift] = max(abs(circshift(temp,1,2)),[],2); %rotate (ie x->y,z->x)
shift_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_shift);
V2 = zeros(size(Tangent_V));
V2(max_pts) = Tangent_V(shift_pts);
V2(shift_pts) = -Tangent_V(max_pts);
N = repmat(sqrt(sum(abs(V2).^2,2)),[1 3]); %repeat vel. magnitude as Nx3
V2 = V2./N;
V3 = cross(Tangent_V,V2); %Third vector that is normalized
% V3,V2,Tangent_V are all orthogonal (i.e. dot( V3(1,:),Tangent_V(1,:) )=0)

end
