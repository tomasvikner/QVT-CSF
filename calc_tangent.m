function [Tangent_V, V2, V3] = calc_tangent(branchList)

d = 2; % dist. behind/ahead of current pt for tangent plane calc (d=2->5pts)
Tangent_V = zeros(0, 3);
for n = 1:max(branchList(:, 4))
    branchActual = branchList(branchList(:, 4) == n, :);
    nPts = size(branchActual, 1);
    if nPts < 1
        continue
    end
    dir_temp = zeros(nPts, 3);
    for i = 1:nPts
        % Indices clamped so short branches / single-point stubs never exceed bounds
        i_lo = max(1, i - d);
        i_hi = min(nPts, i + d);
        if i_hi > i_lo
            dir = branchActual(i_hi, 1:3) - branchActual(i_lo, 1:3);
        elseif nPts >= 2
            % Two points: segment direction
            if i < nPts
                dir = branchActual(i + 1, 1:3) - branchActual(i, 1:3);
            else
                dir = branchActual(i, 1:3) - branchActual(i - 1, 1:3);
            end
        else
            % Single-point stub: arbitrary tangent so downstream geometry does not NaN
            dir = [0, 0, 1];
        end
        nd = norm(dir(:));
        if nd > eps
            dir_temp(i, :) = dir(:).' / nd;
        else
            dir_temp(i, :) = [0, 0, 1];
        end
    end
    Tangent_V = [Tangent_V; dir_temp]; %#ok<*AGROW>
end

if isempty(Tangent_V)
    Tangent_V = zeros(size(branchList, 1), 3);
    Tangent_V(:, 3) = 1;
end

% This will find a normalized vector perpendicular to the tangent vector
[~, idx_max] = max(abs(Tangent_V), [], 2); %get max unit along rows
idx_max(idx_max == 2) = 1; %flatten to 2D
max_pts = sub2ind(size(Tangent_V), (1:size(Tangent_V, 1))', idx_max);
temp = zeros(size(Tangent_V));
temp(max_pts) = 1; %binary matrix of location of max unit vectors
[~, idx_shift] = max(abs(circshift(temp, 1, 2)), [], 2); %rotate (ie x->y,z->x)
shift_pts = sub2ind(size(Tangent_V), (1:size(Tangent_V, 1))', idx_shift);
V2 = zeros(size(Tangent_V));
V2(max_pts) = Tangent_V(shift_pts);
V2(shift_pts) = -Tangent_V(max_pts);
N = repmat(sqrt(sum(abs(V2) .^ 2, 2)), [1 3]); %repeat vel. magnitude as Nx3
N(N < eps) = 1;
V2 = V2 ./ N;
V3 = cross(Tangent_V, V2); %Third vector that is normalized
% V3,V2,Tangent_V are all orthogonal (i.e. dot( V3(1,:),Tangent_V(1,:) )=0)

end
