function [new_coord, mirror_x, mirror_y, swap_xy] = ConvertCoord(x, y, z, num_x, num_y)
x_mid = ceil(num_x/2);
y_mid = ceil(num_y/2);

new_z = z;
mirror_x = 0;
mirror_y = 0;
swap_xy = 0;

if x < x_mid
    x = num_x - x +1;
    mirror_x = 1;
end
if y < y_mid
    y = num_y - y + 1;
    mirror_y = 1;
end

new_x = min(x, y);
new_y = max(x, y);
if new_x~=x || new_y~=y
    swap_xy = 1;
end

new_coord = [new_x, new_y, new_z];
