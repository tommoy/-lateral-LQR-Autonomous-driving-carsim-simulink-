target_path_end=700;
target_path_space=0.1;
target_path_size=target_path_end/target_path_space;
target_path_x=0:target_path_space:target_path_end;
target_path_y=zeros(1,target_path_size+1);
target_path_heading=zeros(1,target_path_size+1);
target_path_curvature=zeros(1,target_path_size+1);
for n=1:target_path_size+1
    target_path_y(n)=road_y_function(target_path_x(n));
    target_path_heading(n)=road_heading_function(target_path_x(n));
    target_path_curvature(n)=road_curvature_function(target_path_x(n));
end
plot(target_path_x,target_path_y);
plot(target_path_x,target_path_heading);
plot(target_path_x,target_path_curvature);