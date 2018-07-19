function f = gaussian_bubble_source_half_radius

f =@(x,y,path_x,path_y)  100*3*190/(pi*.05^2)  * exp(- 3*(x-path_x).^2/.05.^2 - 3*(y-path_y).^2/.05.^2 );

end