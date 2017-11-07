function f = gaussian_bubble_source

f =@(x,y,path_x,path_y)  30000.0  * exp(- (x-path_x).^2/.1.^2 - (y-path_y).^2/.1.^2 );

end