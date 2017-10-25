function f = gaussian_bubble_source

f =@(x,y,path_x,path_y)  300.0  * exp(- (x-path_x).^2/1.0.^2 - (y-path_y).^2/1.0.^2 );

end