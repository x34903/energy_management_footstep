oo_constants

%% Basic Functions used throughout

syms xi yi roti 
Trans_rot = [cos(roti) -sin(roti) 0 xi ;
             sin(roti)  cos(roti) 0 yi ;
                     0          0 1  0 ;
                     0          0 0  1 ] ;

matlabFunction( Trans_rot , 'File', 'symbolic_functions/Translate_rotate', 'Vars', {xi,yi,roti}) ; 


%% no conjugate cross product
avect = sym('avect',[3,1]) ;
bvect = sym('bvect',[3,1]) ;
my_cross = [ avect(2)*bvect(3) - avect(3)*bvect(2) ;
             avect(3)*bvect(1) - avect(1)*bvect(3) ;
             avect(1)*bvect(2) - avect(2)*bvect(1) ] ;
      
matlabFunction( my_cross , 'File', 'symbolic_functions/mycross_noconj', 'Vars', {avect,bvect}) ; 

