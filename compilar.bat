c++ -Wall -std=c++11 fd3d_diffusion_prb.cpp fd3d_diffusion.cpp -o fd3d_diffusion


pgc++ -std:c++11 -acc -ta=multicore,tesla -Minfo=accel  fd3d_diffusion_prb.cpp fd3d_diffusion.cpp -o fd3d_diffusion 


export ACC_DEVICE_TYPE=nvidia

ssh felixmejia@toctoc.grid.uis.edu.co


ssh guane.uis.edu.co

scp * felixmejia@toctoc.grid.uis.edu.co:/home/felixmejia/FD3D





