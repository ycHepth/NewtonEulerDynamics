clc;
clear;

syms q1 q2 m1 m2 L1 L2 real

% modified DH = [alpha, a, d, theta]
dh_params = [0, 0, 0, q1;
             0, L1,0, q2;
             0, L2,0, 0];
         
mass_center= [L1,0,0;
              L2,0,0;];
          
mass = [m1 m2];

inertia_1 = [0 0 0;
             0 0 0;
             0 0 0];
         
inertia_2 = [0 0 0;
             0 0 0;
             0 0 0];
         
inertia_tensor_list(:,:,1) = inertia_1;
inertia_tensor_list(:,:,2) = inertia_2;

f_ext = [0 0 0;
         0 0 0];
     
torque = NewtonEulerDynamics(dh_params, mass, mass_center, inertia_tensor_list, f_ext);
