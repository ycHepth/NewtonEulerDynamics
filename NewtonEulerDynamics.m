function torque = NewtonEulerDynamics(dh_list,mass_list,mass_center_list,inertia_tensor_list,f_external)
% Newton_Euler approach to obtain robot linkage dynamics
% Input:
%       dh_list: modified DH list
%       f_external: force/torque applied to end-linkage
%       inertia_tensor_list: inertia tensor of mass_center coordinate( aligned with linkage coordinate)
% Output:
%       torque_list: [q,dq,ddq]
% 
% Note:
%       1. the DH coordinate must be modified coordinate.
%       2. required inertia tensor before modeling
%       3. for no disturbance case, f_external = 0

[rows,cols] = size(dh_list);
number_of_links = rows-1;
if cols ~= 4
    error('wrong DH paralist')
end

T = sym([]);
R = sym([]);
a = sym([]);
d = sym([]);
alpha = sym([]);
theta = sym([]);

for i = 1:rows
    eval(['syms ','q',num2str(i),' real;']);
    eval(['syms ','dq',num2str(i),' real;']);
    eval(['syms ','ddq',num2str(i),' real;']);
    eval(['q(i)=','q',num2str(i),';']);
    eval(['dq(i)=','dq',num2str(i),';']);
    eval(['ddq(i)=','ddq',num2str(i),';']);
end

for i = 1:rows
    dh = dh_list(i,:);
    alpha(i) = dh(1);
    a(i) = dh(2);
    d(i) = dh(3);
    theta(i) = dh(4);
    
    if i == rows
        q(i) = 0;
    end
    
%   actuated axis transform matrix T(i-i -> i):

   T(:,:,i) =  [cos(q(i))              ,-sin(q(i))              ,  0            ,  a(i);
                sin(q(i))*cos(alpha(i)), cos(q(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i);
                sin(q(i))*sin(alpha(i)), cos(q(i))*sin(alpha(i)),  cos(alpha(i)),  cos(alpha(i))*d(i);
                0                      , 0                      ,  0            ,  1                 ;];
            
   T = T(:,:,i);
   
%    extract inverse matrix of rotate matrix R
   R(:,:,i) = simplify(inv(T(1:3,1:3)));
   P(:,:,i) = T(1:3,4:4);
   
end

   z = [0;0;1];
   syms g real;
   
%    forward
   for i = 0:number_of_links-1
       if i == 0
           wi  = [0 0 0]';
           dwi = [0 0 0]';
           dvi = [0 g 0]';
       else
           wi  = w(:,i);
           dwi = dw(:,i);
           dvi= dv(:,i);
       end
       
       w(:,:,i+1)   = R(:,:,i+1)*wi + dq(i+1)*z;
       dw(:,:,i+1)  = R(:,:,i+1)*dwi + cross(R(:,:,i+1)*wi,dq(i+1)*z) + ddq(i+1)*z;
       dv(:,:,i+1)  = R(:,:,i+1)*(cross(dwi,P(:,:,i+1)) + cross(wi,cross(wi,P(:,:,i+1))) + dvi);
       dvc(:,:,i+1) = cross(dw(:,:,i+1),mass_center_list(i+1,:)')...
                    + cross(w(:,:,i+1),cross(w(:,:,i+1),mass_center_list(i+1,:)'))...
                    + dv(:,:,i+1);
       F(:,:,i+1) = mass_list(i+1)*dvc(:,:,i+1);
       N(:,:,i+1) = inertia_tensor_list(:,:,i+1)*dw(:,:,i+1) + cross(w(:,:,i+1),inertia_tensor_list(:,:,i+1)*w(:,:,i+1));    
   end
   
   f = sym([]);
   n = sym([]);
   
%    backward
    for i = number_of_links:-1:1
        if i == number_of_links  
            f(:,:,i+1) = f_external(1,:)';
            n(:,:,i+1) = f_external(2,:)';
        end
    
    f(:,:,i) = R(:,:,i+1)\f(:,:,i+1) + F(:,:,i);
    f(:,:,i) = simplify(f(:,:,i));
    n(:,:,i) = N(:,:,i) + R(:,:,i+1)\n(:,:,i+1) + cross(mass_center_list(i,:)',F(:,:,i))...
                + cross(P(:,:,i+1),R(:,:,i+1)\f(:,:,i+1));
    n(:,:,i) = simplify(n(:,:,i));
    torque_list(i) = dot(n(:,:,i),z);
    end
torque = torque_list';


    
   
   
       
    
   
   
   
   
   
   
   
            