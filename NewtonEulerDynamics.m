function torque = NewtonEulerDynamics(dh_list,mass_list,mass_center_list,inertia_tensor_list,f_external)
% Newton_Euler approach to obtain robot linkage dynamics
% Input:
%       dh_list: modified DH list
%       f_external: force/torque applied to end-linkage
%       inertia_tensor_list: inertia tensor at origin coordinate.
% Output:
%       torque_list: [q,dq,ddq]
% 
% Note:
%       1. DH parameters based on modified DH approach.
%       2. required inertia tensor of linkage system before modeling
%       3. for no disturbance case, f_external = 0, which consist of force
%       and torque.

[rows,cols] = size(dh_list);
number_of_links = rows-1;
if cols ~= 4
    error('wrong DH paralist')
end

T = sym([]);
R = sym([]);
P = sym([]);
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
   
%    forward : 0 --> links-1
   for i = 0:number_of_links-1
       if i == 0
           wi  = [0 0 0]';
           dwi = [0 0 0]';
           dvi = [0 g 0]';
%      here applied 'g' at y-axis meaning robot side-assembly
%      the specify direction of 'g' associated with robot coordinate
%      establish.
       else
           wi  = w(:,i);
           dwi = dw(:,i);
           dvi= dv(:,i);
       end
%        angle
       w(:,:,i+1)   = R(:,:,i+1)*wi + dq(i+1)*z;    
%        angle velocity
       dw(:,:,i+1)  = R(:,:,i+1)*dwi + cross(R(:,:,i+1)*wi,dq(i+1)*z) + ddq(i+1)*z;
%        origin acceleration
       dv(:,:,i+1)  = R(:,:,i+1)*(cross(dwi,P(:,:,i+1)) + cross(wi,cross(wi,P(:,:,i+1))) + dvi);
%        mass center acceleraion

% ==========   !!! NOTE   ======================
% here the 'dvc' is acceleration of mess center,
%       --- ref : https://zhuanlan.zhihu.com/p/109770180
% and from the example of this project source, the mass_center is 
% actually at end of each linkage, thus equal to origin of each linkage
% =======================================

       dvc(:,:,i+1) = cross(dw(:,:,i+1),mass_center_list(i+1,:)')...
                    + cross(w(:,:,i+1),cross(w(:,:,i+1),mass_center_list(i+1,:)'))...
                    + dv(:,:,i+1);          
%  ===== above 4 parameters independt with inertia tensor =====

%   ===============            
%            sum force   F
%   ===============
%   As the next equation of F, the sum force F is associated with mass
%   center, and show a linear relationship with <mass*mass_center>
% -------------------------------------------------------------------------------
%  F = m*a + (dw)X(m*mass_center) + (w)X [(w) X (m*mass_center) ]
% -------------------------------------------------------------------------------
%  using cross product matrix S(w) form:
% 
%  F = m*a + S(dw)*m*mass_center + S(w)S(w)*(m*mass_center)
% -------------------------------------------------------------------------------
%  expanded inertia tensor I to a 10x1 vector P
%  P = [Ixx Ixy Ixz Iyy Iyz Izz m*xc m*yc m*zc m]'  
%  F = [ 0 | S(dw) + S(w)S(w), a] * P
%  define H = [ 0 | S(dw) + S(w)S(w), a ]
%  then F = HP
%  thus matrix H is independ with inertia parameters which has 3x10
%  dimension( 3 for three direction of force and 10 for parameters )
% -------------------------------------------------------------------------------

       F(:,:,i+1) = mass_list(i+1)*dvc(:,:,i+1);
       
%  =============== 
%          sum torque N
%  ===============
%  same as sum force F
%  N is associated with inertia parameters.
%  notice :
% I*w = [wx wy wz 0 0 0; 0 wx 0 wy wz 0; 0 0 wx 0 wy wz] 
%       *  [Ixx Ixy Ixz Iyy Iyz Izz]'
%   let K(w) = [wx  wy  wz   0     0   0  ]
%                    [0    wx   0    wy  wz  0  ]
%                    [0    0     wx   0   wy  wz]
%  thus N can be simplify as
%  N = [K(dw) + S(w)K(w) | 0] * P
%  let A = [K(dw) + S(w)K(w) | 0] 
%  then N = AP

       N(:,:,i+1) = inertia_tensor_list(:,:,i+1)*dw(:,:,i+1) + cross(w(:,:,i+1),inertia_tensor_list(:,:,i+1)*w(:,:,i+1));    
       
   end
%  ===== above 2 parameters used for seperating inertia parameters =====

   f = sym([]);
   n = sym([]);
   
%    backward  1 <-- links
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

%    ===== above all parameters used for seperating inertia parameters ===== 

    torque_list(i) = dot(n(:,:,i),z);
    end
    
torque = torque_list';

%     define n = YP, where P is inertia parameters vector
%     then for number of linkage = 6 
%     f_6 = [ 0 0 0 0 0 H6]P = Yf_6*P
%     n_6 = [ 0 0 0 0 0 A6]P = Y_6*P
%     f_5 = H5*P5 + R6\f_6 = H5*P5 + R6\Yf_6*P
%     n_5 = [[ 0 0 0 0 H5 0] + R6\Yf_6] P
%     = [[ 0 0 0 0 A5 0] + R6\Y6 + S(Pos_6)R6\Yf_6 + S(center_6)F6] P
%     all in all, we can obtain following equations
%     n_4  =[[0 0 0 A4 0 0] + R5\Y5 + S(Pos_5)R5\Yf_5 + S(center_5)F5] P
%     n_3  =[[0 0 A3 0 0 0] + R4\Y4 + S(Pos_4)R4\Yf_4 + S(center_4)F4] P
%     n_2  =[[0 A2 0 0 0 0] + R3\Y3 + S(Pos_3)R5\Yf_3 + S(center_3)F3] P
%     n_1  =[[A1 0 0 0 0 0] + R2\Y2 + S(Pos_2)R5\Yf_2 + S(center_2)F2] P
%  for recognize matrix Y
%     Y_6 = [ 0 0 0 0 0 A6]
%     Y_5 = [[0 0 0 0 A5 0] + R6\Y6 + S(Pos_6)R6\Yf_6 + S(center_6)F6]
%     Y_4 = [[0 0 0 A4 0 0] + R5\Y5 + S(Pos_5)R5\Yf_5 + S(center_5)F5] 
%     Y_3 = [[0 0 A3 0 0 0] + R4\Y4 + S(Pos_4)R4\Yf_4 + S(center_4)F4]
%     Y_2 = [[0 A2 0 0 0 0] + R3\Y3 + S(Pos_3)R5\Yf_3 + S(center_3)F3]
%     Y_1 = [[A1 0 0 0 0 0] + R2\Y2 + S(Pos_2)R5\Yf_2 + S(center_2)F2]

%    torque = [Z'Y1 Z'Y2 Z'Y3 Z'Y4 Z¡®Y5 Z'Y6]' P 
%    the result of linearization of dynamics:
%    Y = [Z'Y1 Z'Y2 Z'Y3 Z'Y4 Z¡®Y5 Z'Y6]'
   
       
    
   
   
   
   
   
   
   
            