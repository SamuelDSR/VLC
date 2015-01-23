%% Room configuration
room.length = 5;
room.width = 5;
room.height = 2.55;

%% LED Lamps configuration
led.half_angle = 65;  %degree
led.lambert_index = -1/log2(cos(led.half_angle*pi/180));
led.minFC = 100;  %minimal forward current, mA
led.maxFC = 1000; %maximal forward current, mA
led.current_flux_polynomial_coeff = [1 1 1 1];

led.fc = 700;  %forward current, mA
led.flux = 5000;   %lumen, at some FC current
led.power = led.flux/683; %watt, approximately

%% LEDs scenario configuration
%% configuration 1
led.position_x = [1.25, 3.75];
led.position_y = [1.25, 3.75]';
led.position_z = 3;
led.position_matrix_x = repmat(led.position_x,length(led.position_y),1);
led.position_matrix_y = repmat(led.position_y,1,length(led.position_x));

%% configuration 2
% position_matrix_x1 = zeros(15,15);
% position_matrix_y1 = zeros(15,15);
% for i = 1:15
%     position_matrix_x1(i,:) = linspace(0.75,1.75,15);
% end
% for j = 1:15
%     position_matrix_y1(:,j) = linspace(0.75,1.75,15);
% end
% 
% position_matrix_x2= zeros(15,15);
% position_matrix_y2 = zeros(15,15);
% for i = 1:15
%     position_matrix_x2(i,:) = linspace(3.25,4.25,15);
% end
% for j = 1:15
%     position_matrix_y2(:,j) = linspace(0.75,1.75,15);
% end
% 
% position_matrix_x3= zeros(15,15);
% position_matrix_y3 = zeros(15,15);
% for i = 1:15
%     position_matrix_x3(i,:) = linspace(0.75,1.75,15);
% end
% for j = 1:15
%     position_matrix_y3(:,j) = linspace(3.25,4.25,15);
% end
% 
% position_matrix_x4= zeros(15,15);
% position_matrix_y4 = zeros(15,15);
% for i = 1:15
%     position_matrix_x4(i,:) = linspace(3.25,4.25,15);
% end
% for j = 1:15
%     position_matrix_y4(:,j) = linspace(3.25,4.25,15);
% end
% 
% position_matrix_x= [position_matrix_x1,position_matrix_x2;position_matrix_x3,position_matrix_x4];
% position_matrix_y= [position_matrix_y1,position_matrix_y2;position_matrix_y3,position_matrix_y4];

%
%% receiver configuration, point, or plane, or objects
% receiver plan
incre_space = 0.05;
[X, Y] = meshgrid(0:incre_space:room.length, 0:incre_space:room.width);
[r,c] = size(X);
illumination = zeros(r,c);
receiver_height = 0.85;

%% Luminous distribution calculation
m = length(led.position_x);
n = length(led.position_y);
for i = 1:r
    for j = 1:c
        for k = 1:m
            for l = 1:n
                vector_A  = [X(i,j)-led.position_matrix_x(k,l), Y(i,j)-led.position_matrix_y(k,l), receiver_height-room.height];
                cos_theta = abs(vector_A(3))/norm(vector_A);
                cos_yita = cos_theta;
                gain = (led.lambert_index+1)/2/pi*cos_theta^(led.lambert_index)*cos_yita/norm(vector_A)^2;
                illumination(i,j) = illumination(i,j) + led.flux*gain;
%                 res(i,j)
            end
        end
    end
end

surf(X,Y,illumination);
%image(illumination)
colormap gray
colorbar
%figure(2);
%contourf(X,Y,illumination);