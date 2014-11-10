function gain = los_gain(Tx_cordinates, Rx_cordinates,half_angle,area,fov)
    d_square = sum((Tx_cordinates-Rx_cordinates).^2);
    lambertian_order = -log(2)/log(cos(half_angle*pi/180));
   
    incident_cos = abs(Tx_cordinates(3)-Rx_cordinates(3))/d_square^0.5;
    if incident_cos <= cos(fov/2*pi/180)
        visibility = 0;
    else
        visibility = 1;
    end
    
    if visibility == 0
        gain = 0;
    else
        gain = (lambertian_order+1)/2/pi*incident_cos^lambertian_order*area*incident_cos/d_square;
    end
end
    
    