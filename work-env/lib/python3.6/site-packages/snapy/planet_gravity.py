from numpy import sin,sqrt,pi

def jupiter_gravity(pclat):
    grav = 23.3;
    S = sin(pclat * pi / 180.);
    SS = S*S;
    CS = S*sqrt(1 - SS);
    GR = - grav + SS*(-4.26594 + SS*(0.47685 + SS*(-0.100513 + SS*(0.0237067 - 0.00305515*SS))));
    GTH = CS*(-3.42313 +  SS*(0.119119 + SS*(0.00533106 + SS*(-0.00647658 + SS*0.000785945))));
    return sqrt(GR*GR+GTH*GTH);

def saturn_gravity(pclat):
    S = sin(pclat * pi / 180.);
    SS = S * S;
    CS = S*sqrt(1 - SS);
    GR = -9.06656 + SS*(-3.59253 + SS*(0.704538 + SS*(-0.260158 + SS*(0.0923098 - SS*0.0166287))));
    GTH = CS*(-2.25384 + SS*(.152112 + SS*(-.0102391 + SS*(-.00714765 + SS*.000865634))));
    return sqrt(GR*GR+GTH*GTH);


