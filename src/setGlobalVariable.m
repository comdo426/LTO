function setGlobalVariable 
% setGlobalVariable sets gravitation parameters and semimajor axis for
% major planets of interest.
% Unit should be km^3/s^2(mu), km(SMA)

    global MU_SUN
    MU_SUN = 1.32712440018e11;
    global MU_EARTH
    MU_EARTH = 3.986004418e5;
    global MU_MOON 
    MU_MOON = 4.9048695e3;
    global MU_MARS
    MU_MARS = 4.282837e4;
    global SMA_EARTHMOON
    SMA_EARTHMOON = 384400;
    global SMA_SUNEARTH
    SMA_SUNEARTH = 149649952;
    global SMA_SUNMARS
    SMA_SUNMARS = 227953016;
    global SMA_MARSDEIMOS
end