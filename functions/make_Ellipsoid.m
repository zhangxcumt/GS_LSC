function [ellipsoid] = make_Ellipsoid(ellipType)
    %MAKE_ELLIPSOID  Parameters of the WGS84 or GRS80 reference ellipsoid.
    %
    %   ellipsoid = MAKE_ELLIPSOID('WGS84' | 'GRS80')
    %
    %   Fields:
    %     GM           geocentric gravitational constant (m^3/s^2)
    %     omega        angular velocity (rad/s)
    %     gammaEquator normal gravity at the equator (mGal)
    %     gammaPole    normal gravity at the pole (mGal)
    %     a, b         semi-major and semi-minor axis (m)
    %     e            first eccentricity
    %     e2           first eccentricity squared
    %     f            flattening
    %
    %   WGS84 from NIMA (2004), GRS80 from Moritz (2000).

    switch ellipType
        case 'WGS84'
            ellipsoid.GM = 3986004.418E8;
            ellipsoid.omega = 7.292115e-5;
            ellipsoid.gammaEquator = 978032.53359;
            ellipsoid.gammaPole = 983218.49378;
            ellipsoid.a = 6378137.0;
            ellipsoid.b = 6356752.3142;
            ellipsoid.e = 8.1819190842622 * 1E-2; % first eccentricity
            ellipsoid.e2 = 6.69437999014 * 1E-3; % first eccentricity squared
            ellipsoid.f = 1 / 298.257223563; % flattening
            % WGS84 defining parameters, NIMA (2004).
        case 'GRS80'
            ellipsoid.GM = 3986005E8;
            ellipsoid.omega = 7.292115e-5;
            ellipsoid.gammaEquator = 978032.67715;
            ellipsoid.gammaPole = 983218.63685;
            ellipsoid.a = 6378137;
            ellipsoid.b = 6356752.3142;
            ellipsoid.e = 0.081819190842622; % first eccentricity
            ellipsoid.e2 = 0.00669437999014; % first eccentricity squared
            ellipsoid.f = 0.00335281066474; % flattening
            % GRS80 defining parameters, Moritz (2000).
        otherwise
            warning('This ellipsoid is not supported')
            return
    end
end
