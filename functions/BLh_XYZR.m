function XYZR = BLh_XYZR(BLh)
    %BLH_XYZR  Geodetic to geocentric cartesian coordinates.
    %
    %   XYZR = BLH_XYZR(BLh) with BLh = [latitude longitude height] in
    %   [deg deg m], returns [X Y Z r] with r = sqrt(X^2 + Y^2 + Z^2).
    %
    %   Ellipsoid: GRS80 (a = 6378137 m, 1/f = 298.257222101).
    a = 6378137;
    f = 1/298.257222101;
    e2 = (2 - f) * f;
    h = BLh(:,3);

    slat = sind(BLh(:,1));
    clat = cosd(BLh(:,1));
    slon = sind(BLh(:,2));
    clon = cosd(BLh(:,2));

    v = a ./ sqrt(1 - e2 * slat.^2);

    XYZR(:,1) = (v + h) .* clat .* clon;
    XYZR(:,2) = (v + h) .* clat .* slon;
    XYZR(:,3) = (v * (1 - e2) + h) .* slat;
    XYZR(:,4) = sqrt(XYZR(:,1).^2 + XYZR(:,2).^2 + XYZR(:,3).^2);
end
