function [gamma] = normal_gravity_confocal(ityp,Ellipsoid,lat,height)
%NORMAL_GRAVITY_CONFOCAL  Normal gravity or gravitation of a reference ellipsoid.
%
%   gamma = NORMAL_GRAVITY_CONFOCAL(ityp, Ellipsoid, lat, height)
%
%   ityp = 0 for normal gravitation, 1 for normal gravity; lat in degrees,
%   geodetic height in metres; gamma in mGal.  Closed expressions of the
%   confocal-ellipsoid method are used, so the result is precise at any
%   geodetic height.
%
%Translated from NGA's Fortran by Theresa Damiani, NGS, November 2013
%Simplified for geodetic coordinate input only.
%
%From the FORTRAN subroutine "radgrav" in the hsynth_WGS84.f
%online codes. Written originally by Nikos K. Pavlis and Simon Holmes. 
%Available:
%http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%
%Documentation as follows is directly from the original subroutine, 
%modified only to remove mentions of non-geodetic coordinates:
%"C Subroutine to perform coordinate transformations on the meridian
% C plane and to compute the components of normal gravity/gravitation
% C vector (in ellipsoidal, geocentric, and geodetic coordinates),
% C and the magnitude of the normal gravity/gravitation vector.
% C
% C An equipotential ellipsoid of revolution and its associated
% C Somigliana-Pizzetti normal gravity field underlie all computations.
% C
% C The location of the computation point is input through its
% C geodetic latitude and geodetic (ellipsoidal) height.
% C Due to rotational symmetry, no longitude information is necessary.
% C
% C This subroutine uses closed expressions that yield precise results
% C regardless of the point's geodetic height."
%
%NOTE: GPS provides geodetic coordinates
%
%INPUT
%ityp output type: 0 for normal gravitation, 1 for normal gravity
%Ellipsoid structure ellipsoidal parameters
%lat latitude in degrees (geodetic)
%height ellipsoidal height in meters
%
%OUTPUT
%gamma full-field gravity/gravitation of ellipsoid at location (mGal)
%
%FUNCTION CALL
%For gravity (i.e. Mass and Rotational effects of Earth):
%[gamma]=normal_gravity(1,'wgs84',lat,h);
%
%For gravitation (i.e. Mass effect of Earth only):
%[gamma]=normal_gravity(0,'wgs84',lat,h);
%
%A disturbance at alitude is g_disturbance=ffg_measured-gamma;
%where ffg_measured is the output of NGS' Newton software and is the
%gravity value of our GRAV-D products released on the web.
%
%REFERENCES
%Heiskanen, W. and H. Mortiz (1967). Physical Geodesy. W.H. Freedman and
%Company: San Francisco. 354pp.
%
%Moritz, H. (2000). Geodetic Reference System 1980 (vol 66, pg 2, 1992). 
%[Correction]. Journal of Geodesy, 74(1), 128-133.
%
%NIMA (2004). Department of Defense World Geodetic System 1984: Its 
%Definition and Relationshops with Local Geodetic Systems. NIMA Technical 
%Report (Vol. Stock No. DMATR83502WGS84; NSN 7643-01-402-0347, pp. 175).

%% Ellipsoidal parameters
a=Ellipsoid.a;
b=Ellipsoid.b;
e=Ellipsoid.e;
e2=Ellipsoid.e2;
GM=Ellipsoid.GM;
omega=Ellipsoid.omega;

%% Other parameters
e21=1-(e2);

dtr=pi./180;
asqr=a.^2;
omega2=omega.^2;
Esqr=asqr.*e2;
E=sqrt(Esqr); %linear eccentricity, H&M p. 74
e_prime=E./b; %second eccentricity, H&M p. 74
ep2=E.^2./b.^2;

q0=0.5.*((1+3./(ep2)).*atan(e_prime)-(3./(e_prime))); %H&M p.66, 2-58

%Convert to radians
phi=dtr.*lat;

%% Calculate Coordinates psi and beta (geodetic coordinates)
sphi = sin(phi);
cphi = cos(phi);
sphi2 = sphi.*sphi;
dn=a./sqrt(1-e2.*sphi2);
p=(dn+height).*cphi;
z=(dn.*e21+height).*sphi;
z2=z.^2;

%Compute geocentric distance (m)
r=sqrt(p.*p+z2);

%Test for critical locations (poles/equator)
beta=atan(sqrt(e21).*sin(phi)./cos(phi));

%% Compute semi-minor axis of confocal ellipsoid (m)
dummy = r.*r-Esqr;
usqr = .5.*dummy.*(1+sqrt(1+4*Esqr.*z2./(dummy.*dummy)));
u = sqrt(usqr);

uE = usqr + Esqr;
uEi = 1./uE;
uErt = sqrt(uE);
uErti = 1./uErt;
uEb = usqr + Esqr.*sin(beta).^2;
uEbrt = sqrt(uEb);

Ei = 1./E;
uoE = u./E;
uoE2 = uoE.*uoE;
Eou = E./u;
atEu = atan(Eou);
const = omega2.*a.*a.*E./q0;
const2 = const.*Ei;
fnb = (0.5.*sin(beta).^2 - 1./6);

q = 0.5.*((1 + 3.*uoE2).*atEu - 3.*uoE);
q1 = 3*(1 + uoE2).*(1 - uoE.*atEu) - 1;
w = uEbrt./uErt;
wi = 1./w;
q1uEi = q1.*uEi;

fu = GM.*uEi + const.*fnb.*q1uEi - ityp.*omega2.*cos(beta).^2.*u;
gamu = 1E5.*(-wi.*fu);

fb = (-const2.*q.*uErti + ityp.*omega2.*uErt); 
gamb = 1E5.*(-wi.*fb.*cos(beta).*sin(beta));

%% Compute normal gravity/gravitation (mGal)
%Ellipsoidal System
gamma = sqrt(gamu.*gamu + gamb.*gamb);

end

