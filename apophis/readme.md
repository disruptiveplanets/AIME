# Apophis

This folder contains code designed specifically to analyze the encounter of Apophis. I'm taking the known shape data of Apophis from [Lee et al 2022](https://arxiv.org/pdf/2204.02540.pdf), which includes Apophis as non-tumbling. So I have to alter my simulation to take this possibility into account.

## Data

Earth-centric orbital elements at perigee [JPL Ephemeris](https://ssd.jpl.nasa.gov/horizons/app.html#/):
    EC= 4.253901873212790E+00 QR= 3.801148649129103E+04 IN= 1.628561543134850E+02
    OM= 1.520019490904037E+02 W = 3.256762207107416E+01 Tp=  2462240.407091934700
    N = 2.865008496911329E-02 MA=-3.650928506076285E-01 TA= 3.598574271262584E+02
    A =-1.168181708373394E+04 AD= 9.999999999999998E+99 PR= 9.999999999999998E+99
    JDTDB    Julian Day Number, Barycentric Dynamical Time
      EC     Eccentricity, e
      QR     Periapsis distance, q (km)
      IN     Inclination w.r.t X-Y plane, i (degrees)
      OM     Longitude of Ascending Node, OMEGA, (degrees)
      W      Argument of Perifocus, w (degrees)
      Tp     Time of periapsis (Julian Day Number)
      N      Mean motion, n (degrees/sec)
      MA     Mean anomaly, M (degrees)
      TA     True anomaly, nu (degrees)
      A      Semi-major axis, a (km)
      AD     Apoapsis distance (km)
      PR     Sidereal orbit period (sec)

Physical properties [Lee et al 2022](https://arxiv.org/pdf/2204.02540.pdf):
    Iabc=0.64, 0.97, 1
    Klm=0.0316092, -0.0747126
    L in ecliptic coords: (275, -85 deg)
    rotation period: 264.178 hr
    precession period: 27.38547 hr

Location of North [Wikipedia](https://en.wikipedia.org/wiki/Orbital_pole)

Example of Apophis encounter map: https://arxiv.org/pdf/2111.08144.pdf