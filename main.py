import math

# Constants for celestial bodies
CELESTIAL_BODIES = {
    "Mercury": {"G": 6.67430e-11, "M": 0.33011e24, "R": 2439.7e3},
    "Venus": {"G": 6.67430e-11, "M": 4.8675e24, "R": 6051.8e3},
    "Earth": {"G": 6.67430e-11, "M": 5.972e24, "R": 6371e3},
    "Moon": {"G": 6.67430e-11, "M": 0.07346e24, "R": 1737.4e3},
    "Mars": {"G": 6.67430e-11, "M": 0.64171e24, "R": 3389.5e3},
    "Jupiter": {"G": 6.67430e-11, "M": 1898.19e24, "R": 69911e3},
    "Saturn": {"G": 6.67430e-11, "M": 568.34e24, "R": 58232e3},
    "Uranus": {"G": 6.67430e-11, "M": 86.813e24, "R": 25362e3},
    "Neptune": {"G": 6.67430e-11, "M": 102.413e24, "R": 24622e3},
    "Sun": {"G": 6.67430e-11, "M": 1.989e30, "R": 696340e3}
}

AU_TO_METERS = 1.496e11  # 1 AU in meters

def orbital_velocity(G, M, R, altitude):
    """Calculate the orbital velocity at a given altitude above the celestial body's surface."""
    r = R + altitude
    v = math.sqrt(G * M / r)
    return v

def orbital_period(G, M, R, altitude):
    """Calculate the orbital period at a given altitude above the celestial body's surface."""
    r = R + altitude
    v = orbital_velocity(G, M, R, altitude)
    T = 2 * math.pi * r / v
    return T

def escape_velocity(G, M, R, altitude):
    """Calculate the escape velocity at a given altitude above the celestial body's surface."""
    r = R + altitude
    v_escape = math.sqrt(2 * G * M / r)
    return v_escape

def gravitational_force(G, M, R, mass, altitude):
    """Calculate the gravitational force on a mass at a given altitude above the celestial body's surface."""
    r = R + altitude
    F = G * M * mass / r**2
    return F

def potential_energy(G, M, R, mass, altitude):
    """Calculate the gravitational potential energy of a mass at a given altitude above the celestial body's surface."""
    r = R + altitude
    U = -G * M * mass / r
    return U

def specific_orbital_energy(G, M, R, altitude):
    """Calculate the specific orbital energy at a given altitude above the celestial body's surface."""
    r = R + altitude
    v = orbital_velocity(G, M, R, altitude)
    epsilon = v**2 / 2 - G * M / r
    return epsilon

def hohmann_transfer_delta_v(G, M, r1, r2):
    """Calculate the delta-v required for a Hohmann transfer between two circular orbits."""
    v1 = math.sqrt(G * M / r1)
    v2 = math.sqrt(G * M / r2)
    v_transfer1 = math.sqrt(2 * G * M * r2 / (r1 * (r1 + r2)))
    v_transfer2 = math.sqrt(2 * G * M * r1 / (r2 * (r1 + r2)))
    delta_v1 = v_transfer1 - v1
    delta_v2 = v2 - v_transfer2
    return delta_v1, delta_v2

def is_geostationary_orbit(altitude, celestial_body):
    """Check if a given altitude corresponds to a geostationary orbit for Earth."""
    if celestial_body == "Earth":
        GEOSTATIONARY_ALTITUDE = 35786e3  # Altitude of geostationary orbit for Earth, m
        return math.isclose(altitude, GEOSTATIONARY_ALTITUDE, rel_tol=1e-5)
    return False

def inclination_change_delta_v(v, delta_i):
    """Calculate the delta-v required to change the inclination of an orbit."""
    delta_v = 2 * v * math.sin(math.radians(delta_i) / 2)
    return delta_v

def semi_major_axis(R, altitude):
    """Calculate the semi-major axis of an orbit at a given altitude above the celestial body's surface."""
    return R + altitude

def eccentricity(R, perigee, apogee):
    """Calculate the eccentricity of an orbit given the perigee and apogee altitudes."""
    r_perigee = R + perigee
    r_apogee = R + apogee
    e = (r_apogee - r_perigee) / (r_apogee + r_perigee)
    return e

def true_anomaly(eccentricity, mean_anomaly):
    """Calculate the true anomaly from the mean anomaly and eccentricity."""
    E = mean_anomaly  # Initial guess for eccentric anomaly
    for _ in range(100):  # Iterate to solve Kepler's equation
        E = mean_anomaly + eccentricity * math.sin(E)
    theta = 2 * math.atan2(math.sqrt(1 + eccentricity) * math.sin(E / 2), math.sqrt(1 - eccentricity) * math.cos(E / 2))
    return math.degrees(theta)

def apoapsis_velocity(G, M, R, perigee, apogee):
    """Calculate the velocity at apoapsis for an elliptical orbit."""
    r_perigee = R + perigee
    r_apogee = R + apogee
    v_apogee = math.sqrt(G * M * (2 / r_apogee - 1 / ((r_perigee + r_apogee) / 2)))
    return v_apogee

def periapsis_velocity(G, M, R, perigee, apogee):
    """Calculate the velocity at periapsis for an elliptical orbit."""
    r_perigee = R + perigee
    r_apogee = R + apogee
    v_perigee = math.sqrt(G * M * (2 / r_perigee - 1 / ((r_perigee + r_apogee) / 2)))
    return v_perigee

def synodic_period(T1, T2):
    """Calculate the synodic period for two orbiting bodies."""
    return abs(1 / (1 / T1 - 1 / T2))

def keplerian_elements(R, perigee, apogee, inclination, raan, arg_perigee, true_anomaly):
    """Calculate the Keplerian orbital elements."""
    a = (R + perigee + R + apogee) / 2  # Semi-major axis
    e = eccentricity(R, perigee, apogee)  # Eccentricity
    i = inclination  # Inclination
    Ω = raan  # Right ascension of the ascending node
    ω = arg_perigee  # Argument of perigee
    ν = true_anomaly  # True anomaly
    return a, e, i, Ω, ω, ν

def main():
    try:
        print("Available celestial bodies: Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune, Sun")
        celestial_body = input("Choose the celestial body: ")
        if celestial_body not in CELESTIAL_BODIES:
            raise ValueError("Invalid celestial body choice.")
        
        G = CELESTIAL_BODIES[celestial_body]["G"]
        M = CELESTIAL_BODIES[celestial_body]["M"]
        R = CELESTIAL_BODIES[celestial_body]["R"]
        
        unit_choice = input("Choose units: (1) Meters and Kilograms (2) Kilometers and Kilograms (3) Astronomical Units and Kilograms: ")
        if unit_choice not in ['1', '2', '3']:
            raise ValueError("Invalid unit choice.")
        
        altitude = float(input("Enter the altitude above the celestial body's surface: "))
        mass = float(input("Enter the mass of the satellite: "))
        
        if unit_choice == '2':
            altitude *= 1000  # Convert kilometers to meters
        elif unit_choice == '3':
            altitude *= AU_TO_METERS  # Convert AU to meters
        
        if altitude < 0 or mass <= 0:
            raise ValueError("Altitude must be non-negative and mass must be positive.")
        
        v_orbital = orbital_velocity(G, M, R, altitude)
        T_orbital = orbital_period(G, M, R, altitude)
        v_escape = escape_velocity(G, M, R, altitude)
        F_gravitational = gravitational_force(G, M, R, mass, altitude)
        U_potential = potential_energy(G, M, R, mass, altitude)
        epsilon = specific_orbital_energy(G, M, R, altitude)
        a_semi_major = semi_major_axis(R, altitude)
        
        if unit_choice == '3':
            altitude_display = altitude / AU_TO_METERS
            unit_display = 'AU'
        else:
            altitude_display = altitude / 1000 if unit_choice == '2' else altitude
            unit_display = 'km' if unit_choice == '2' else 'm'
        
        print(f"Orbital velocity at {altitude_display} {unit_display}: {v_orbital:.2f} m/s")
        print(f"Orbital period at {altitude_display} {unit_display}: {T_orbital:.2f} seconds")
        print(f"Escape velocity at {altitude_display} {unit_display}: {v_escape:.2f} m/s")
        print(f"Gravitational force on a {mass} kg satellite at {altitude_display} {unit_display}: {F_gravitational:.2f} N")
        print(f"Gravitational potential energy of a {mass} kg satellite at {altitude_display} {unit_display}: {U_potential:.2f} J")
        print(f"Specific orbital energy at {altitude_display} {unit_display}: {epsilon:.2f} J/kg")
        print(f"Semi-major axis of the orbit at {altitude_display} {unit_display}: {a_semi_major / 1000 if unit_choice == '2' else a_semi_major:.2f} {'km' if unit_choice == '2' else 'm'}")
        
        if is_geostationary_orbit(altitude, celestial_body):
            print(f"The altitude {altitude_display} {unit_display} corresponds to a geostationary orbit.")
        else:
            print(f"The altitude {altitude_display} {unit_display} does not correspond to a geostationary orbit.")
        
        r1 = R + altitude
        r2 = R + (35786e3 if celestial_body == "Earth" else altitude)  # Use geostationary altitude for Earth
        delta_v1, delta_v2 = hohmann_transfer_delta_v(G, M, r1, r2)
        print(f"Delta-v for Hohmann transfer from {altitude_display} {unit_display} to geostationary orbit: {delta_v1:.2f} m/s and {delta_v2:.2f} m/s")
        
        delta_i = float(input("Enter the inclination change (in degrees): "))
        delta_v_inclination = inclination_change_delta_v(v_orbital, delta_i)
        print(f"Delta-v required for inclination change of {delta_i} degrees: {delta_v_inclination:.2f} m/s")
        
        perigee = float(input("Enter the perigee altitude: "))
        apogee = float(input("Enter the apogee altitude: "))
        
        if unit_choice == '2':
            perigee *= 1000  # Convert kilometers to meters
            apogee *= 1000  # Convert kilometers to meters
        elif unit_choice == '3':
            perigee *= AU_TO_METERS  # Convert AU to meters
            apogee *= AU_TO_METERS  # Convert AU to meters
        
        e_orbit = eccentricity(R, perigee, apogee)
        v_apogee = apoapsis_velocity(G, M, R, perigee, apogee)
        v_perigee = periapsis_velocity(G, M, R, perigee, apogee)
        print(f"Eccentricity of the orbit with perigee {perigee / 1000 if unit_choice == '2' else perigee} {'km' if unit_choice == '2' else 'm'} and apogee {apogee / 1000 if unit_choice == '2' else apogee} {'km' if unit_choice == '2' else 'm'}: {e_orbit:.2f}")
        print(f"Velocity at apoapsis: {v_apogee:.2f} m/s")
        print(f"Velocity at periapsis: {v_perigee:.2f} m/s")
        
        mean_anomaly = float(input("Enter the mean anomaly (in degrees): "))
        true_anomaly_value = true_anomaly(e_orbit, math.radians(mean_anomaly))
        print(f"True anomaly for mean anomaly {mean_anomaly} degrees and eccentricity {e_orbit:.2f}: {true_anomaly_value:.2f} degrees")
        
        # New functionality for Hohmann transfer to another orbit
        target_altitude = float(input("Enter the target altitude for Hohmann transfer: "))
        
        if unit_choice == '2':
            target_altitude *= 1000  # Convert kilometers to meters
        elif unit_choice == '3':
            target_altitude *= AU_TO_METERS  # Convert AU to meters
        
        if target_altitude < 0:
            raise ValueError("Target altitude must be non-negative.")
        
        r_target = R + target_altitude
        delta_v1, delta_v2 = hohmann_transfer_delta_v(G, M, r1, r_target)
        print(f"Delta-v for Hohmann transfer from {altitude_display} {unit_display} to {target_altitude / 1000 if unit_choice == '2' else target_altitude} {'km' if unit_choice == '2' else 'm'}: {delta_v1:.2f} m/s and {delta_v2:.2f} m/s")
    
    except ValueError as e:
        print(f"Invalid input: {e}")

if __name__ == "__main__":
    main()
    def orbital_inclination(launch_latitude, target_inclination):
        """Calculate the orbital inclination based on launch latitude and target inclination."""
        if abs(launch_latitude) > abs(target_inclination):
            raise ValueError("Launch latitude cannot be greater than target inclination.")
        inclination = math.degrees(math.acos(math.cos(math.radians(target_inclination)) / math.cos(math.radians(launch_latitude))))
        return inclination

    # Add functionality to calculate orbital inclination based on launch latitude
    launch_latitude = float(input("Enter the launch latitude (in degrees): "))
    target_inclination = float(input("Enter the target orbital inclination (in degrees): "))

    try:
        inclination = orbital_inclination(launch_latitude, target_inclination)
        print(f"Orbital inclination based on launch latitude {launch_latitude} degrees and target inclination {target_inclination} degrees: {inclination:.2f} degrees")
    except ValueError as e:
        print(f"Error calculating orbital inclination: {e}")
