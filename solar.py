import requests
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import webbrowser
import simplekml
import time

def get_geographical_coordinates(city_name):
    """
    Fetches geographical coordinates (latitude and longitude) for a given city name
    using the OpenStreetMap Nominatim API.
    
    Args:
        city_name (str): The name of the city.

    Returns:
        tuple: A tuple containing latitude and longitude as floats.
    """
    print("🌍 Getting geographical coordinates...")
    try:
        geo_api = f'https://nominatim.openstreetmap.org/search?q={city_name.replace(" ", "%20")}&format=json&limit=1'
        response = requests.get(geo_api)
        response.raise_for_status() # Raise an exception for bad status codes
        geo_data = response.json()
        
        if not geo_data:
            raise ValueError(f"City '{city_name}' not found.")
        
        latitude = float(geo_data[0]['lat'])
        longitude = float(geo_data[0]['lon'])
        
        print(f"📍 Location Found: {city_name} (Lat: {latitude:.4f}, Lon: {longitude:.4f})\n")
        return latitude, longitude
    except requests.exceptions.RequestException as e:
        print(f"Error fetching geographical coordinates: {e}")
        return None, None
    except (ValueError, IndexError) as e:
        print(f"Error processing city data: {e}")
        return None, None

def get_solar_radiation_data(latitude, longitude, start_date, end_date):
    """
    Fetches daily solar radiation data from the NASA POWER API.

    Args:
        latitude (float): Latitude of the location.
        longitude (float): Longitude of the location.
        start_date (str): Start date in 'YYYYMMDD' format.
        end_date (str): End date in 'YYYYMMDD' format.
    
    Returns:
        tuple: A tuple containing a list of dates and a list of solar values.
    """
    print("☀️ Fetching solar radiation data...")
    try:
        api_url = (
            'https://power.larc.nasa.gov/api/temporal/daily/point?'
            f'parameters=ALLSKY_SFC_SW_DWN&community=RE&longitude={longitude:.4f}&latitude={latitude:.4f}'
            f'&start={start_date}&end={end_date}&format=JSON'
        )
        response = requests.get(api_url)
        response.raise_for_status()
        data = response.json()
        
        solar_radiation = data['properties']['parameter']['ALLSKY_SFC_SW_DWN']
        dates = list(solar_radiation.keys())
        solar_values = list(solar_radiation.values())

        if not dates or not solar_values:
            print("No solar data found for the specified dates.")
            return [], []

        # Convert date strings (e.g., '20250101') to datetime objects
        date_array = [datetime.strptime(d, '%Y%m%d') for d in dates]

        return date_array, np.array(solar_values)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching solar radiation data: {e}")
        return [], []
    except (KeyError, TypeError) as e:
        print(f"Error parsing NASA POWER API response: {e}")
        return [], []
    except ValueError as e:
        print(f"Error converting date format: {e}")
        return [], []

def plot_solar_radiation(date_array, solar_values, city_name):
    """
    Plots past solar radiation data.
    """
    plt.figure(figsize=(12, 6))
    plt.plot(date_array, solar_values, '-o', linewidth=2)
    plt.xlabel('Date')
    plt.ylabel('Solar Radiation (kWh/m²/day)')
    plt.title(f'Past Solar Radiation Data for {city_name} (NASA POWER API)')
    plt.grid(True)
    plt.show()

def shading_analysis(date_array, solar_values, city_name, shade_threshold=0.3):
    """
    Performs a basic shading analysis and plots the results.
    """
    normalized_values = solar_values / np.max(solar_values)
    shaded_days = [date_array[i] for i, value in enumerate(normalized_values) if value < shade_threshold]
    
    if shaded_days:
        print('\n⚠️ Warning: Shading detected on these days:')
        for day in shaded_days:
            print(day.strftime('%Y-%m-%d'))
    else:
        print('\n✅ No major shading issues detected!')

    plt.figure(figsize=(12, 6))
    plt.bar(date_array, normalized_values, color='skyblue', label='Normalized Radiation')
    
    # Plot shaded days in a different color
    if shaded_days:
        shaded_values = [normalized_values[i] for i, value in enumerate(normalized_values) if value < shade_threshold]
        plt.bar(shaded_days, shaded_values, color='red', label='Shaded Days')
        
    plt.xlabel('Date')
    plt.ylabel('Normalized Solar Radiation')
    plt.title(f'Past Shading Analysis for {city_name} (Red = Shaded Days)')
    plt.grid(True)
    plt.legend()
    plt.show()

def get_elevation_data(latitude, longitude, grid_size=5):
    """
    Fetches elevation data for a grid of points around a given location.
    
    Returns:
        tuple: A tuple containing a grid of longitudes, latitudes, and elevation data.
    """
    print("\nFetching Elevation Data...")
    try:
        latitudes = np.linspace(latitude - 0.01, latitude + 0.01, grid_size)
        longitudes = np.linspace(longitude - 0.01, longitude + 0.01, grid_size)
        lon_grid, lat_grid = np.meshgrid(longitudes, latitudes)
        
        locations = ';'.join([f"{lat:.6f},{lon:.6f}" for lat, lon in zip(lat_grid.flatten(), lon_grid.flatten())])
        api_url = f'https://api.open-elevation.com/api/v1/lookup?locations={locations}'
        
        response = requests.get(api_url, timeout=20)
        response.raise_for_status()
        elev_response = response.json()
        
        elevation_data = np.array([result['elevation'] for result in elev_response['results']]).reshape(grid_size, grid_size)
        
        print('✅ Elevation Data Fetched Successfully!')
        return lon_grid, lat_grid, elevation_data
    except requests.exceptions.RequestException as e:
        print(f"Error fetching elevation data: {e}")
        return None, None, None
    except (KeyError, TypeError) as e:
        print(f"Error parsing elevation API response: {e}")
        return None, None, None

def plot_3d_terrain(lon_grid, lat_grid, elevation_data, city_name):
    """
    Generates a 3D terrain plot.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(lon_grid, lat_grid, elevation_data, cmap='jet', edgecolor='none')
    
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Elevation (m)')
    ax.set_title(f'3D Terrain Model for {city_name}')
    ax.grid(True)
    plt.show()

def generate_kml_file(current_lat, current_lon, elevation, file_name='SolarPanel_3D_Map.kml'):
    """
    Generates a KML file for Google Earth.
    """
    print(f"\nCreating KML file '{file_name}' for Google Earth...")
    kml = simplekml.Kml()
    pnt = kml.newpoint(name='Solar Panel')
    pnt.coords = [(current_lon, current_lat, elevation)] # (lon, lat, alt)
    kml.save(file_name)
    print(f"✅ KML file '{file_name}' generated successfully!")

def main():
    """
    Main function to run the solar analysis script.
    """
    # 1. Prompt for city name
    city_name = input('Enter city name: ')
    
    # 2. Get geographical coordinates
    latitude, longitude = get_geographical_coordinates(city_name)
    if latitude is None:
        return
    
    # 3. Fetch past solar radiation data
    start_date = input('Enter the start date for past data (YYYYMMDD): ')
    end_date = input('Enter the end date for past data (YYYYMMDD): ')
    date_array, solar_values = get_solar_radiation_data(latitude, longitude, start_date, end_date)
    if not date_array:
        return

    # 4. Plot solar radiation and shading analysis
    plot_solar_radiation(date_array, solar_values, city_name)
    shading_analysis(date_array, solar_values, city_name)
    
    # 5. Ask for corrected coordinates and implementation date
    print('\n🌐 Please enter more precise coordinates from Google Maps.')
    webbrowser.open('https://www.google.com/maps')
    current_lat = float(input('Enter your implementation place latitude: '))
    current_lon = float(input('Enter your implementation place longitude: '))
    print(f'📍 Using Selected Location: Lat: {current_lat:.15f}, Lon: {current_lon:.15f}\n')

    impl_date = input('Enter the date of implementation (YYYYMMDD): ')
    impl_date_obj = datetime.strptime(impl_date, '%Y%m%d')
    
    # Fetch data for the implementation date
    impl_dates, impl_values = get_solar_radiation_data(current_lat, current_lon, impl_date, impl_date)
    if not impl_dates:
        return
    
    # 6. Predict best location (basic heuristic)
    average_radiation = np.mean(impl_values)
    if average_radiation > 2.0: # Using a kWh/m2/day value as a more realistic threshold
        print('✅ The selected location is suitable for solar panel implementation!')
    else:
        print('⚠️ Warning: The selected location may not be optimal for solar panel implementation. Consider alternative locations.')

    # 7. Solar panel calculations
    required_energy = float(input('Enter the required energy output in kWh: '))
    use_rotation = input('Do you want rotatable solar panels? (y/n): ').lower() == 'y'

    panel_efficiency = 0.2
    panel_area = 1.6 # in m^2

    # Calculate required number of panels
    energy_per_panel = np.mean(impl_values) * panel_efficiency * panel_area
    num_panels = np.ceil(required_energy / energy_per_panel)

    print(f'\n🔋 Estimated Number of Solar Panels Required: {int(num_panels)}')
    print(f'⚡ Estimated Energy Output per Panel: {energy_per_panel:.2f} kWh/day')
    print(f'⚡ Total Estimated Energy Output: {num_panels * energy_per_panel:.2f} kWh/day')

    # 8. Optimal tilt and rotation
    # MATLAB's sind, cosd etc use degrees, so we convert here.
    
    # Calculate optimal tilt angle
    day_of_year = impl_date_obj.timetuple().tm_yday
    declination = 23.45 * np.sin(np.deg2rad(360 * (day_of_year - 81) / 365))
    optimal_tilt_angle = np.abs(np.rad2deg(np.deg2rad(current_lat) - np.deg2rad(declination)))

    print(f'\n📐 Optimal Tilt Angle for Maximum Energy Generation: {optimal_tilt_angle:.2f}°')

    if use_rotation:
        print('\n✅ Rotatable Solar Panels Selected! Calculating optimal rotation angles...')
        
        solar_time = np.linspace(6, 18, 100)
        hour_angle = (solar_time - 12) * 15
        
        # Calculate altitude angle and optimal rotation angle throughout the day
        altitude_angle = np.rad2deg(
            np.arcsin(
                np.sin(np.deg2rad(current_lat)) * np.sin(np.deg2rad(declination)) +
                np.cos(np.deg2rad(current_lat)) * np.cos(np.deg2rad(declination)) * np.cos(np.deg2rad(hour_angle))
            )
        )
        optimal_rotation_angles = np.rad2deg(np.arctan2(np.sin(np.deg2rad(hour_angle)), 
                                                   np.cos(np.deg2rad(hour_angle)) * np.sin(np.deg2rad(current_lat))))
        
        print('\n✅ Optimal Panel Rotation Angles Calculated!')
        for i, t in enumerate(solar_time[::10]): # Print every 10th value for brevity
            print(f'⏳ At {t:.2f} Hours, Rotate Panel to {optimal_rotation_angles[i*10]:.2f}°')
            
        optimal_tilt = np.max(altitude_angle)
        print(f'\n📐 Maximum Optimal Tilt Angle of the Day: {optimal_tilt:.2f}°')
        
        # Plot optimal rotation angles
        plt.figure(figsize=(10, 6))
        plt.plot(solar_time, optimal_rotation_angles, '-o', linewidth=2)
        plt.xlabel('Time (Hours)')
        plt.ylabel('Optimal Rotation Angle (Degrees)')
        plt.title('Optimal Panel Rotation Angle Throughout the Day')
        plt.grid(True)
        plt.show()

    # 9. 3D Terrain Visualization
    lon_grid, lat_grid, elevation_data = get_elevation_data(current_lat, current_lon)
    if lon_grid is not None:
        plot_3d_terrain(lon_grid, lat_grid, elevation_data, city_name)
        
    # 10. KML file generation
    # Fetch a single elevation point for the KML file
    try:
        elev_api = f'https://api.open-elevation.com/api/v1/lookup?locations={current_lat:.6f},{current_lon:.6f}'
        response = requests.get(elev_api, timeout=20)
        response.raise_for_status()
        elev_data = response.json()
        elevation = elev_data['results'][0]['elevation']
        generate_kml_file(current_lat, current_lon, elevation)
    except Exception as e:
        print(f"Could not get elevation for KML file: {e}")

if __name__ == '__main__':
    main()
