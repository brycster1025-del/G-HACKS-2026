import novatel_edie as ne
import numpy as np

np.set_printoptions(precision=15)

SEMI_MAJOR = 6378137
INV_FLAT = 298.257222101

flat = 1 / INV_FLAT
square_eccentricity = flat * (2 - flat)

def geo_to_cartesian(coords):
    lat = np.radians(coords[0])
    long = np.radians(coords[1])
    height = coords[2]
    N = SEMI_MAJOR / np.sqrt(1 - square_eccentricity * np.sin(lat)**2)
    x = (N + height) * np.cos(lat) * np.cos(long)
    y = (N + height) * np.cos(lat) * np.sin(long)
    z = ((1 - square_eccentricity) * N + height) * np.sin(lat)
    return np.array([x, y, z])

def cartesian_to_geo(coords, precision=10**-20):
    x = coords[0]
    y = coords[1]
    z = coords[2]

    long = np.arctan2(y, x)
    p = np.sqrt(x**2 + y**2)

    lat = np.arctan2(z, ((1 - square_eccentricity) * p))
    lat_prev = 0

    while np.abs(lat - lat_prev) > precision:
        lat_prev = lat
        N = SEMI_MAJOR / np.sqrt(1 - square_eccentricity * np.sin(lat_prev)**2)
        h = p / np.cos(lat_prev) - N
        lat = np.arctan2(z, ((1 - square_eccentricity * N / (N + h)) * p))

    lat = np.degrees(lat)
    long = np.degrees(long)
    return np.array([lat, long, h])

def cartesian_to_local(origin, coords):
    lat = np.radians(origin[0])
    long = np.radians(origin[1])
    rot = np.array([[-np.sin(lat) * np.cos(long), -np.sin(lat) * np.sin(long), np.cos(lat)],
                        [-np.sin(long), np.cos(long), 0],
                        [np.cos(lat) * np.cos(long), np.cos(lat) * np.sin(long), np.sin(lat)]])
    return rot @ coords

def local_to_cartesian(origin, coords):
    lat = np.radians(origin[0])
    long = np.radians(origin[1])
    rot = np.array([[-np.sin(lat) * np.cos(long), -np.sin(lat) * np.sin(long), np.cos(lat)],
                        [-np.sin(long), np.cos(long), 0],
                        [np.cos(lat) * np.cos(long), np.cos(lat) * np.sin(long), np.sin(lat)]])
    return np.transpose(rot) @ coords

def local_change(bearing, distance, height_change):
    bearing = np.radians(bearing)
    north_change = distance * np.cos(bearing)
    east_change = distance * np.sin(bearing)
    change = np.array([north_change, east_change, height_change])
    return change

def similarity_transform(offset, coord1_legacy, coord1_current, coord2_legacy, coord2_current):
    offset_12_legacy = coord2_legacy - coord1_legacy
    offset_12_current = coord2_current - coord1_current
    bearing_12_legacy = np.arctan2(offset_12_legacy[1], offset_12_legacy[0])
    bearing_12_current = np.arctan2(offset_12_current[1], offset_12_current[0])
    theta = -(bearing_12_legacy - bearing_12_current)
    scale = np.sqrt((offset_12_current[0]**2 + offset_12_current[1]**2) / (offset_12_legacy[0]**2 + offset_12_legacy[1]**2))
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                   [np.sin(theta), np.cos(theta)]])
    return scale * rot @ offset

def circle_intersect(pos1, rad1, pos2, rad2):  # https://math.stackexchange.com/a/1033561
    distance = np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)
    length = (rad1**2 - rad2**2 + distance**2) / (2 * distance)
    height = np.sqrt(rad1**2 - length**2)
    parallel = (pos2 - pos1) / distance
    perpendicular = np.array([parallel[1], -parallel[0]])
    intersect1 = length * parallel + height * perpendicular + pos1
    intersect2 = length * perpendicular - height * perpendicular + pos1
    return intersect1, intersect2

geo_A = np.array([51.07899666397, -114.13251409834, 1114.3291 - 1.4557])

cartesian_A = geo_to_cartesian(geo_A)
local_B = local_change(153.371958, 128.5967, -3.710)
cartesian_AB = local_to_cartesian(geo_A, local_B)
cartesian_B = cartesian_A + cartesian_AB
geo_B = cartesian_to_geo(cartesian_B)

cartesian_BC = np.array([-74.7377, -69.2944, -76.4427])
cartesian_C = cartesian_B + cartesian_BC
geo_C = cartesian_to_geo(cartesian_C)

local_C_D = np.array([-131.4161, 15.3513, -1.5663])
cartesian_CD = local_to_cartesian(geo_C, local_C_D)
cartesian_D = cartesian_C + cartesian_CD
geo_D = cartesian_to_geo(cartesian_D)

bearing_DE = 246.102478
bearing_EP = 139.361144
height_change_DE = 1.1277
local_P = np.array([-505.2385, 36.3253, 0])
local_D = cartesian_to_local(geo_A, cartesian_D - cartesian_A)
local_DP = local_P - local_D
bearing_DP = np.degrees(np.arctan2(local_DP[1], local_DP[0]))
distance_DP = np.sqrt(local_DP[0]**2 + local_DP[1]**2)
distance_DE = distance_DP * np.sin(np.radians(bearing_DP - bearing_EP)) / np.sin(np.radians(bearing_EP - bearing_DE + 180))
local_DE = local_change(bearing_DE, distance_DE, height_change_DE)
local_E = local_D + local_DE
cartesian_AE = local_to_cartesian(geo_A, local_E)
cartesian_E = cartesian_A + cartesian_AE
geo_E = cartesian_to_geo(cartesian_E)

local_1_legacy = np.array([225, 150])
local_1_current = np.array([489.4949, 123.2447])
local_2_legacy = np.array([600, 450])
local_2_current = np.array([1063.8712, 170.0614])
local_EF_legacy = np.array([-151.2606, 33.6304])
height_change_EF = -3.6288
local_EF = np.append(similarity_transform(local_EF_legacy, local_1_legacy, local_1_current, local_2_legacy, local_2_current), height_change_EF)
local_F = local_E + local_EF
cartesian_AF = local_to_cartesian(geo_A, local_F)
cartesian_F = cartesian_A + cartesian_AF
geo_F = cartesian_to_geo(cartesian_F)

local_Gs = circle_intersect(local_F[:2], 126.6901, local_D[:2], 168.5144)
local_G1 = np.append(local_Gs[0], 3.8126)
local_G2 = np.append(local_Gs[1], 3.8126)

print(geo_A)
print(geo_B)
print(geo_C)
print(geo_D)
print(geo_E)
print(geo_F)

print(cartesian_to_geo(cartesian_A + local_to_cartesian(geo_A, local_G1)))
print(cartesian_to_geo(cartesian_A + local_to_cartesian(geo_A, local_G2)))
