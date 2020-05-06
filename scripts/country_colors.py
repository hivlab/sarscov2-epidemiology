import seaborn as sns
from colour import Color
import pandas as pd
import geopy.geocoders
from geopy.geocoders import Nominatim


metadata = pd.read_csv(snakemake.input[0], sep="\t")
metadata_na = metadata.fillna(value = {"country": "unknown"})
countries = list(set(metadata_na["country"]))
countries = [i.lower() for i in countries]

# Get country coordinates
geopy.geocoders.options.default_timeout = None
geolocator = Nominatim(user_agent="sarscov2")
locations = []
for country in countries:
    location = geolocator.geocode(country)
    locations.append(("country", country, location.latitude, location.longitude))

loc = pd.DataFrame(locations)
loc.columns = ["key", "country", "lat", "lon"]
loc_sorted = loc.sort_values(by=["lon", "lat"])
loc_sorted.to_csv(snakemake.output.loc, sep="\t", index=False, header=False)

# Setup color palette
# Color Estonia black
countries_sorted = list(loc_sorted["country"])
n = len(countries_sorted)
rgb_pal = sns.hls_palette(n, l=0.4)
hex_pal = [Color(rgb=i).hex for i in rgb_pal]
est_black = ["#000000" if x == "estonia" else y for x, y in zip(countries_sorted, hex_pal)]

# Write contries and colors to file
cols = pd.DataFrame(list(zip(["country"] * n, countries_sorted, est_black)))
cols.to_csv(snakemake.output.col, sep="\t", index=False, header=False)

