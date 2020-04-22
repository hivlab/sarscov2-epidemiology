import seaborn as sns
from colour import Color
import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim


metadata = pd.read_csv(snakemake.input[0], sep="\t")
countries = list(set(metadata["country"]))
print(countries)
countries = ["Unknown" if i is np.nan else i for i in countries]
countires = countries.sort()
countries = [i.lower() for i in countries]
countries = countries.append(countries.pop(countries.index("Unknown")))
n = len(countries)

# Setup color palette
rgb_pal = sns.hls_palette(n, l=0.4)
hex_pal = [Color(rgb=i).hex for i in rgb_pal]

# Write colors to file
df = pd.DataFrame(list(zip(["country"] * n, countries, hex_pal)))
df.to_csv(snakemake.output.col, sep="\t", index=False, header=False)

# Get country coordinates
geolocator = Nominatim(user_agent="sars-cov-2")
locations = []
for country in countries:
    location = geolocator.geocode(country)
    locations.append(("country", country, location.latitude, location.longitude))

df = pd.DataFrame(locations)
df.to_csv(snakemake.output.loc, sep="\t", index=False, header=False)
