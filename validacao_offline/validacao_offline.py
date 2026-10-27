"""Validacao offline (sem Earth Engine) dos resultados publicados da v3.
Reproduz: auditoria dos 417 segmentos ANP, distancias reais cluster->duto, checagem dos rotulos de estado.
Executar: python validacao_offline.py  (a partir desta pasta)"""
import os, sys, warnings
import numpy as np, pandas as pd, geopandas as gpd
from shapely.ops import transform, nearest_points
from shapely.geometry import Point
from pyproj import Geod
warnings.filterwarnings("ignore")
PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "anp"))
V3 = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "v3_tabelas_publicadas"))
SP = os.path.dirname(os.path.abspath(__file__))
geod = Geod(ellps="WGS84")

GJ = os.path.join(SP, "..", "data", "limites", "brazil-states.geojson")
if not os.path.exists(GJ):
    import urllib.request
    urllib.request.urlretrieve("https://raw.githubusercontent.com/codeforamerica/click_that_hood/master/public/data/brazil-states.geojson", GJ)
states = gpd.read_file(GJ).to_crs(4326)
NE = ['MA', 'PI', 'CE', 'RN', 'PB', 'PE', 'AL', 'SE', 'BA']
ne = states[states.sigla.isin(NE)].reset_index(drop=True)
ne_union = ne.union_all()
study_limit = ne_union.buffer(1.0)
brazil = states.union_all()

files = {'gasoduto-escoamento.shp': ('Gathering Pipeline', 'Nome_Dut_1'),
         'gasoduto-distribuicao.shp': ('Distribution Pipeline', 'Distrib'),
         'gasoduto-transporte.shp': ('Transmission Pipeline', 'Nome_Dut_1'),
         'oleoduto.shp': ('Oil Pipeline', 'Nome_Dut_1')}
parts = []
for fn, (cls, col) in files.items():
    g = gpd.read_file(os.path.join(PROJ, fn)).to_crs(4326)
    g['ASSET_NAME'] = g[col].fillna('Not Identified'); g['CLASS'] = cls; g['SRC'] = fn
    g = g[['ASSET_NAME', 'CLASS', 'SRC', 'geometry']]
    g = g[g.intersects(study_limit)].copy()
    g['geometry'] = g.geometry.map(lambda geom: transform(lambda x, y, z=None: (x, y), geom))
    parts.append(g)
infra = pd.concat(parts, ignore_index=True)
infra = gpd.GeoDataFrame(infra, crs=4326)
print("== Legacy selection (1-degree buffer):", len(infra), infra.CLASS.value_counts().to_dict())

# Classify each segment: fraction of length onshore inside NE, offshore, other states
metric = infra.to_crs(5880)  # SIRGAS 2000 / Brazil Polyconic
ne_m = gpd.GeoSeries([ne_union], crs=4326).to_crs(5880).iloc[0]
br_m = gpd.GeoSeries([brazil], crs=4326).to_crs(5880).iloc[0]
L = metric.length
L_ne = metric.intersection(ne_m).length
L_land = metric.intersection(br_m).length
infra['len_km'] = L / 1000
infra['frac_in_NE_land'] = (L_ne / L).round(3)
infra['frac_offshore'] = (1 - L_land / L).round(3)
infra['frac_other_state'] = ((L_land - L_ne) / L).round(3)
def cat(r):
    if r.frac_in_NE_land >= 0.5: return 'onshore_NE'
    if r.frac_offshore >= 0.5: return 'offshore'
    return 'other_state'
infra['location'] = infra.apply(cat, axis=1)
print(pd.crosstab(infra.CLASS, infra.location, margins=True))
print("Length (km) by location:", infra.groupby('location').len_km.sum().round(0).to_dict())
oth = infra[infra.location == 'other_state']
print("Other-state segments by centroid lat (sample):", oth.geometry.centroid.y.describe().round(2).to_dict())
infra.drop(columns='geometry').to_csv(os.path.join(SP, 'auditoria_417_segmentos.csv'), index=False)

# ---- Clusters: 20 from Table4 + the GASALP one from Table2 ----
t4 = pd.read_csv(os.path.join(V3, 'Table4_Critical_Targets_Ranking.csv'), sep=';', encoding='utf-8-sig')
t2 = pd.read_csv(os.path.join(V3, 'Table2_Attribution_Impact_Matrix.csv'), sep=';', encoding='utf-8-sig')
cl = t4[['Rank', 'Latitude', 'Longitude', 'CH4_Peak_ppb']].copy()
for _, r in t2.iterrows():
    d = np.hypot(cl.Latitude - r.Latitude, cl.Longitude - r.Longitude)
    if d.min() > 1e-4:
        cl = pd.concat([cl, pd.DataFrame([{'Rank': 21, 'Latitude': r.Latitude, 'Longitude': r.Longitude, 'CH4_Peak_ppb': np.nan}])])
cl = cl.reset_index(drop=True)
pts = gpd.GeoDataFrame(cl, geometry=gpd.points_from_xy(cl.Longitude, cl.Latitude), crs=4326)
pts = gpd.sjoin(pts, ne[['sigla', 'geometry']], how='left', predicate='within').drop(columns='index_right')

def geodesic_nearest(pt, gdf):
    """Exact geodesic distance (km) from point to nearest line of gdf (candidate prefilter in metric CRS)."""
    best = (np.inf, None)
    pm = gpd.GeoSeries([pt], crs=4326).to_crs(5880).iloc[0]
    gm = gdf.to_crs(5880)
    dm = gm.distance(pm)
    for idx in dm.nsmallest(5).index:
        q = nearest_points(pt, gdf.geometry.loc[idx])[1]  # planar nearest in lon/lat, refine geodesically below
        # refine: densify candidate line and take geodesic min
        line = gdf.geometry.loc[idx]
        segs = [line] if line.geom_type == 'LineString' else list(line.geoms)
        dmin = np.inf
        for s in segs:
            s2 = s.segmentize(0.001)
            xs, ys = np.array(s2.coords).T[:2]
            _, _, d = geod.inv(np.full(xs.size, pt.x), np.full(ys.size, pt.y), xs, ys)
            dmin = min(dmin, d.min())
        if dmin < best[0]: best = (dmin, idx)
    return best[0] / 1000, best[1]

coast_line = gpd.GeoSeries([brazil.boundary], crs=4326)
rows = []
for i, r in pts.iterrows():
    d_all, j_all = geodesic_nearest(r.geometry, infra)
    on = infra[infra.location == 'onshore_NE']
    d_on, j_on = geodesic_nearest(r.geometry, on)
    pm = gpd.GeoSeries([r.geometry], crs=4326).to_crs(5880).iloc[0]
    d_coast = gpd.GeoSeries([brazil.boundary], crs=4326).to_crs(5880).distance(pm).iloc[0] / 1000
    rows.append(dict(Rank=int(r.Rank), Lat=round(r.Latitude, 4), Lon=round(r.Longitude, 4), UF_pip=r.sigla,
                     Peak=r.CH4_Peak_ppb, dist_nearest_any_km=round(d_all, 1),
                     nearest_any=f"{infra.loc[j_all,'ASSET_NAME']} [{infra.loc[j_all,'CLASS']}; {infra.loc[j_all,'location']}]",
                     dist_nearest_onshoreNE_km=round(d_on, 1), dist_coast_km=round(d_coast, 1)))
res = pd.DataFrame(rows).sort_values('Rank')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(res.to_string(index=False))
res.to_csv(os.path.join(SP, 'distancias_reais_clusters_dutos.csv'), index=False)
print("\nClusters (centroids) within 5/10/20 km of any pipeline:",
      [(k, int((res.dist_nearest_any_km <= k).sum())) for k in (5, 10, 20)])

# ---- State relabel validation ----
t1 = pd.read_csv(os.path.join(V3, 'Table1_Regional_Inventory_Final.csv'), sep=';', encoding='utf-8-sig')
wrong_to_idx = {'Maranhão': 0, 'Piauí': 1, 'Ceará': 2, 'Rio Grande do Norte': 3, 'Paraíba': 4, 'Pernambuco': 5,
                'Alagoas': 6, 'Sergipe': 7, 'Bahia': 8}
alpha = ['AL', 'BA', 'CE', 'MA', 'PB', 'PE', 'PI', 'RN', 'SE']
t1['UF_true'] = t1.State.map(wrong_to_idx).map(lambda k: alpha[k])
maxes = t1.set_index('UF_true').Max_ppb
wrong_uf = {'Maranhão': 'MA', 'Piauí': 'PI', 'Ceará': 'CE', 'Rio Grande do Norte': 'RN', 'Paraíba': 'PB',
            'Pernambuco': 'PE', 'Alagoas': 'AL', 'Sergipe': 'SE', 'Bahia': 'BA'}
max_wrong = t1.assign(UFw=t1.State.map(wrong_uf)).set_index('UFw').Max_ppb
ok_true = ok_wrong = n = 0
for _, r in res.dropna(subset=['UF_pip', 'Peak']).iterrows():
    n += 1
    ok_true += r.Peak <= maxes[r.UF_pip] + 1e-6
    ok_wrong += r.Peak <= max_wrong[r.UF_pip] + 1e-6
print(f"\nRelabel check (cluster peak <= its state's max): corrected labels {ok_true}/{n}, original labels {ok_wrong}/{n}")
viol = [(int(r.Rank), r.UF_pip, r.Peak, maxes[r.UF_pip]) for _, r in res.dropna(subset=['UF_pip','Peak']).iterrows() if r.Peak > maxes[r.UF_pip] + 1e-6]
print("Violations with corrected labels:", viol)
print(t1[['State', 'UF_true', 'Mean_ppb', 'Max_ppb', 'Std_Dev_ppb']].round(2).to_string(index=False))
