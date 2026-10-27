"""
Tests that use ONLY real inputs: the ANP pipeline layers in data/anp, the
pinned state boundaries in data/limites and the tables published with v3
(validacao_offline/v3_tabelas_publicadas, produced from real TROPOMI data).
No synthetic or invented data are used anywhere in this repository.

Run from the repository root:  python -m pytest tests -q
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pyproj import Geod
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ch4ne.boundaries import load_states, ne_states  # noqa: E402
from ch4ne.config import Config  # noqa: E402
from ch4ne.detection import nanmedian_filter  # noqa: E402
from ch4ne.infra import inventory_table, load_infrastructure, onshore_lines  # noqa: E402
from ch4ne.provenance import verify  # noqa: E402
from ch4ne.raster import distance_raster, distance_to_lines_km, grid_from_config, rasterize_states  # noqa: E402

V3 = ROOT / "validacao_offline" / "v3_tabelas_publicadas"
CFG = Config(base_dir=ROOT)


@pytest.fixture(scope="module")
def geo():
    states = load_states(CFG)
    ne = ne_states(states, CFG)
    infra = load_infrastructure(CFG, states, ne)
    lines = onshore_lines(infra, states, CFG)
    return states, ne, infra, lines


def test_states_come_from_attributes(geo):
    _, ne, _, _ = geo
    assert list(ne.UF) == sorted(CFG.ne_states)
    names = dict(zip(ne.UF, ne.state_name))
    assert names["MA"] == "Maranhão" and names["RN"] == "Rio Grande do Norte" and names["BA"] == "Bahia"


def test_inventory_reproduces_v3_selection(geo):
    """v3 integrated 417 segments; 262 are onshore NE, 127 offshore, 28 in other states."""
    _, _, infra, _ = geo
    inv = inventory_table(infra).set_index("infrastructure_type")
    assert int(inv.loc["Total", "segments_total"]) == 417
    assert int(inv.loc["Total", "segments_onshore_NE"]) == 262
    assert int(inv.loc["Total", "segments_offshore"]) == 127
    assert int(inv.loc["Total", "segments_other_state"]) == 28


def test_v3_state_labels_bug_and_fix(geo):
    """Published Table 1 used an IBGE-order index map on alphabetically ordered
    features. With the corrected labels 19/20 published cluster maxima fit under
    the maximum of the state that contains them (11/20 with the original labels)."""
    _, ne, _, _ = geo
    from shapely.geometry import Point
    t1 = pd.read_csv(V3 / "Table1_Regional_Inventory_Final.csv", sep=";", encoding="utf-8-sig")
    t4 = pd.read_csv(V3 / "Table4_Critical_Targets_Ranking.csv", sep=";", encoding="utf-8-sig")
    ibge_order = ["Maranhão", "Piauí", "Ceará", "Rio Grande do Norte", "Paraíba", "Pernambuco",
                  "Alagoas", "Sergipe", "Bahia"]
    alpha_uf = sorted(CFG.ne_states)                       # order of the source GeoJSON
    wrong_uf = dict(zip(ibge_order, ["MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA"]))
    true_uf = {name: alpha_uf[i] for i, name in enumerate(ibge_order)}
    mx_true = {true_uf[s]: m for s, m in zip(t1.State, t1.Max_ppb)}
    mx_wrong = {wrong_uf[s]: m for s, m in zip(t1.State, t1.Max_ppb)}
    ok_true = ok_wrong = 0
    for _, r in t4.iterrows():
        uf = ne[ne.contains(Point(r.Longitude, r.Latitude))].UF.iloc[0]
        ok_true += r.CH4_Peak_ppb <= mx_true[uf] + 1e-6
        ok_wrong += r.CH4_Peak_ppb <= mx_wrong[uf] + 1e-6
    assert (ok_true, ok_wrong) == (19, 11)


def test_v3_zero_distances_are_not_real(geo):
    """The four 'direct intersections (0 km)' of v3 are 4.9-19.2 km away from the
    nearest onshore pipeline; the pipeline's UTM distance agrees with an
    independent geodesic computation to better than 0.2 km."""
    _, _, _, lines = geo
    t2 = pd.read_csv(V3 / "Table2_Attribution_Impact_Matrix.csv", sep=";", encoding="utf-8-sig")
    assert (t2.Distance_km == 0).all()
    d, _ = distance_to_lines_km(lines, t2.Longitude.values, t2.Latitude.values)
    assert np.all(d > 4.0), d
    geod = Geod(ellps="WGS84")
    for (lon, lat), dk in zip(zip(t2.Longitude, t2.Latitude), d):
        best = np.inf
        for geom in lines.geometry:
            xs, ys = np.array(geom.segmentize(0.002).coords).T[:2]
            _, _, dist = geod.inv(np.full(xs.size, lon), np.full(ys.size, lat), xs, ys)
            best = min(best, dist.min() / 1000)
        assert abs(best - dk) < 0.2, (lon, lat, best, dk)
    assert sorted(np.round(d, 1)) == pytest.approx([4.9, 12.4, 14.0, 19.2], abs=0.3)


def test_nanmedian_filter_matches_scipy_on_real_field(geo):
    """Fast moving median equals scipy's generic_filter on a real field
    (distance to the ANP network over NE land; NaN over the sea)."""
    _, ne, _, lines = geo
    g = grid_from_config(CFG)
    inside = rasterize_states(ne, g) > 0
    d, _ = distance_raster(lines, g, inside)
    sub = d[120:200, 150:230]
    ref = ndimage.generic_filter(sub, np.nanmedian, size=9, mode="constant", cval=np.nan)
    got = nanmedian_filter(sub, 9)
    assert np.array_equal(np.isnan(ref), np.isnan(got))
    assert np.nanmax(np.abs(ref - got)) == 0.0


def test_analysis_refuses_non_earth_engine_inputs(tmp_path):
    with pytest.raises(RuntimeError, match="PROVENANCE"):
        verify(tmp_path)


def test_cell_areas_are_wgs84_ellipsoidal():
    """Analysis-grid cell areas equal pyproj's WGS84 geodesic polygon areas."""
    from shapely.geometry import box
    g = grid_from_config(CFG)
    geod = Geod(ellps="WGS84")
    A = g.pixel_area_km2
    for row in (0, 100, 200, g.height - 1):
        top = g.transform.f + row * g.transform.e
        ref = abs(geod.geometry_area_perimeter(box(-40, top - g.res, -40 + g.res, top))[0]) / 1e6
        assert abs(A[row, 0] / ref - 1) < 1e-6


def test_pipeline_lengths_are_geodesic(geo):
    """Segment lengths match pyproj's WGS84 geodesic length of the real ANP lines."""
    _, _, infra, _ = geo
    geod = Geod(ellps="WGS84")
    ref = infra.geometry.map(lambda g: geod.geometry_length(g) / 1000)
    assert np.allclose(infra["length_km"], ref, rtol=1e-9)
    inv = inventory_table(infra).set_index("infrastructure_type")
    assert inv.loc["Total", "length_on_NE_land_km"] == pytest.approx(9880.0, abs=0.5)

