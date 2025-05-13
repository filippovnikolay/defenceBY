import streamlit as st
import osmnx as ox
from shapely.geometry import shape, LineString, Polygon
from shapely.ops import unary_union
import pandas as pd
import geopandas as gpd
import folium
from folium.plugins import Draw
from streamlit_folium import st_folium

# Set page configuration
st.set_page_config(layout="wide")

left_col, right_col = st.columns([1, 1])

# Initialize session state variables to store the forest and lake data
if 'forest_lake_data' not in st.session_state:
    st.session_state.forest_lake_data = None
if 'last_polygon' not in st.session_state:
    st.session_state.last_polygon = None

# Define UTM CRS for Belarus (zone 35N)
utm_crs = "EPSG:32635"

# --- Function to normalize line direction (A→B == B→A) ---
def normalize_linestring(geom):
    if isinstance(geom, LineString):
        coords = list(geom.coords)
        rev_coords = list(reversed(coords))
        return min(LineString(coords).wkt, LineString(rev_coords).wkt)
    return geom.wkt

# --- Function to split large polygons into smaller chunks ---
def split_polygon(polygon, num_chunks=4):
    minx, miny, maxx, maxy = polygon.bounds
    width = (maxx - minx) / num_chunks
    height = (maxy - miny) / num_chunks

    polygons = []
    for i in range(num_chunks):
        for j in range(num_chunks):
            x_min = minx + i * width
            x_max = minx + (i + 1) * width
            y_min = miny + j * height
            y_max = miny + (j + 1) * height
            polygons.append(Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)]))

    return polygons

# Left column: drawing map
with left_col:
    m = folium.Map(location=[53.8989, 27.5592], zoom_start=14, control_scale=True)

    # Add draw plugin
    draw = Draw(
        export=False,
        draw_options={
            "polygon": {
                "shapeOptions": {
                    "color": "#3388ff"
                }
            },
            "polyline": {
                "shapeOptions": {
                    "color": "red",
                    "dashArray": "5, 10",
                    "weight": 3
                }
            },
            "circle": False,
            "rectangle": False,
            "marker": True,
            "circlemarker": False,
        }
    )
    draw.add_to(m)
    map_data = st_folium(m, width="100%", height=1024)

# Right column: results
with right_col:
    if map_data and map_data.get("all_drawings"):
        drawings = map_data["all_drawings"]
        polygon_geom = None
        line_geom = None

        for drawing in drawings:
            if drawing["geometry"]["type"] == "Polygon":
                polygon_geom = shape(drawing["geometry"])
            elif drawing["geometry"]["type"] == "LineString":
                line_geom = shape(drawing["geometry"])

        # Check if the polygon is different from the last one
        if polygon_geom != st.session_state.last_polygon:
            st.session_state.forest_lake_data = None  # Clear previous forest and lake data
            st.session_state.last_polygon = polygon_geom  # Save the new polygon

        if polygon_geom:
            st.subheader("🛣️ Анализ дорог")

            # Inside the 'with st.spinner' block where you're loading and processing road data
            with st.spinner("🔄 Загрузка данных дорог..."):
                try:
                    # If polygon is too large, split it into smaller chunks
                    if polygon_geom.area > 1000000:  # Threshold can be adjusted based on your needs
                        smaller_polygons = split_polygon(polygon_geom)
                        roads = []
                        for small_polygon in smaller_polygons:
                            try:
                                G = ox.graph_from_polygon(small_polygon, network_type="drive", simplify=True)
                                G_proj = ox.project_graph(G)
                                edges = ox.graph_to_gdfs(G_proj, nodes=False)
                                roads.append(edges)
                            except Exception as e:
                                continue
                        if roads:
                            edges = pd.concat(roads, ignore_index=True)
                        else:
                            raise ValueError("Дороги не обнаружены внутри выбранного полигона.")
                    else:
                        # Fetch roads for the current polygon
                        G = ox.graph_from_polygon(polygon_geom, network_type="drive", simplify=True)
                        G_proj = ox.project_graph(G)
                        edges = ox.graph_to_gdfs(G_proj, nodes=False)

                    if edges.empty:
                        raise ValueError("Дороги не обнаружены внутри выбранного полигона.")

                    if "ref" not in edges.columns:
                        edges["ref"] = ""

                    edges["ref"] = edges["ref"].astype(str).fillna('').str.strip()

                    # Filter roads where 'ref' is not empty and not 'nan'
                    filtered_roads = edges[edges["ref"].notna() & (edges["ref"] != "") & (edges["ref"] != "nan")]

                    if filtered_roads.empty:
                        st.warning("⚠️Индексированные дороги не найдены.")
                    else:
                        # Optional: Clip roads strictly within polygon
                        poly_proj = gpd.GeoSeries([polygon_geom], crs="EPSG:4326").to_crs(G_proj.graph['crs']).iloc[0]
                        filtered_roads = gpd.clip(filtered_roads, poly_proj)

                        # Normalize geometries to eliminate reverse duplicates
                        filtered_roads["normalized_wkt"] = filtered_roads.geometry.apply(normalize_linestring)
                        filtered_roads = filtered_roads.drop_duplicates(subset=["ref", "normalized_wkt"])

                        # Compute length in kilometers (first calculation)
                        filtered_roads["length_km"] = filtered_roads.geometry.length / 1000

                        # ⚠️ Adjust length for highways starting with "М" or "M" after length calculation
                        filtered_roads["adjusted_length_km"] = filtered_roads.apply(
                            lambda row: row["length_km"] / 2 if str(row["ref"]).strip().upper().startswith("М") or str(
                                row["ref"]).strip().upper().startswith("M")
                            else row["length_km"],
                            axis=1
                        )

                        # Summarize
                        summary_table = (
                            filtered_roads.groupby("ref")["adjusted_length_km"]
                                .sum()
                                .reset_index()
                                .rename(columns={"ref": "Наименование дороги", "adjusted_length_km": "Длина (км)"})
                        )

                        total_length = summary_table["Длина (км)"].sum()

                        summary_table = summary_table.sort_values(by="Длина (км)", ascending=False)

                        # Reset index and remove old one
                        summary_table = summary_table.reset_index(drop=True)
                        roads_count = summary_table["Длина (км)"].count()
                        # Display the table without the index column
                        st.dataframe(summary_table.style.hide(axis="index"), use_container_width=True)

                        col1, col2 = st.columns(2)
                        with col1:

                            st.markdown(r'$n_{дор} = $' + f" {round(roads_count, 0)}")
                        with col2:
                            st.markdown(r'$L_{дор} = $' + f" {round(total_length, 1)} км")


                except ValueError as e:
                    st.error(f"Попробуйте увеличить размер полигона или убедитесь, что в области есть дороги.")
                except Exception as e:
                    st.error(f"Неизвестная ошибка: {str(e)}")

            # Load forest and lake data for the drawn polygon (only once)
            if st.session_state.forest_lake_data is None:
                st.subheader("🌳 Леса и водоемы в области")
                with st.spinner("🔄 Загрузка данных о лесах и водоемах..."):
                    tags = {"natural": ["wood", "water"], "landuse": "forest"}
                    gdf = ox.features_from_polygon(polygon_geom, tags)
                    gdf = gdf.to_crs(utm_crs)
                    st.session_state.forest_lake_data = gdf
                    st.success("✅ Данные о лесах и водоемах загружены!")

        if line_geom and polygon_geom:
            st.subheader("🌳 Анализ протяженности танкодоступной местности")

            if st.session_state.forest_lake_data is not None:
                with st.spinner("🔄 Анализ протяженности танкодоступной местности..."):
                    line_proj = gpd.GeoSeries([line_geom], crs="EPSG:4326").to_crs(utm_crs).geometry[0]
                    poly_proj = gpd.GeoSeries([polygon_geom], crs="EPSG:4326").to_crs(utm_crs).geometry[0]

                    water_forest = st.session_state.forest_lake_data[st.session_state.forest_lake_data.geometry.intersects(poly_proj)]
                    combined_obstacles = unary_union(water_forest.geometry.values)

                    line_length = line_proj.length
                    line_km = line_length / 1000  # Convert meters to kilometers

                    accessible_geom = line_proj.difference(combined_obstacles)
                    accessible_km = accessible_geom.length / 1000  # Convert meters to kilometers

                    # Show lengths in two columns
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown(r'$L= $' + f" {round(line_km, 1)} км")
                    with col2:
                        st.markdown(r'$L_{ТД}= $' + f" {round(accessible_km, 1)} км")

                    st.subheader("📝 Расчетные Коэффициенты")

                    # Create two input columns
                    col1, col2 = st.columns(2)

                    with col1:
                        st.markdown(r'$N_{Т}^{Пр}$')
                        tanks_count = st.number_input("Количество танков противника", min_value=1, step=1,
                                                      label_visibility="collapsed")

                        st.markdown(r'$g_{M}$')
                        gm = st.number_input("Дм", min_value=0.05, max_value=0.08, step=0.01,
                                             label_visibility="collapsed")

                        st.markdown(r'$К_{ИБП}^{МВЗ}$')
                        Kmvz = st.number_input("Kмвз", min_value=0.75, max_value=1.00, step=0.01,
                                               label_visibility="collapsed")

                    with col2:
                        st.markdown(r'$K_{у}^{Пр}$')
                        k_y_pr = st.number_input("Kупр", value=0.95, label_visibility="collapsed", disabled=True)

                        st.markdown(r'$Н_{З}$')
                        Nz = st.number_input("Нз", min_value=0.75, max_value=1.00, step=0.01,
                                             label_visibility="collapsed")

                        st.markdown(r'$К_{ИБП}^{УЗ}$')
                        Kuz = st.number_input("Kуз", min_value=0.40, max_value=0.80, step=0.01,
                                              label_visibility="collapsed")


                    st.subheader("🛡️ Оценка Эффективности")

                    n_uz_treb = round(0.01*summary_table["Длина (км)"].count()*summary_table["Длина (км)"].sum(), 1)
                    st.markdown(r'$n_{уз}^{Треб}=0,01 * n_{дор} * L_{дор} = $' + " " + str(n_uz_treb))

                    N_ibp_treb = round((Kuz * n_uz_treb + Kmvz * round(accessible_km, 1)) * 1.1, 2)
                    st.markdown(r'$N_{ИБП}^{Треб}=(К_{ИБП}^{УЗ} * n_{уз}^{Треб} + К_{ИБП}^{МВЗ} * L_{МВЗ}^{Треб}) * 1,1= $' + " " + str(N_ibp_treb))

                    P_treb = round((tanks_count * k_y_pr * gm) / (round(accessible_km, 1) * Nz), 2)
                    st.markdown(
                        r'$П_{Треб} = (N_{Т}^{Пр} * K_{у}^{Пр} * g_{M}) / (L_{тд} * Н_{з}) = $' + " " + str(P_treb))

                    N_poter = round(P_treb * round(accessible_km, 1) * Nz, 2)
                    st.markdown(
                        r'$N_{Т}^{Потерь} = П_{Треб} * L_{тд} * Н_{з} = $' + " " + str(N_poter))


        elif not line_geom:
            st.info("ℹ️ Нарисуйте линию, чтобы проанализировать протяженность танкодоступной местности.")
    else:
        st.info("📝 Нарисуйте полигон для анализа дорог на участке местности.")
