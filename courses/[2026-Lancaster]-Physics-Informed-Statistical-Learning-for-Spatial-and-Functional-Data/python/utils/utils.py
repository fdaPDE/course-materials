## Utilities functions

## Libraries
import numpy as np
import folium
import geopandas as gpd
from IPython.display import HTML, display
import math
import uuid
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import matplotlib
import branca.colormap as bcm
import matplotlib.colors as mcolors
from shapely.geometry import Polygon, Point
import rasterio
import rasterio.warp


# Function to plot mesh
def plot_mesh(
    triangulation,
    boundary_nodes=None,
    domain_shape=None,
    map_location=None,
    zoom_start=7,
    tiles="CartoDB positron",
    domain_layer_name="domain",
    domain_color="black",
    domain_weight=2,
    domain_fill=True,
    domain_fill_color="grey",
    domain_fill_opacity=0.4,
    domain_show=True,
    mesh_layer_name="mesh",
    mesh_color="black",
    mesh_weight=0.8,
    mesh_opacity=0.8,
    mesh_show=True,
    map_height="500px"   
):
    """
    Plot a spatial domain and a triangular mesh on an interactive Folium map.

    The domain and the mesh are shown as separate layers that can be toggled
    using the layer control.

    Exactly **one** of `boundary_nodes` or `domain_shape` must be provided.

    Parameters
    ----------
    triangulation : object
        Object with attributes:
        - nodes : array-like (N, 2) in (lon, lat)
        - cells : array-like (M, 3) with triangle vertex indices.
    boundary_nodes : array-like (K, 2), optional
        Domain boundary nodes as (lon, lat). Used if `domain_shape` is None.
    domain_shape : geopandas.GeoDataFrame, optional
        GeoDataFrame containing the domain geometry.
        Used if `boundary_nodes` is None.
    map_location : list or None, default=None
        Map center [lat, lon]. If None, the centroid of the mesh nodes is used.
    zoom_start : int, default=7
        Initial zoom level.
    tiles : str, default="CartoDB positron"
        Folium basemap.

    Domain style parameters
    -----------------------
    domain_layer_name : str
        Name of the domain layer.
    domain_color : str
        Domain border color.
    domain_weight : float
        Domain border width.
    domain_fill : bool
        Whether to fill the domain polygon (boundary_nodes case only).
    domain_fill_color : str
        Domain fill color.
    domain_fill_opacity : float
        Domain fill opacity.
    domain_show : bool
        Whether the domain layer is shown by default.

    Mesh style parameters
    --------------------
    mesh_layer_name : str
        Name of the mesh layer.
    mesh_color : str
        Mesh edge color.
    mesh_weight : float
        Mesh line width.
    mesh_opacity : float
        Mesh line opacity.
    mesh_show : bool
        Whether the mesh layer is shown by default.
    map_height: str
        Map height.

    Returns
    -------
    folium.Map
        Folium map object with domain and mesh layers.
    """

    # --- Check input
    if boundary_nodes is None and domain_shape is None:
        raise ValueError("Provide boundary_nodes or domain_shape.")
    if boundary_nodes is not None and domain_shape is not None:
        raise ValueError("Provide only one of boundary_nodes or domain_shape.")

    # --- Extract mesh
    nodes = np.asarray(triangulation.nodes)
    cells = np.asarray(triangulation.cells)

    if map_location is None:
        map_location = [nodes[:, 1].mean(), nodes[:, 0].mean()]

    # --- Create map
    m = folium.Map(
        location=map_location,
        zoom_start=zoom_start,
        tiles=tiles,
        width="100%",
        height=map_height
    )

    # --- Domain layer ---
    domain_fg = folium.FeatureGroup(name=domain_layer_name, show=domain_show)
    if boundary_nodes is not None:
        boundary_nodes = np.asarray(boundary_nodes)
        boundary_nodes = boundary_nodes[:, [1,0]]  # lat, lon
        folium.Polygon(
            locations=boundary_nodes,
            color=domain_color,
            weight=domain_weight,
            fill=domain_fill,
            fill_color=domain_fill_color,
            fill_opacity=domain_fill_opacity
        ).add_to(domain_fg)
    else:
        folium.GeoJson(
            domain_shape,
            style_function=lambda x: {
                "fillColor": domain_fill_color,
                "color": domain_color,
                "weight": domain_weight,
                "fillOpacity": domain_fill_opacity
            }
        ).add_to(domain_fg)
    domain_fg.add_to(m)

    # --- Mesh layer ---
    mesh_fg = folium.FeatureGroup(name=mesh_layer_name, show=mesh_show)
    for tri in cells:
        points = [[nodes[tri[i],1], nodes[tri[i],0]] for i in range(3)]
        points.append(points[0])
        folium.PolyLine(
            points, color=mesh_color, weight=mesh_weight, opacity=mesh_opacity
        ).add_to(mesh_fg)
    mesh_fg.add_to(m)

    folium.LayerControl(collapsed=False).add_to(m)

    # --- Display forcing height ---
    iframe_html = m._repr_html_()
    display(HTML(f"""
    <div style="width:100%; height:{map_height}; margin:0; padding:0;">
        <style>
            iframe {{
                width: 100% !important;
                height: {map_height} !important;
            }}
        </style>
        {iframe_html}
    </div>
    """))

# Function to plot areal data
def plot_areal_data(
    domain_shape,
    value_column,
    map_location=None,
    epsg_map=4326,
    epsg_centroid=32617,
    zoom_start=6,
    tiles="CartoDB positron",
    layer_name="Choropleth",
    cmap_name="viridis",     
    palette=None,            
    border_color="black",
    border_weight=1,
    fill_opacity=0.7,
    show_layer=True,
    show_layer_control=True,
    legend_name=None
):
    """
    Plot a choropleth map of a GeoDataFrame using Folium.
    """

    # Check inputs
    if value_column not in domain_shape.columns:
        raise ValueError(f"Column '{value_column}' not found in GeoDataFrame.")
    if not np.issubdtype(domain_shape[value_column].dtype, np.number):
        raise ValueError(f"Column '{value_column}' must be numeric.")

    domain = domain_shape.copy()

    if domain.crs is None:
        domain = domain.set_crs(epsg=epsg_map)
    domain = domain.to_crs(epsg=epsg_map)

    # Map center
    if map_location is None:
        domain_utm = domain.to_crs(epsg=epsg_centroid)
        centroid = domain_utm.geometry.centroid.to_crs(epsg=epsg_map)
        map_location = [centroid.y.mean(), centroid.x.mean()]

    m = folium.Map(
        location=map_location,
        zoom_start=zoom_start,
        tiles=tiles,
        width="100%",
        height="600px"
    )

    # Colormap / palette
    values = domain[value_column].values
    vmin = np.nanmin(values)
    vmax = np.nanmax(values)

    if palette is not None:
        # crea colormap continua dalla lista di colori
        colormap = bcm.LinearColormap(
            colors=palette,
            vmin=vmin,
            vmax=vmax
        )
    else:
        # usa colormap matplotlib
        cmap = matplotlib.colormaps[cmap_name]
        colormap = bcm.LinearColormap(
            colors=[matplotlib.colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, 256)],
            vmin=vmin,
            vmax=vmax
        )

    if legend_name is None:
        legend_name = value_column
    colormap.caption = legend_name

    # Feature group
    fg = folium.FeatureGroup(name=layer_name, show=show_layer)

    # Style function
    def style_function(feature):
        value = feature["properties"][value_column]
        return {
            "fillColor": colormap(value) if value is not None else "transparent",
            "color": border_color,
            "weight": border_weight,
            "fillOpacity": fill_opacity,
        }

    folium.GeoJson(
        domain,
        style_function=style_function,
        name=layer_name,
        tooltip=folium.GeoJsonTooltip(
            fields=[value_column],
            aliases=[f"{value_column}:"],
            localize=True
        )
    ).add_to(fg)

    fg.add_to(m)
    colormap.add_to(m)

    if show_layer_control:
        folium.LayerControl(collapsed=False).add_to(m)

    html = m.get_root()._repr_html_()
    HTML(f"""
    <div style="
        width: 100%;
        height: 600px;
        margin: 0;
        padding: 0;
        overflow: hidden;
    ">
        {html}
    </div>
    """)

    return m

# Function to plot the domain
def plot_domain(
    boundary_nodes=None,
    domain_shape=None,
    map_location=None,
    epsg_map=4326,
    epsg_centroid=32617,
    zoom_start=6,
    tiles="CartoDB positron",
    layer_name="Domain",
    border_color="black",
    border_weight=2,
    fill=True,
    fill_color="gray",
    fill_opacity=0.25,
    show_layer=True,
    show_layer_control=True,
    popup_label = None,
    map_height="500px" 
):
    """
    Plot a spatial domain on an interactive Folium map.

    The function supports two input types:
    1. `domain_shape` : a GeoDataFrame containing the domain geometry.
    2. `boundary_nodes` : an array-like of boundary nodes (lon, lat).

    Exactly **one** of `boundary_nodes` or `domain_shape` must be provided.

    Parameters
    ----------
    boundary_nodes : array-like (N, 2), optional
        Boundary nodes as (lon, lat). Only used if `domain_shape` is None.
    domain_shape : geopandas.GeoDataFrame, optional
        GeoDataFrame containing the domain geometry. Only used if `boundary_nodes` is None.
    map_location : list or None, default=None
        Map center [lat, lon]. If None, uses centroid of boundary_nodes. Only used if `boundary_nodes` is provided.
    epsg_map : int, default=4326
        CRS used for map visualization (WGS84). Only used if `domain_shape` is provided.
    epsg_centroid : int, default=32617
        Metric CRS used to compute centroid safely. Only used if `domain_shape` is provided.
    zoom_start : int, default=6
        Initial zoom level.
    tiles : str, default="CartoDB positron"
        Folium basemap.
    layer_name : str, default="Domain"
        Name of the layer.
    border_color : str, default="black"
        Border color of the polygon(s).
    border_weight : float, default=2
        Width of polygon border.
    fill : bool, default=True
        Whether to fill the polygon. Only used for boundary_nodes.
    fill_color : str, default="gray"
        Fill color of the polygon(s).
    fill_opacity : float, default=0.25
        Polygon fill opacity.
    show_layer : bool, default=True
        Whether the layer is visible by default.
    show_layer_control : bool, default=True
        Whether to add a layer control.
    popup_label : str or None, optional
        Name of the GeoDataFrame column to be shown in a popup.
        If None, no popup is added. Only used if `domain_shape` is provided.
    map_height: str
        Map height

    Returns
    -------
    folium.Map
        Folium map object with the domain plotted.
    """

    # Check input
    if boundary_nodes is None and domain_shape is None:
        raise ValueError("You must provide either boundary_nodes or domain_shape.")

    if boundary_nodes is not None and domain_shape is not None:
        raise ValueError("Provide only one of boundary_nodes or domain_shape, not both.")

    # Case 1: boundary_nodes provided
    if boundary_nodes is not None:
        boundary_nodes = np.asarray(boundary_nodes)

        if boundary_nodes.shape[1] != 2:
            raise ValueError("Boundary nodes must have shape (N, 2) as (lon, lat).")

        # Convert to (lat, lon) for Folium
        boundary_nodes = boundary_nodes[:, [1, 0]]

        # Map center
        if map_location is None:
            map_location = [
                boundary_nodes[:, 0].mean(),
                boundary_nodes[:, 1].mean()
            ]

        # Create map
        m = folium.Map(location=map_location, zoom_start=zoom_start, tiles=tiles, width="100%", height="600px" ) 


        # Feature group
        domain_fg = folium.FeatureGroup(name=layer_name, show=show_layer)

        # Add polygon
        folium.Polygon(
            locations=boundary_nodes,
            color=border_color,
            weight=border_weight,
            fill=fill,
            fill_color=fill_color,
            fill_opacity=fill_opacity
        ).add_to(domain_fg)

        domain_fg.add_to(m)

    # Case 2: domain_shape provided
    else:
        domain = domain_shape.copy()

        # Set CRS if missing
        if domain.crs is None:
            domain = domain.set_crs(epsg=epsg_map)

        # Reproject for visualization
        domain = domain.to_crs(epsg=epsg_map)

        # Compute centroid in metric CRS
        domain_utm = domain.to_crs(epsg=epsg_centroid)
        centroid = domain_utm.geometry.centroid.to_crs(epsg=epsg_map)

        # Map center
        center = [centroid.y.mean(), centroid.x.mean()]

        # Create map
        m = folium.Map(location=center, zoom_start=zoom_start, tiles=tiles)

        # Feature group
        domain_fg = folium.FeatureGroup(name=layer_name, show=show_layer)

        # Add GeoDataFrame polygon(s)
        # Prepare popup (only if requested)
        popup = None
        if popup_label is not None:
            if popup_label not in domain.columns:
                raise ValueError(
                    f"Column '{popup_label}' not found in domain_shape GeoDataFrame."
                )
            popup = folium.GeoJsonPopup(
                fields=[popup_label],
                aliases=[f"{popup_label}:"],
                localize=True
            )

        # Add GeoDataFrame polygon(s)
        folium.GeoJson(
            domain_shape,
            name=layer_name,
            style_function=lambda x: {
                "fillColor": fill_color,
                "color": border_color,
                "weight": border_weight,
                "fillOpacity": fill_opacity
            },
            popup=popup
        ).add_to(domain_fg)

        domain_fg.add_to(m)

    # Layer control
    if show_layer_control:
        folium.LayerControl(collapsed=False).add_to(m)

    # --- Display forcing height ---
    iframe_html = m._repr_html_()
    display(HTML(f"""
    <div style="width:100%; height:{map_height}; margin:0; padding:0;">
        <style>
            iframe {{
                width: 100% !important;
                height: {map_height} !important;
            }}
        </style>
        {iframe_html}
    </div>
    """))

# Function to plot locations
def plot_locations(
    locations,
    boundary_nodes=None,
    domain_shape=None,
    zoom_start=7,
    tiles="cartodb positron",

    # --- marker style ---
    marker_radius=6,
    marker_weight=1,
    marker_fill_opacity=0.6,
    marker_fill_color="red",
    marker_border_color="black",

    # --- domain style ---
    domain_color="black",
    domain_weight=2,
    domain_fill=True,
    domain_fill_color="grey",
    domain_fill_opacity=0.5,

    # --- CRS options (used if domain_shape provided) ---
    epsg_map=4326,
    epsg_centroid=32617 
):
    """
    Plot spatial point locations on a Folium map.

    The function supports two domain definitions:
    1. boundary_nodes : array-like (N, 2) in (lon, lat)
    2. domain_shape   : GeoDataFrame containing polygon geometry

    Exactly one of boundary_nodes or domain_shape must be provided.

    Parameters
    ----------
    locations : array-like (M, 2)
        Point coordinates given as (lon, lat).
    boundary_nodes : array-like (N, 2), optional
        Domain boundary nodes in (lon, lat).
    domain_shape : GeoDataFrame, optional
        GeoDataFrame describing the spatial domain.
    zoom_start : int
        Initial zoom level.
    tiles : str
        Folium basemap.
    marker_radius : float
        Circle marker radius.
    marker_weight : float
        Marker border width.
    marker_fill_opacity : float
        Marker fill opacity.
    marker_fill_color : str
        Marker fill color.
    marker_border_color : str
        Marker border color.
    domain_color : str
        Domain border color.
    domain_weight : float
        Domain border width.
    domain_fill : bool
        Whether to fill the domain polygon.
    domain_fill_color : str
        Domain fill color.
    domain_fill_opacity : float
        Domain fill opacity.
    epsg_map : int
        CRS used for visualization (default WGS84).
    epsg_centroid : int
        Metric CRS used to compute centroid safely.

    Returns
    -------
    folium.Map
        Interactive map with domain and locations.
    """


    # Input checks
    if boundary_nodes is None and domain_shape is None:
        raise ValueError("Provide either boundary_nodes or domain_shape.")

    if boundary_nodes is not None and domain_shape is not None:
        raise ValueError("Provide only one of boundary_nodes or domain_shape.")

    coords = np.asarray(locations)

    if coords.shape[1] != 2:
        raise ValueError("locations must have shape (M, 2) as (lon, lat).")

    # --------------------------------------------------
    # Map center
    if domain_shape is not None:

        domain = domain_shape.copy()

        if domain.crs is None:
            domain = domain.set_crs(epsg=epsg_map)

        domain = domain.to_crs(epsg=epsg_map)

        domain_utm = domain.to_crs(epsg=epsg_centroid)
        centroid = domain_utm.geometry.centroid.to_crs(epsg=epsg_map)

        map_center = [centroid.y.mean(), centroid.x.mean()]

    else:
        map_center = [coords[:, 1].mean(), coords[:, 0].mean()]

    # --------------------------------------------------
    # Create map
    m = folium.Map(
        location=map_center,
        zoom_start=zoom_start,
        tiles=tiles
    )

    # --------------------------------------------------
    # DOMAIN
    domain_fg = folium.FeatureGroup(name="domain", show=True)

    if boundary_nodes is not None:

        boundary_nodes = np.asarray(boundary_nodes)

        if boundary_nodes.shape[1] != 2:
            raise ValueError("boundary_nodes must have shape (N, 2).")

        folium.Polygon(
            locations=boundary_nodes[:, [1, 0]],  # lat, lon
            color=domain_color,
            weight=domain_weight,
            fill=domain_fill,
            fill_color=domain_fill_color,
            fill_opacity=domain_fill_opacity
        ).add_to(domain_fg)

    else:
        folium.GeoJson(
            domain,
            style_function=lambda x: {
                "fillColor": domain_fill_color if domain_fill else "none",
                "color": domain_color,
                "weight": domain_weight,
                "fillOpacity": domain_fill_opacity if domain_fill else 0
            }
        ).add_to(domain_fg)

    domain_fg.add_to(m)

    # --------------------------------------------------
    # LOCATIONS
    data_fg = folium.FeatureGroup(name="locations", show=True)

    for i in range(coords.shape[0]):
        folium.CircleMarker(
            location=[coords[i, 1], coords[i, 0]],  # lat, lon
            radius=marker_radius,
            color=marker_border_color,
            weight=marker_weight,
            fill=True,
            fill_color=marker_fill_color,
            fill_opacity=marker_fill_opacity
        ).add_to(data_fg)

    data_fg.add_to(m)

    # --------------------------------------------------
    folium.LayerControl(collapsed=False).add_to(m)

    html = m.get_root()._repr_html_()
    HTML(f"""
    <div style="
        width: 100%;
        height: 600px;
        margin: 0;
        padding: 0;
        overflow: hidden;
    ">
        {html}
    </div>
    """)

    return m

# Function to plot areal estimate 
def plot_areal_estimate(
    domain_shape,
    values,
    map_location=None,
    epsg_map=4326,
    epsg_centroid=32617,
    zoom_start=6,
    tiles="CartoDB positron",
    layer_name="Choropleth",
    cmap_name="viridis",
    palette=None,
    border_color="black",
    border_weight=1,
    fill_opacity=0.7,
    show_layer=True,
    show_layer_control=True,
    legend_name=None
):
    """
    Plot a choropleth map of a GeoDataFrame using a numeric NumPy array.

    Parameters
    ----------
    domain_shape : geopandas.GeoDataFrame
        GeoDataFrame containing polygon geometries.
    values : numpy.ndarray
        1D numeric array with length equal to len(domain_shape).
        Each element corresponds to the polygon at the same index.
    map_location : list [lat, lon], optional
        Initial centre of the map. If None, computed from centroids.
    epsg_map : int, default=4326
        EPSG code for map projection (typically WGS84).
    epsg_centroid : int, default=32617
        EPSG code used to compute centroids more accurately.
    zoom_start : int, default=6
        Initial zoom level.
    tiles : str, default="CartoDB positron"
        Folium tile style.
    layer_name : str, default="Choropleth"
        Name of the map layer.
    cmap_name : str, default="viridis"
        Matplotlib colormap name (used if palette is None).
    palette : list of str, optional
        List of hex colours defining a custom continuous palette.
    border_color : str, default="black"
        Polygon border colour.
    border_weight : int, default=1
        Polygon border thickness.
    fill_opacity : float, default=0.7
        Fill transparency.
    show_layer : bool, default=True
        Whether the layer is visible at load.
    show_layer_control : bool, default=True
        Whether to show layer control.
    legend_name : str, optional
        Legend title. If None, defaults to "ESTIMATE".

    Returns
    -------
    folium.Map
        Interactive Folium map object.
    """

    # --- Input checks ---
    if not isinstance(values, np.ndarray):
        raise TypeError("`values` must be a numpy.ndarray.")

    if values.ndim != 1:
        raise ValueError("`values` must be a 1D numpy array.")

    if len(values) != len(domain_shape):
        raise ValueError(
            "`values` length must match number of geometries in domain_shape."
        )

    if not np.issubdtype(values.dtype, np.number):
        raise ValueError("`values` must contain numeric data.")

    domain = domain_shape.copy()
    domain["_values_"] = values

    # --- CRS handling ---
    if domain.crs is None:
        domain = domain.set_crs(epsg=epsg_map)
    domain = domain.to_crs(epsg=epsg_map)

    # --- Map centre ---
    if map_location is None:
        domain_utm = domain.to_crs(epsg=epsg_centroid)
        centroid = domain_utm.geometry.centroid.to_crs(epsg=epsg_map)
        map_location = [centroid.y.mean(), centroid.x.mean()]

    m = folium.Map(
        location=map_location,
        zoom_start=zoom_start,
        tiles=tiles,
        width="100%",
        height="600px"
    )

    # --- Colormap ---
    vmin = np.nanmin(values)
    vmax = np.nanmax(values)

    if palette is not None:
        colormap = bcm.LinearColormap(
            colors=palette,
            vmin=vmin,
            vmax=vmax
        )
    else:
        cmap = matplotlib.colormaps[cmap_name]
        colormap = bcm.LinearColormap(
            colors=[
                matplotlib.colors.rgb2hex(cmap(i))
                for i in np.linspace(0, 1, 256)
            ],
            vmin=vmin,
            vmax=vmax
        )

    if legend_name is None:
        legend_name = "ESTIMATE"

    colormap.caption = legend_name

    fg = folium.FeatureGroup(name=layer_name, show=show_layer)

    def style_function(feature):
        value = feature["properties"]["_values_"]
        return {
            "fillColor": colormap(value) if value is not None else "transparent",
            "color": border_color,
            "weight": border_weight,
            "fillOpacity": fill_opacity,
        }

    folium.GeoJson(
        domain,
        style_function=style_function,
        name=layer_name,
        tooltip=folium.GeoJsonTooltip(
            fields=["_values_"],
            aliases=[f"{legend_name}:"],
            localize=True
        )
    ).add_to(fg)

    fg.add_to(m)
    colormap.add_to(m)

    if show_layer_control:
        folium.LayerControl(collapsed=False).add_to(m)

    html = m.get_root()._repr_html_()
    HTML(f"""
    <div style="
        width: 100%;
        height: 600px;
        margin: 0;
        padding: 0;
        overflow: hidden;
    ">
        {html}
    </div>
    """)

    return m

# Function to plot satellite data on folium map
def display_satellite_data(
    raster_path,
    layer_name,
    boundary_nodes=None,
    cmap_name="magma",
    palette=None,
    opacity=0.8,
    zoom_start=6,
    tiles="cartodb positron",
    map_location=None,
    show_layer=True,
    border_color="black",
    border_weight=2,
    fill=True,
    fill_color="gray",
    fill_opacity=0.00,
):


    with rasterio.open(raster_path) as src:
        data = src.read(1, masked=True) - 273.15  
        bounds = src.bounds
        crs = src.crs

    vmin = np.nanmin(data)
    vmax = np.nanmax(data)

    if palette is not None:
        colormap = bc.LinearColormap(
            colors=palette, vmin=vmin, vmax=vmax, caption=layer_name
        )
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        rgba_img = np.zeros((*data.shape, 4))

        vals = data.filled(np.nan)
        mask = np.isnan(vals)

        rgba_img[~mask] = [
            mcolors.to_rgba(colormap(v))
            for v in vals[~mask]
        ]
        rgba_img[mask] = [0, 0, 0, 0]

        rgba_img = (rgba_img * 255).astype(np.uint8)

    else:
        cmap = matplotlib.colormaps[cmap_name]
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        rgba_img = cmap(norm(data))
        rgba_img = (rgba_img * 255).astype(np.uint8)

        colormap = bc.LinearColormap(
            colors=[cmap(i) for i in np.linspace(0, 1, 256)],
            vmin=vmin,
            vmax=vmax,
            caption=layer_name,
        )
        

    # Map center
    if map_location is None:
        center_lat = (bounds.bottom + bounds.top) / 2
        center_lon = (bounds.left + bounds.right) / 2
        map_location = [center_lat, center_lon]
            
    m = folium.Map(
        location= map_location,
        zoom_start=zoom_start,
        tiles=tiles,
    )

    if boundary_nodes is not None:
        boundary_nodes = np.asarray(boundary_nodes)

        if boundary_nodes.shape[1] != 2:
            raise ValueError("Boundary nodes must have shape (N, 2) as (lon, lat).")

        # Convert to (lat, lon) for Folium
        boundary_nodes = boundary_nodes[:, [1, 0]]

        # Feature group
        domain_fg = folium.FeatureGroup(name=layer_name, show=show_layer)

        # Add polygon
        folium.Polygon(
            locations=boundary_nodes,
            color=border_color,
            weight=border_weight,
            fill=fill,
            fill_color=fill_color,
            fill_opacity=fill_opacity
        ).add_to(domain_fg)

        domain_fg.add_to(m)

    folium.raster_layers.ImageOverlay(
        image=rgba_img,
        bounds=[
            [bounds.bottom, bounds.left],
            [bounds.top, bounds.right],
        ],
        opacity=opacity,
        name=layer_name,
        interactive=True,
    ).add_to(m)

    colormap.add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)

    html = m.get_root()._repr_html_()
    HTML(f"""
    <div style="
        width: 100%;
        height: 600px;
        margin: 0;
        padding: 0;
        overflow: hidden;
    ">
    {html}
    </div>
    """)

    return m

# Function to plot maps according to a user defined grid imposing the cell height
def display_folium_map(
    maps,
    nrows = 1,
    ncols = 1,
    cell_height="550px",  
    gap="10px"
):
    if len(maps) > nrows * ncols:
        raise ValueError(
            f"Too many maps ({len(maps)}) for a {nrows}x{ncols} grid."
        )

    map_divs = []
    for m in maps:
        map_divs.append(
            f"""
            <div style="
                width: 100%;
                height: {cell_height};
                position: relative;
            ">
                <style>
                    iframe {{
                        width: 100% !important;
                        height: {cell_height} !important;
                    }}
                </style>
                {m.get_root()._repr_html_()}
            </div>
            """
        )

    # celle vuote
    for _ in range(nrows * ncols - len(maps)):
        map_divs.append(
            f'<div style="width:100%; height:{cell_height};"></div>'
        )

    grid_html = f"""
    <div style="
        display: grid;
        grid-template-columns: repeat({ncols}, 1fr);
        gap: {gap};
        width: 100%;
    ">
        {''.join(map_divs)}
    </div>
    """

    display(HTML(grid_html))

# Auxiliary function
def _mpl_fig_to_html(fig, dpi=150):
    """Convert a matplotlib figure to an HTML <img> tag with correct aspect ratio."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")

    # Larghezza al 100% della cella, altezza auto, centrata
    html = f"""
    <div style="display:flex; justify-content:center; align-items:center; width:100%; height:100%;">
        <img src="data:image/png;base64,{img_base64}"
             style="max-width:100%; max-height:100%; height:auto; width:auto;">
    </div>
    """
    return html

