import inspect
import json
import logging
from functools import partial
from logging import config
from pathlib import Path
from typing import Optional
from typing import Union

import fiona
import pyproj
import shapely.geometry as geo
import shapely.ops as ops
from custom_inherit import doc_inherit
from rasterio.crs import CRS

from miranda.scripting import LOGGING_CONFIG

config.dictConfig(LOGGING_CONFIG)
WGS84 = "+init=epsg:4326"
WGS84_PROJ4 = "+proj=longlat +datum=WGS84 +no_defs"

__all__ = [
    "feature_envelope",
    "layer_envelope",
    "feature_convex_hull",
    "layer_convex_hull",
    "feature_centroid",
    "layer_centroid",
    "feature_contained",
    "layer_contained",
]


def _geom_transform(geom, source_crs=WGS84_PROJ4, target_crs=None):
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : shapely.geometry
      Source geometry.
    source_crs : Union[str, CRS]
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    shapely.geometry
      Reprojected geometry.
    """
    try:
        reprojected = ops.transform(
            partial(pyproj.transform, pyproj.Proj(source_crs), pyproj.Proj(target_crs)),
            geom,
        )
        return reprojected
    except Exception as e:
        msg = "{}: Failed to reproject geometry".format(e)
        logging.error(msg)
        raise Exception(msg)


def _layer_operation(func):
    @doc_inherit(parent=func, style="numpy_napoleon")
    def wrapper(
        vector,
        output,
        source_crs: Optional[Union[str, CRS]] = None,
        target_crs: Optional[Union[str, CRS]] = None,
        **kwargs,
    ):
        """Parameters
        ----------
        vector : Union[str, Path]
          Path to a file containing a valid vector layer.
        output: Union[str, Path]
          Path to a file to be written.
        source_crs : Optional[Union[str, CRS]]
          Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
        target_crs : Optional[Union[str, CRS]]
          Projection identifier (proj4) for the target geometry.

        Returns
        -------
        None
        """
        if target_crs is None:
            msg = "No target CRS is defined. No vector transform will occur"
            logging.warning(msg)

        try:
            count = len(fiona.listlayers(vector))
        except TypeError:
            count = len(fiona.listlayers(str(vector)))

        for i in range(count):
            with fiona.open(vector, "r", layer=i) as src, open(output, "w") as sink:
                output_dict = {"type": "FeatureCollection", "features": []}
                for feature in src:
                    try:
                        geom = geo.shape(feature["geometry"])
                        # Perform vector reprojection using Shapely on each feature
                        if target_crs:
                            transformed = _geom_transform(geom, source_crs, target_crs)
                        else:
                            transformed = geom
                        feat_properties = feature["properties"]
                        func_geometry, func_properties = func(
                            vector=vector,
                            output=output,
                            source_crs=source_crs,
                            target_crs=target_crs,
                            geom=transformed,
                            prop=feat_properties,
                            **kwargs,
                        )
                        feature["geometry"] = geo.mapping(func_geometry)
                        feature["properties"] = func_properties
                        output_dict["features"].append(feature)
                    except Exception as e:
                        msg = "{}: Unable to process feature {}".format(e, feature)
                        logging.exception(msg)

                sink.write("{}".format(json.dumps(output_dict)))
        return

    sig = inspect.signature(wrapper)
    sig = sig.replace(parameters=tuple(sig.parameters.values())[:-1])
    wrapper.__signature__ = sig
    return wrapper


def _feature_operation(func):
    @doc_inherit(parent=func, style="numpy_napoleon")
    def wrapper(
        vector,
        output,
        source_crs: Optional[Union[str, CRS]] = None,
        target_crs: Optional[Union[str, CRS]] = None,
        **kwargs,
    ):
        """Parameters
        ----------
        vector : Union[str, Path]
          Path to a file containing a valid vector layer.
        output: Union[str, Path]
          Path to a folder to be written to.
        source_crs : Optional[Union[str, CRS]]
          Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
        target_crs : Optional[Union[str, CRS]]
          Projection identifier (proj4) for the target geometry.

        Returns
        -------
        None
        """
        if target_crs is None:
            msg = "No target CRS is defined. No vector transform will occur."
            logging.warning(msg)

        if not isinstance(output, Path):
            output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        try:
            count = len(fiona.listlayers(vector))
        except TypeError:
            count = len(fiona.listlayers(str(vector)))

        for i in range(count):
            with fiona.open(vector, "r", layer=i) as src:
                for num, feature in enumerate(src):
                    projected = output.joinpath(
                        "{}_{}.json".format(Path(vector).stem, num)
                    )

                    geom = geo.shape(feature["geometry"])
                    # Perform vector reprojection using Shapely on each feature
                    if target_crs:
                        transformed = _geom_transform(geom, source_crs, target_crs)
                    else:
                        transformed = geom
                    feat_properties = feature["properties"]
                    func_geometry, func_properties = func(
                        vector=vector,
                        output=output,
                        source_crs=source_crs,
                        target_crs=target_crs,
                        geom=transformed,
                        prop=feat_properties,
                        **kwargs,
                    )
                    feature["geometry"] = geo.mapping(func_geometry)
                    feature["properties"] = func_properties

                    with open(projected, "w") as sink:
                        outfile = {"type": "FeatureCollection", "features": []}
                        try:
                            outfile["features"].append(feature)
                        except Exception as e:
                            msg = "{}: Unable to process feature {}".format(
                                e, feature["id"]
                            )
                            logging.exception(msg)
                        sink.write("{}".format(json.dumps(outfile)))
        return kwargs

    sig = inspect.signature(wrapper)
    sig = sig.replace(parameters=tuple(sig.parameters.values())[:-1])
    wrapper.__signature__ = sig
    return wrapper


@_feature_operation
def feature_convex_hull(**kwargs):
    """Create convex hulls for all features within a single-layer vector file and return multiple GeoJSON files."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.convex_hull, prop


@_layer_operation
def layer_convex_hull(**kwargs):
    """Create convex hulls for all features within a single-layer vector file and a single GeoJSON file."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.convex_hull, prop


@_feature_operation
def feature_envelope(**kwargs):
    """Create envelopes for all features within a single-layer vector file and return multiple GeoJSON files."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.envelope, prop


@_layer_operation
def layer_envelope(**kwargs):
    """Create envelopes for all features within a single-layer vector file and return a single GeoJSON file."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.envelope, prop


@_feature_operation
def feature_centroid(**kwargs):
    """Calculate centroids for all features within a single-layer vector file and return multiple GeoJSON files."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.centroid, prop


@_layer_operation
def layer_centroid(**kwargs):
    """Calculate centroids for all features within a single-layer vector file and return a single GeoJSON file."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        return geom.centroid, prop


@_feature_operation
def feature_contained(boundary: geo, **kwargs):
    """Analyse containment for all features within a single-layer vector file according to a Geometry
    and return multiple GeoJSON files."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        prop["valid"] = boundary.contains(geom)
        return geom, prop


@_layer_operation
def layer_contained(boundary: geo, **kwargs):
    """Analyse containment for all features within a single-layer vector file according to a Geometry
    and return a single GeoJSON file."""
    geom, prop = kwargs["geom"], kwargs["prop"]
    if isinstance(geom, geo.Polygon):
        prop["valid"] = boundary.contains(geom)
        return geom, prop


if __name__ == "__main__":
    source = Path("/home/tjs/Downloads/GIS_Data/hydro_s")
    vec = source.joinpath("hydro_s.shp")
    output_folder = source.joinpath("output")
    output_file = source.joinpath("output.json")

    bound = geo.asPolygon(
        [(-71.5, 49), (-71.5, 54), (-63.5, 54), (-63.5, 49), (-71.5, 49)]
    )

    # feature_convex_hull(vector=vec, output=output_folder)
    feature_contained(boundary=bound, vector=vec, output=output_folder)
    print(feature_envelope.__doc__, "\n")
    print(inspect.signature(feature_envelope))
    layer_contained(boundary=bound, vector=vec, output=output_file)
    print(layer_envelope.__doc__, "\n")
    print(inspect.signature(layer_envelope))
