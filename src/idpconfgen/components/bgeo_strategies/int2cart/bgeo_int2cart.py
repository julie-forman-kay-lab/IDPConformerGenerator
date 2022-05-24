"""Using the machine learning Int2Cart algorithm to provide bond geometries."""
import os
from pathlib import Path

import torch
import yaml
from modelling.models.builder import BackboneBuilder
from modelling.utils.predict import predict


_folder = Path(__file__).parent


class BGEO_Int2Cart:
    """Prepara Int2Cart module."""

    def __init__(self) -> None:
        model_config = os.path.join(_folder, "int2cart.yml")
        model_addr = os.path.join(_folder, "model.tar")

        with open(model_config, "r") as fin:
            settings = yaml.safe_load(fin)

        builder = BackboneBuilder(settings)

        _ = torch.load(model_addr, map_location=torch.device('cuda'))
        model_state = _['model_state_dict']

        builder.load_predictor_weights(model_state)
        self.builder = builder

    def get_internal_coords(self, sequence, torsions):
        """Get internal coords."""
        predictions = predict(
            self.builder,
            sequence,
            torsions,
            build=False,
            nits="radian",
            )
        d1 = predictions['d1'][0, -1, 0]
        d2 = predictions['d2'][0, -1, 0]
        d3 = predictions['d3'][0, -1, 0]
        theta1 = predictions['theta1'][0, -1, 0]
        theta2 = predictions['theta2'][0, -1, 0]
        theta3 = predictions['theta3'][0, -1, 0]
        return (d1, d2, d3, theta1, theta2, theta3)
