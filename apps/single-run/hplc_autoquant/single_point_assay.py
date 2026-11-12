"""Single-point external-standard assay for one selected UV peak."""
from __future__ import annotations

import math


def calculate_single_point(standard_area, standard_concentration, unknown_area,
                           unknown_dilution_factor=1.0, response_factor=1.0):
    """Return unknown concentration from one standard and one unknown peak.

    The result is in the same concentration units entered for the standard.
    ``unknown_dilution_factor`` converts the measured unknown solution back to
    the requested original-sample basis.  A response factor is optional and
    defaults to one because no relative response calibration is assumed.
    """
    values = {
        'standard_area': standard_area,
        'standard_concentration': standard_concentration,
        'unknown_area': unknown_area,
        'unknown_dilution_factor': unknown_dilution_factor,
        'response_factor': response_factor,
    }
    if any(not isinstance(v, (int, float)) or not math.isfinite(float(v)) for v in values.values()):
        raise ValueError('Single-point assay inputs must be finite numbers')
    if standard_area <= 0 or standard_concentration <= 0 or unknown_area < 0:
        raise ValueError('Areas must be positive (unknown may be zero) and standard concentration must be positive')
    if unknown_dilution_factor <= 0 or response_factor <= 0:
        raise ValueError('Dilution and response factors must be positive')
    measured = (unknown_area / standard_area) * standard_concentration / response_factor
    return {
        'standard_area': float(standard_area),
        'standard_concentration': float(standard_concentration),
        'unknown_area': float(unknown_area),
        'unknown_dilution_factor': float(unknown_dilution_factor),
        'response_factor': float(response_factor),
        'unknown_concentration': float(measured * unknown_dilution_factor),
        'formula': '(unknown area / standard area) × standard concentration × unknown dilution factor ÷ response factor',
        'method': 'Single-point external standard',
        'provisional': True,
        'notes': [
            'Uses one standard injection and one unknown injection.',
            'Assumes equal detector response between standard and unknown unless a response factor is entered.',
            'This is not a calibration curve and does not establish linearity, accuracy, precision, identity or purity.',
        ],
    }
