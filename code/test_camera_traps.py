#!/usr/bin/env python3
"""
test_camera_traps.py
====================
Unit tests for the helper and logic functions in code/process_camera_traps.py.
Runs quickly without requiring any real camera trap datasets.
"""

import sys
import unittest
import numpy as np
import pandas as pd
from pathlib import Path

# Add the script's directory to python path to allow importing process_camera_traps
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR))

try:
    import process_camera_traps
except ImportError as e:
    print(f"Error importing process_camera_traps: {e}")
    raise


class TestCameraTraps(unittest.TestCase):
    def test_safe_str(self):
        self.assertEqual(process_camera_traps.safe_str(" Hello "), "Hello")
        self.assertEqual(process_camera_traps.safe_str(np.nan), "")
        self.assertEqual(process_camera_traps.safe_str(None), "")
        self.assertEqual(process_camera_traps.safe_str(123), "123")

    def test_assign_taxon_quality(self):
        # 1. Human detections
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"common_name": "Human"})), "human")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Homo", "species": "sapiens"})), "human")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"common_name": "Trapper interaction"})), "human")

        # 2. Blanks / Setup / Unknown
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"common_name": "Blank"})), "blank")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Setup"})), "blank")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"common_name": ""})), "blank")

        # 3. Domestic species
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Bos", "species": "taurus", "common_name": "Domestic Cattle"})), "domestic")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Canis", "species": "familiaris", "common_name": "Domestic Dog"})), "domestic")

        # 4. Valid species / genus / family
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Panthera", "species": "onca", "common_name": "Jaguar"})), "species")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Panthera", "species": "sp.", "common_name": "Panthera species"})), "genus")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "Panthera", "species": "", "common_name": "Panthera species"})), "genus")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "", "family": "Felidae", "common_name": "Cats"})), "family")
        self.assertEqual(process_camera_traps.assign_taxon_quality(pd.Series({"genus": "", "family": "", "common_name": "Unknown mammal"})), "higher")

    def test_haversine(self):
        # Known distances:
        # Distance between (0, 0) and (0, 1) latitude degree is approx 111 km
        d1 = process_camera_traps.haversine(0.0, 0.0, 0.0, 1.0)
        self.assertAlmostEqual(d1, 111.19, places=1)

        # Distance between (0, 0) and (1, 0) longitude degree on equator is approx 111 km
        d2 = process_camera_traps.haversine(0.0, 0.0, 1.0, 0.0)
        self.assertAlmostEqual(d2, 111.19, places=1)

        # Identical points should have 0 distance
        self.assertEqual(process_camera_traps.haversine(12.34, -56.78, 12.34, -56.78), 0.0)

    def test_parse_wi_timestamp(self):
        ts = pd.Series([
            "2023-05-15 14:30:00 (GMT+1)",
            "2023-05-15T14:30:00Z",
            None
        ])
        parsed = process_camera_traps._parse_wi_timestamp(ts)
        self.assertEqual(parsed.iloc[0].year, 2023)
        self.assertEqual(parsed.iloc[1].year, 2023)
        self.assertTrue(pd.isna(parsed.iloc[2]))

    def test_collapse_to_independent_events(self):
        # Create image records for 2 different deployments/taxa
        # Deployment 1: 3 images of species A within 15 minutes, 1 image 45 minutes later
        # Deployment 2: 1 image of species B
        df = pd.DataFrame([
            {"project_id": "P1", "deployment_id": "D1", "genus": "Panthera", "species": "onca", "common_name": "Jaguar", "number_of_objects": 1, "timestamp": "2023-01-01 12:00:00"},
            {"project_id": "P1", "deployment_id": "D1", "genus": "Panthera", "species": "onca", "common_name": "Jaguar", "number_of_objects": 2, "timestamp": "2023-01-01 12:10:00"},
            {"project_id": "P1", "deployment_id": "D1", "genus": "Panthera", "species": "onca", "common_name": "Jaguar", "number_of_objects": 1, "timestamp": "2023-01-01 12:15:00"},
            {"project_id": "P1", "deployment_id": "D1", "genus": "Panthera", "species": "onca", "common_name": "Jaguar", "number_of_objects": 1, "timestamp": "2023-01-01 13:00:00"},
            {"project_id": "P1", "deployment_id": "D2", "genus": "Tapirus", "species": "terrestris", "common_name": "Tapir", "number_of_objects": 1, "timestamp": "2023-01-01 12:00:00"},
        ])

        collapsed = process_camera_traps._collapse_to_independent_events(df, threshold_minutes=30.0)
        
        # Expected:
        # Deployment D1 Panthera onca:
        #   - Event 1: max objects = 2, timestamp = 12:00:00
        #   - Event 2: max objects = 1, timestamp = 13:00:00
        # Deployment D2 Tapirus terrestris:
        #   - Event 1: max objects = 1, timestamp = 12:00:00
        # Total independent events = 3
        self.assertEqual(len(collapsed), 3)
        
        d1_events = collapsed[collapsed["deployment_id"] == "D1"]
        self.assertEqual(len(d1_events), 2)
        self.assertEqual(d1_events["number_of_objects"].max(), 2)

    def test_match_body_mass(self):
        # Mock traits database (species median and genus median)
        traits = pd.DataFrame([
            {"genus": "Panthera", "species": "onca", "body_mass_kg": 80.0},
            {"genus": "Panthera", "species": "pardus", "body_mass_kg": 60.0},
            {"genus": "Tapirus", "species": "terrestris", "body_mass_kg": 200.0},
        ])

        # Test detections
        det = pd.DataFrame([
            # 1. Exact species match
            {"genus": "Panthera", "species": "onca", "family": "Felidae"},
            # 2. Genus fallback (species isn't in traits, matches genus median = (80+60)/2 = 70.0)
            {"genus": "Panthera", "species": "leo", "family": "Felidae"},
            # 3. Family fallback (not in traits, matches FALLBACK_FAMILY_MASS_KG["Viverridae"] = 2.5)
            {"genus": "Civettictis", "species": "civetta", "family": "Viverridae"},
            # 4. Unmatched
            {"genus": "Unknown", "species": "unknown", "family": "UnknownFamily"}
        ])

        matched = process_camera_traps.match_body_mass(det, traits)
        self.assertEqual(matched.loc[0, "body_mass_kg"], 80.0)
        self.assertEqual(matched.loc[0, "mass_match_level"], "species")

        self.assertEqual(matched.loc[1, "body_mass_kg"], 70.0)
        self.assertEqual(matched.loc[1, "mass_match_level"], "genus")

        self.assertEqual(matched.loc[2, "body_mass_kg"], 2.5)
        self.assertEqual(matched.loc[2, "mass_match_level"], "family_fallback")

        self.assertTrue(pd.isna(matched.loc[3, "body_mass_kg"]))
        self.assertEqual(matched.loc[3, "mass_match_level"], "unmatched")

    def test_compute_site_heterotroph_metrics(self):
        # Create sample detection rate DataFrame
        det = pd.DataFrame([
            {
                "project_id": "P1", "deployment_id": "D1", "longitude": -60.0, "latitude": -3.0,
                "trap_days": 100.0, "taxon_key": "Panthera_onca", "class": "Mammalia", "order": "Carnivora",
                "family": "Felidae", "genus": "Panthera", "species": "onca", "taxon_quality": "species",
                "n_detections": 10.0, "RAI": 10.0, "body_mass_kg": 100.0, "common_name": "Jaguar"
            },
            {
                "project_id": "P1", "deployment_id": "D1", "longitude": -60.0, "latitude": -3.0,
                "trap_days": 100.0, "taxon_key": "Blank_None", "class": "", "order": "",
                "family": "", "genus": "", "species": "", "taxon_quality": "blank",
                "n_detections": 0.0, "RAI": 0.0, "body_mass_kg": np.nan, "common_name": ""
            }
        ])

        metrics = process_camera_traps.compute_site_heterotroph_metrics(det)
        self.assertEqual(len(metrics), 1)
        self.assertEqual(metrics.loc[0, "n_species"], 1)
        self.assertEqual(metrics.loc[0, "n_detections_total"], 10)
        self.assertTrue(metrics.loc[0, "B_H_index"] > 0)
        self.assertTrue(metrics.loc[0, "M_H_index"] > 0)
        # Jaguar is > 50 and > 100 kg, so B_H_gt50 and B_H_gt100 should be equal to B_H_index, fraction = 100
        self.assertAlmostEqual(metrics.loc[0, "megafauna_fraction"], 100.0)

    def test_calculate_temporal_weights_py(self):
        # D1 is set 1 year before GEDI start (2019-04-17) -> start_date = 2018-04-17
        # D2 is set 8 years before GEDI start -> start_date = 2011-04-17
        det = pd.DataFrame([
            {
                "region": "Amazon", "longitude": -60.0, "latitude": -3.0,
                "project_id": "P1", "deployment_id": "D1", "trap_days": 100.0,
                "start_date": "2018-04-17", "end_date": "2018-07-25"
            },
            {
                "region": "Amazon", "longitude": -61.0, "latitude": -3.0,
                "project_id": "P1", "deployment_id": "D2", "trap_days": 100.0,
                "start_date": "2011-04-17", "end_date": "2011-07-25"
            }
        ])
        cluster_map = {
            ("Amazon", -60.0, -3.0): "C1",
            ("Amazon", -61.0, -3.0): "C2"
        }

        weights = process_camera_traps.calculate_temporal_weights_py(det, cluster_map)
        # Years before GEDI:
        # D1: ~1.0 yr -> bracket[0] = 1.0 -> weight = 1.0
        # D2: ~8.0 yr -> between 6.0 and 11.0 -> weight = 0.25
        self.assertAlmostEqual(weights["C1"], 1.0, places=2)
        self.assertAlmostEqual(weights["C2"], 0.25, places=2)

    def test_calculate_taxonomic_keep_proportions_py(self):
        det = pd.DataFrame([
            # Wild species (kept)
            {"region": "Congo", "longitude": 15.0, "latitude": -1.0, "taxon_quality": "species", "n_detections": 8},
            # Higher resolution (not in species/genus/family -> excluded/dropped wild)
            {"region": "Congo", "longitude": 15.0, "latitude": -1.0, "taxon_quality": "higher", "n_detections": 2},
            # Non-wild (domestic/blank/human -> completely ignored in denom & numer)
            {"region": "Congo", "longitude": 15.0, "latitude": -1.0, "taxon_quality": "human", "n_detections": 5},
        ])
        cluster_map = {
            ("Congo", 15.0, -1.0): "C1"
        }
        p_keep = process_camera_traps.calculate_taxonomic_keep_proportions_py(det, cluster_map)
        # Total wild = 8 (species) + 2 (higher) = 10
        # Kept wild = 8 (species)
        # Expected: 8 / 10 = 0.8
        self.assertAlmostEqual(p_keep["C1"], 0.8)


if __name__ == "__main__":
    unittest.main()
