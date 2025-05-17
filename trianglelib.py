import numpy as np
import warnings
from typing import List, Dict


# Represents a mass spectrum with associated isolation parameters and MS data
class Spectrum:
    def __init__(self, isolation_center: float, isolation_width: float):
        self.isolation_center: float = isolation_center
        self.isolation_width: float = isolation_width
        self.scan_number: int
        self.ms_data: np.nparray = np.array([[], []])

    def add_ms_data(self, mass_list: List[float], intensity_list: List[float]):
        self.ms_data = np.append(self.ms_data, [mass_list, intensity_list], axis=1)

    def get_ms_data(self):
        return self.ms_data

    def filter_data(self, low_intensity_cuttoff: float = 750):
        mass_list, intensity_list = self.ms_data
        data_filter = intensity_list >= low_intensity_cuttoff
        data_filter &= mass_list > 0
        filtered_mass_list, filtered_intensity_list = [mass_list[data_filter], intensity_list[data_filter]]
        self.ms_data = np.array([[], []])
        self.add_ms_data(filtered_mass_list, filtered_intensity_list)

# Models a triangular intensity profile around m/z found in consecutive scans
class TriangleProfile:
    def __init__(self, observed_mass: float, isolation_width: float):
        self.observed_mass: float = observed_mass
        self.isolation_width: float = isolation_width
        self.profile_data: np.nparray = np.array([[], []])

    def add_pofile_data(self, isolation_center_list: List[float], intensity_list: List[float]):
        self.profile_data = np.append(self.profile_data, [isolation_center_list, intensity_list], axis=1)

    def internal_calibrant_regression(self, knownMass):
        isolation_center_list, intensity_list = self.profile_data
        isolation_center_list = (isolation_center_list - knownMass) #/ (self.isolation_width / 2.0)
        intensity_list = intensity_list / np.max(intensity_list)
        rising_data_bool_list = isolation_center_list <= 0
        m_r,b_r = np.polyfit(isolation_center_list[rising_data_bool_list],intensity_list[rising_data_bool_list],1)
        m_f, b_f = np.polyfit(isolation_center_list[~rising_data_bool_list],
                                        intensity_list[~rising_data_bool_list], 1)
        #x_intersection = ((b2/m1)-(b1/m1))/(1-(m2/m1))
        x_intersection = ((b_f/m_r)-(b_r/m_r))/(1-(m_f/m_r))
        left_x_intercept = -b_r/m_r + (self.isolation_width / 2.0)
        right_x_intercept = -b_f/m_f - (self.isolation_width / 2.0)
        slope_ratio = -m_f/m_r
        return (x_intersection,left_x_intercept,right_x_intercept,slope_ratio)

    def regression(self, correction_factors=[0,0,1]):
        data_points = len(self.profile_data[0])
        if data_points <= 1:
            warnings.warn("Triagnle profile of mass "+str(self.observed_mass)+" and width "+str(self.isolation_width)+" only has "+data_points+"data points: ;regression could not be preformed")
        elif data_points == 2:
            if np.isclose(506.2756,self.observed_mass,rtol=1e-4):
                print("here")
            isolation_center_list, intensity_list = self.profile_data
            left_Area, right_Area = intensity_list
            left_Center, right_Center = isolation_center_list
            width = self.isolation_width
            slope = -(left_Area + right_Area) / (right_Center - width - left_Center)
            y_intersection = (slope * right_Center + left_Area - slope * left_Center + right_Area) / 2
            x_intersection = ((y_intersection - left_Area) / slope) + left_Center
            return (x_intersection, y_intersection)
        elif data_points <= 3:
            #TODO
            x=1
        elif data_points <= 5:
            #TODO 
            x=1
        else:
            x = 1
        return (0,0)

# Holds a set of spectra from an experiment and processes profiles across them.
class ExperimentSpectra:
    def __init__(self, window_step: float):
        self.window_step: float = window_step
        self.spectra: List[Spectrum] = []

    def addSpectrum(self, spectrum: Spectrum):
        spectrum.filter_data()
        self.spectra.append(spectrum)

    def extract_profiles_from_experiment(self, min_data_points: int = 5,
                                         rel_mass_tol: float = 7.5e-6) -> \
            List[TriangleProfile]:
        """
        Extracts triangle intensity profiles for masses observed across spectra.
        """
        triangle_profiles: List[TriangleProfile] = []
        spectra_isolation_centers: np.nparray = np.array([], dtype='f')
        previously_observed_masses: np.nparray = np.array([], dtype='f')
        previously_observed_listof_intensitiesby_mass: List[List[float]] = []
        for spectrum in self.spectra:
            spectra_isolation_centers = np.append(spectra_isolation_centers, spectrum.isolation_center)
            current_spectrum_masses, current_spectrum_intensities = spectrum.ms_data
            if len(previously_observed_masses) == 0:
                previously_observed_masses = current_spectrum_masses
                previously_observed_listof_intensitiesby_mass = current_spectrum_intensities.reshape(-1, 1).tolist()
            else:
                accounted_current_mass_bool_list = np.full(len(current_spectrum_masses), False)
                continued_observed_masses: np.nparray = np.array([], dtype='f')
                continued_observed_listof_intensitiesby_mass: List[List[float]] = []
                for previously_observed_mass, previously_observed_intensities in zip(previously_observed_masses,
                                                                                 previously_observed_listof_intensitiesby_mass):
                    bool_mached_mass_list = np.isclose(previously_observed_mass, current_spectrum_masses, rtol=rel_mass_tol)
                    if any(bool_mached_mass_list):
                        accounted_current_mass_bool_list = np.logical_or(accounted_current_mass_bool_list, bool_mached_mass_list)
                        average_mass = np.average(current_spectrum_masses[bool_mached_mass_list])
                        weight = len(previously_observed_intensities)
                        weight_list = [weight / (weight + 1), 1 / (weight + 1)]
                        weighted_average_mass = np.average([previously_observed_mass, average_mass], weights=weight_list)
                        continued_observed_masses = np.append(continued_observed_masses, weighted_average_mass)
                        summed_intensity = np.sum(current_spectrum_intensities[bool_mached_mass_list])
                        intensitiy_list = np.append(previously_observed_intensities, summed_intensity)
                        continued_observed_listof_intensitiesby_mass.append(intensitiy_list.tolist())
                    elif len(previously_observed_intensities) >= min_data_points:
                        triangle_profile = TriangleProfile(previously_observed_mass, spectrum.isolation_width)
                        observed_isolation_centers = spectra_isolation_centers[-len(previously_observed_intensities) - 1:-1]
                        triangle_profile.add_pofile_data(observed_isolation_centers, previously_observed_intensities)
                        triangle_profiles.append(triangle_profile)
                unaccounted_current_mass_bool_list = ~accounted_current_mass_bool_list
                remaining_spectrum_intensities = current_spectrum_intensities[unaccounted_current_mass_bool_list].reshape(-1,
                                                                                                                  1).tolist()
                for intensity in remaining_spectrum_intensities: continued_observed_listof_intensitiesby_mass.append(
                    intensity)
                continued_observed_masses = np.append(continued_observed_masses,
                                                    current_spectrum_masses[unaccounted_current_mass_bool_list])
                previously_observed_masses = continued_observed_masses
                previously_observed_listof_intensitiesby_mass = continued_observed_listof_intensitiesby_mass
        for previously_observed_mass, previously_observed_intensities in zip(previously_observed_masses,
                                                                         previously_observed_listof_intensitiesby_mass):
            if len(previously_observed_intensities) >= min_data_points:
                triangle_profile = TriangleProfile(previously_observed_mass, isolation_width=spectrum.isolation_width)
                observed_isolation_centers = spectra_isolation_centers[-len(previously_observed_intensities):]
                triangle_profile.add_pofile_data(observed_isolation_centers, previously_observed_intensities)
                triangle_profiles.append(triangle_profile)
        return triangle_profiles