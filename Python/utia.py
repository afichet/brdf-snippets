#!/usr/bin/env python3

"""
:mod: `utia` -- BTF-UTIA BRDF Python support
============================================

.. module:: utia 
    :synopsis: This module implements the support for UTIA BRDF material
               http://btf.utia.cas.cz/

.. moduleauthor:: Alban Fichet <alban.fichet@inria.fr>
"""

import struct
import math
from enum import Enum, unique

class Utia:
    @unique
    class FileType(Enum):
        binary = 1
        exr = 2
        png = 3

    @unique
    class ColorFormat(Enum):
        XYZ = 1
        sRGB = 2
        RGB = 3
        
    def __init__(self, utia_file,
                 nti = 6, ntv = 6, npi = 48, npv = 48,
                 file_type = FileType.binary,
                 color_format = ColorFormat.sRGB):
        """
        Initialize and load an UTIA BRDF file in binary format.
        Default parameters are suitable for the 150 materials in 
        'BRDF Database'. For 'BRDF Dense', these should be set as follow:
        - nti = 44
        - ntv = 44
        - npi = 180
        - npo = 180
        
        :param utia_file: The path of the file to load
        :param nti: number of theta_i values captured
        :param ntv: number of theta_v values captured
        :param npi: number of phi_i values captured
        :param npv: number of phi_v values captured
        :param file_type: Specifies the type of file to be read
        :param color_format: Specified the format of color used in the input 
        file
        """
        self.nti = nti
        self.ntv = ntv
        self.npi = npi
        self.npv = npv
        
        with open(utia_file, 'rb') as f:
            data = f.read()
            self.brdf = struct.unpack(str(3*nti*ntv*npi*npv) + 'd', data)

        if color_format is self.ColorFormat.sRGB:
            self.brdf = [self.sRGB_to_RGB(sRGB) for sRGB in self.brdf]
    
    def eval_raw(self, theta_i, phi_i, theta_v, phi_v):
        """
        Lookup the BRDF value for given incoming and outcoming angles.

        :param theta_i: Incoming elevation angle in radians
        :param phi_i: Incoming azimuthal angle in radians
        :param theta_v: Outgoing elevation angle in radians
        :param phi_v: Outgoing azimuthal angle in radians
        :return: A list of 3 elements giving the BRDF value for R, G, B in
        linear RGB D65)
        """
        theta_i_p, phi_i_p = self.positive_rad(theta_i, phi_i)
        theta_v_p, phi_v_p = self.positive_rad(theta_v, phi_v)

        if theta_i_p > math.pi/2 or theta_v_p > math.pi/2:
            return [0]*3
        
        p_max = 2*math.pi
        t_max = math.pi/2
        
        idx_ti = min(self.nti - 1, round(self.nti * theta_i_p / t_max))
        idx_pi = min(self.npi - 1, round(self.npi * phi_i_p / p_max))
        idx_tv = min(self.ntv - 1, round(self.ntv * theta_v_p / t_max))
        idx_pv = min(self.npv - 1, round(self.npv * phi_v_p / p_max))
        
        return self.__eval_idx(idx_ti, idx_pi, idx_tv, idx_pv)

    def eval_interpolated(self, theta_i, phi_i, theta_v, phi_v):
        """
        Lookup the BRDF value for given incoming and outcoming angles and
        perform an interpolation over theta_i, phi_i, theta_v, phi_v.

        :param theta_i: Incoming elevation angle in radians
        :param phi_i: Incoming azimuthal angle in radians
        :param theta_v: Outgoing elevation angle in radians
        :param phi_v: Outgoing azimuthal angle in radians
        :return: A list of 3 elements giving the BRDF value for R, G, B in
        linear RGB D65)
        """
        theta_i_p, phi_i_p = self.positive_rad(theta_i, phi_i)
        theta_v_p, phi_v_p = self.positive_rad(theta_v, phi_v)

        if theta_i_p > math.pi/2 or theta_v_p > math.pi/2:
            return [0]*3
        
        t_max = math.pi/2
        p_max = 2*math.pi

        idx_ti_b = self.__theta_i_idx(theta_i_p)
        idx_pi_b = self.__phi_i_idx(phi_i_p)
        idx_tv_b = self.__theta_v_idx(theta_v_p)
        idx_pv_b = self.__phi_v_idx(phi_v_p)

        # Calculate the indexes for interpolation
        idx_ti_b = idx_ti_b if idx_ti_b < self.nti - 1 else self.nti - 2
        idx_tv_b = idx_tv_b if idx_tv_b < self.ntv - 1 else self.ntv - 2
    
        idx_ti = [idx_ti_b, idx_ti_b + 1]
        idx_pi = [idx_pi_b, idx_pi_b + 1]
        idx_tv = [idx_tv_b, idx_tv_b + 1]
        idx_pv = [idx_pv_b, idx_pv_b + 1]

        # Calculate the weights
        weight_ti = [abs(x/self.nti - theta_i_p/t_max) for x in idx_ti]
        weight_pi = [abs(x/self.npi - phi_i_p/p_max) for x in idx_pi]
        weight_tv = [abs(x/self.ntv - theta_v_p/t_max) for x in idx_tv]
        weight_pv = [abs(x/self.npv - phi_v_p/p_max) for x in idx_pv]
        
        # Normalize the weigths
        weight_ti = [1 - x / sum(weight_ti) for x in weight_ti]
        weight_pi = [1 - x / sum(weight_pi) for x in weight_pi]
        weight_tv = [1 - x / sum(weight_tv) for x in weight_tv]
        weight_pv = [1 - x / sum(weight_pv) for x in weight_pv]

        idx_pi[1] = idx_pi[1] if idx_pi[1] < self.npi else 0
        idx_pv[1] = idx_pv[1] if idx_pv[1] < self.npv else 0
        
        ret_val = [0]*3
        
        for iti, wti in zip(idx_ti, weight_ti):
            for ipi, wpi in zip(idx_pi, weight_pi):
                for itv, wtv in zip(idx_tv, weight_tv):
                    for ipv, wpv in zip(idx_pv, weight_pv):
                        ret_val = [r + x * wti * wpi * wtv * wpv
                                   for r, x in
                                   zip(ret_val,
                                       self.__eval_idx(iti, ipi, itv, ipv))]

        return ret_val
    
    def __eval_idx(self, iti, ipi, itv, ipv):
        """
        Lookup a BRDF value knowing the index values

        :param iti: theta_i index
        :param ipi: phi_i index
        :param itv: theta_v index
        :param ipv: phi_v index
        :return: A list of 3 elements giving the BRDF value for R, G, B in
        linear RGB D65)
        """
        nc = self.npv * self.ntv
        nr = self.npi * self.nti

        idx = nc * (self.npi * iti + ipi) + self.npv * itv + ipv
        
        return [self.brdf[c*nc*nr + idx] for c in range(0,3)]

    def __theta_i_idx(self, theta_i):
        return max(0, min(self.nti-1, math.floor(self.nti * theta_i / (math.pi/2))))

    def __phi_i_idx(self, phi_i):
        return max(0, min(self.npi-1, math.floor(self.npi * phi_i / (2*math.pi))))
    
    def __theta_v_idx(self, theta_v):
        return max(0, min(self.ntv-1, math.floor(self.ntv * theta_v / (math.pi/2))))

    def __phi_v_idx(self, phi_v):
        return max(0, min(self.npv-1, math.floor(self.npv * phi_v / (2*math.pi))))
    
    
    def positive_rad(self, theta, phi):
        """
        Method to ensure a given set of angles in radians is returned as a positive 
        value.

        :param theta: The elevation angle in radian
        :param phi: The azimuthal angle in radian
        :return: The (theta, phi) angle set given in a range of [0:2*pi]
        """ 
        theta_bound = math.fmod(theta, 2*math.pi)
        phi_bound = phi
        
        if theta_bound < 0:
            theta_bound = abs(theta_bound)
            phi_bound += math.pi

        phi_bound = math.fmod(phi_bound, 2*math.pi)
        while phi_bound < 0:
            phi_bound += 2 * math.pi

        return (theta_bound, phi_bound)

    def sRGB_to_RGB(self, sRGB):
        """
        Convert sRGB color value to RGB
        
        :param sRGB: sRGB value
        :return: A linear RGB value corresponding to the sRGB input value
        """
        a = 0.055
        return (sRGB / 12.92 if sRGB < 0.04045
                else math.pow((sRGB + a) / (1 + a), 2.4))