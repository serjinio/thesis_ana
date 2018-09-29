#!/bin/python
# returns Ay for inverse kinematics p - 4He elastic scattering 
# from direct kinematics data by Moss
# 
# Given: Ay(theta1_dir_kin)
#
# Solution:
#  Build theta1(theta); theta(theta2)
#  Ay(theta2) = Ay(theta1(theta(theta2)))
#
# theta_inv_kin(theta2) - inv. kin. to CMS
# theta_dir_kin(theta1) - dir. kin. to CMS
# theta1_dir_kin(theta) - from inv. to dir. kin transform
#
# Ay(theta2) = Ay(theta1_dir_kin(theta_inv_kin(theta2_inv_kin)))


import math
import numpy as np


def lab_theta1_dir_kin(cms_theta):
    """Returns lab theta1 given CMS angle"""
    return np.arctan(np.sin(cms_theta) /
                     (np.cos(cms_theta) + 1./4.))

def cms_theta_inv_kin(lab_theta2):
    """Returns CMS theta angle given theta2 (residual particle) 
    lab angle"""
    return np.pi - 2 * lab_theta2

_cms_theta_angles = np.radians(np.arange(0,180))
_lab_theta1_angles = lab_theta1_dir_kin(_cms_theta_angles);

def dir_kin_theta1_from_theta(theta):
    """Returns theta1 (scattered particle) angle for dir. kin. case
    given CM theta angle."""
    return np.interp(theta, _cms_theta_angles, _lab_theta1_angles)

def moss_lab_cs_to_cms(cms_theta, lab_cs_value):
    """Return cross section in CM frame by the given CM theta angle
    (for direct kinematics case) and lab CS value."""
    m1, m2 = 1., 4.
    num = 1. + m1/m2*math.cos(cms_theta);
    den = (1. + m1**2/m2**2 + 2 * m1/m2*math.cos(cms_theta))**(3./2.)
    cms_cs_value = lab_cs_value * (num / den)
    return cms_cs_value

raw_moss_ay = [0.39, 0.208, -0.013, -0.184, -0.304, -0.455,
               -0.597, -0.68, -0.783, -0.927, -1.029, -0.957,
               -0.703, -0.320,-0.042,0.025, 0.076, 0.118, 0.065]
raw_moss_cs = [6.23, 4.81, 3.93, 3.08, 2.38, 1.83, 1.44, 1.08,
               0.782, 0.442, 0.249, 0.156, 0.102, 0.0791, 0.0559,
               0.0473, 0.0431, 0.0348, 0.0267]
raw_moss_theta1 = np.radians([29., 31., 33., 35., 37., 39.,
                              41., 43., 45., 49., 53., 57.,
                              61., 65., 73., 77., 80., 85., 90.])

assert(len(raw_moss_ay) == len(raw_moss_cs) == len(raw_moss_theta1))


# measurement range for S13 experiment in CMS angle
s13_theta2 = np.radians(np.arange(55, 71))
# CMS angle from inv. kin.
s13_theta = cms_theta_inv_kin(s13_theta2);
# theta1 angle for dir. kin. from CMS
s13_theta1 = dir_kin_theta1_from_theta(s13_theta);

moss_cs_lab = np.interp(s13_theta1, raw_moss_theta1, raw_moss_cs,
                        left=-9999, right=-9999)
moss_cs_cms = [moss_lab_cs_to_cms(cms_theta, cs_val) for cms_theta, cs_val
               in zip(s13_theta, moss_cs_lab)]


def print_ay():
    s13_ay = np.interp(s13_theta1, raw_moss_theta1, raw_moss_ay,
                       left=-9999, right=-9999);
    s13_cs = np.interp(s13_theta1, raw_moss_theta1, raw_moss_cs,
                       left=-9999, right=-9999);

    print "S13 angles:\n"
    print "{theta2},{theta1},{theta},{ay},{cs}".format(
            theta2="A_inv", theta1="A_dir", theta="A_cms",
            ay="Ay", cs="CS")
    for theta1, theta2, ay, cs, theta in zip(s13_theta1, s13_theta2, s13_ay, s13_cs, s13_theta):
        print "{theta2:.2f},{theta1:.2f},{theta:.2f},{ay:.2f},{cs:.2f}".format(
            theta2=np.degrees(theta2), theta1=np.degrees(theta1), theta=np.degrees(theta),
            ay=ay, cs=cs)

def print_angles():
    print "{s13_th2},{s13_th1},{moss_theta1},{cm}".format(
        s13_th2="S13 A2", s13_th1="S13 A1", moss_theta1="Moss A1", cm="theta CM")
    for s13_th2, s13_th1, moss_th1, cm_th in \
        zip(s13_theta2, s13_theta1, raw_moss_theta1, s13_theta):
        print "{s13_th2:.2f},{s13_th1:.2f},{moss_th1:.2f},{cm:.2f}".format(
            s13_th2=np.degrees(s13_th2), s13_th1=np.degrees(s13_th1),
            moss_th1=np.degrees(moss_th1), cm=np.degrees(cm_th))

def print_moss_cms_cs():
    thetas, cs_values = [], []
    for theta, cs_cms in zip(s13_theta, moss_cs_cms):
        thetas.append(theta);
        cs_values.append(cs_cms);
    thetas.reverse();
    cs_values.reverse();

    print "{theta},{cs_c}".format(theta="A_cms", cs_c="CS (CM)")
    for theta, cs_cms in zip(thetas, cs_values):
        print "{theta:.2f},{cs_c:.2f}".format(
            theta=np.degrees(theta), cs_c=cs_cms)

def print_moss_cms_lab():
    print "{lab_theta1},{cs}".format(lab_theta1="A1_lab", cs="CS (lab)")
    for theta, cs_lab in zip(s13_theta1, moss_cs_lab):
        print "{theta:.2f},{cs_l:.2f}".format(
            theta=np.degrees(theta), cs_l=cs_lab)


if __name__ == "__main__":
    # print_angles()
    print_moss_cms_lab()
    print_moss_cms_cs()
    print_ay()
