#import numpy as np
import autograd.numpy as np
import math as m

def Rxyz(rotations):
    theta_x, theta_y, theta_z = rotations
    return np.dot(Rx(theta_x), np.dot(Ry(theta_y), Rz(theta_z)))

def Rx1D(theta):
    if theta % 360 == 0:
        sigma = 1
    else:
        sigma = -1
    return np.array([sigma])

def Rx2D(theta):
    """ Rotation around x axis
    :param theta angle in degrees"""
    theta = theta * np.pi / 180
    return np.array([[m.cos(theta), -m.sin(theta)],
                     [m.sin(theta), m.cos(theta)]])

def Rx(theta):
    """ Rotation around x axis
    :param theta angle in degrees"""
    theta = theta*np.pi/180
    return np.array([[1, 0, 0],
                      [0, m.cos(theta), -m.sin(theta)],
                      [0, m.sin(theta), m.cos(theta)]])

def Ry(theta):
    """ Rotation around y axis
    :param theta angle in degrees"""
    theta = theta * np.pi / 180
    return np.array([[m.cos(theta), 0, m.sin(theta)],
                      [0, 1, 0],
                      [-m.sin(theta), 0, m.cos(theta)]])

def Rz(theta):
    """ Rotation around z axis
    :param theta angle in degrees"""
    theta = theta * np.pi / 180
    return np.array([[m.cos(theta), -m.sin(theta), 0],
                      [m.sin(theta), m.cos(theta), 0],
                      [0, 0, 1]])