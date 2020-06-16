#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:18:18 2020

@author: maindl
"""

import numpy as np
import json
import yaml

class Parameter(object):

    def __init__(self,
                 name=None,
                 lo=None,
                 hi=None,
                 step=None,
                 fixed=False,
                 value=None,
                 minstep=None,
                 text=None,
                 sformat="%g"):
        self.name = name
        self.lo = lo
        self.hi = hi
        self.step = step
        self.fixed = fixed
        self.value = value
        self.minstep = minstep
        self.text = text
        self.sformat = sformat

    # def validate(self, p):
    #     """
    #     Validates p against parameter limits and returns:
    #         p clipped to [Parameter.lo,Parameter.hi]
    #         None if abs(Parameter.value - p) > Parameter.minstep and Parameter.fixed == False

    #     Parameters
    #     ----------
    #     p : float or int

    #     Returns
    #     -------
    #     p clipped to Parameter range (float) or None

    #     """
    #     if abs(self.value - p) > self.minstep and not self.fixed:
    #         return None
    #     return np.clip(p,self.lo,self.up)


class ConfigurationReaderJson(object):
    """
    Reads the configuration file
    """
    def __init__(self, filename=None):
        if filename:
            with open(filename, 'r') as f:
                par = json.load(f)
        else:
            raise FileNotFoundError(filename)
        self.params={}
        for p in par:
            self.params.update(p)

class ConfigurationReaderYaml(object):
    """
    Reads the configuration file
    """
    def __init__(self, filename=None):
        if filename:
            with open(filename, 'r') as f:
                par = yaml.safe_load(f)
        else:
            raise FileNotFoundError(filename)
        self.params={}
        for p in par:
            self.params.update(p)
        