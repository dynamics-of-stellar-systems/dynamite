#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:14:17 2020

@author: maindl
"""

import yaml

class ConfigurationReaderYaml(object):
    """
    Reads the configuration file
    """
    def __init__(self, filename=None): # instantiate the objects here. instead of the dict, self.system will be a System object 
        if filename:
            with open(filename, 'r') as f:
                par = yaml.safe_load(f)
        else:
            raise FileNotFoundError(filename)
        self.params=par
        # for p in par:
        #     self.params.update(p)

