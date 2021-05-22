#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging

def data_prep_function_test(argument='none specified'):
    if argument == 'hello':
        print('Hello, world!')
        logging.warning('"Hello, world!" printed.')
    else:
        print(argument)
        logging.warning(f'"{argument}" printed.')

class TestClass():
    def __init__(self, text=None):
        self.text = text if text else 'Selftext'
    def printout(self):
        print(self.text)
        logging.warning(f'{__name__}.{__class__.__name__}: "{self.text}" '
                        'printed.')
