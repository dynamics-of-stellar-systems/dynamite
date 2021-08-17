#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Use the following if you want to address the function, classes, etc.
# contained in dynamite/data_prep/data_prep_test.py as
# data_prep_test.yourfunction(...), data_prep_test.YourClass(...), etc.
from dynamite.data_prep import data_prep_test

data_prep_test.data_prep_function_test()

t_none = data_prep_test.TestClass()
t_text = data_prep_test.TestClass('Vader')
t_none.printout()
t_text.printout()

# Use the following if you want to address the function, classes, etc.
# contained in dynamite/data_prep/data_prep_test.py as
# yourfunction(...), YourClass(...), etc.
# You can also - carefully - use:
# from dynamite.data_prep.data_prep_test import *
from dynamite.data_prep.data_prep_test import data_prep_function_test

data_prep_function_test('Grogu')

